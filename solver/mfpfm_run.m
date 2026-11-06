function out = mfpfm_run(p)
% MFPFM_RUN  Multi-flake phase-field model: phase field phi, adatom concentration c and
% orientation theta on a periodic square grid. The equations, discretisation and
% Crank-Nicolson step for c are those of KWC_EulerCN_standalone_v2.m; this version solves
% the implicit c step with FFTs, which gives the same result to round-off and runs faster.
%
% Required fields of p:
%   n, g        symmetry order and anisotropy strength, kappa = k^2 (1 + g cos n psi)
%   lambda      coupling constant
%   conc_eq     equilibrium adatom concentration
%   tau_s       adatom lifetime (desorption)
%   tau_phi     relaxation time of phi
%   D, F        adatom diffusivity and deposition flux
%   tend        end time
% Optional fields (defaults in brackets):
%   L [256], dx [1], dt [1e-3], k [sqrt(2)], lambda_theta [0.1], theta_bg [3]
%   wrap ['2pi']        orientation differences modulo 2*pi; 'sym' = modulo 2*pi/n,
%                       'sym3' = modulo 2*pi/3, 'raw' = no reduction (original script)
%   theta_weight ['eta'] weight of the orientation term in the phi equation, (1+phi)/2;
%                       'phi' = weight phi (original script)
%   seeds [none]        rows [t, x, y, theta, radius] of scheduled seeds (x, y from the centre)
%   seed_on_bare [true] skip a scheduled seed whose centre is already covered
%   rgrain [false]      random nucleation, probability 1/1000 per step until nuc_tmax [15]
%   rseed [2.5]         radius of random nuclei
%   reject_covered [true] place a random nucleus only on bare substrate (phi < 0)
%   rng_seed [0]        seed of the random-number generator
%   snap_times []       times at which phi, c and theta are stored in out.snaps
%   snap_every [0]      if > 0, phi is also stored every snap_every time units in out.phis
%   stop_area [0]       if > 0, stop once the area with phi > 0 reaches stop_area and t >= t_min
%   diag_every [1000]   steps between entries of out.At = [t, area, max phi]
%   survive_dt [5]      a random nucleus counts as survived if phi > 0 at its site survive_dt later
%   progfile ['']       text file updated every 10000 steps with t, tend, area and wall time
%   verbose [false]
% Output: final fields out.phi, out.conc, out.theta (single), snapshots out.snaps and
% out.phis, area history out.At, nucleation log out.NUC = [t, i, j, theta, c, placed,
% survived], seed log out.SEEDLOG, the parameters out.p, the largest mass-balance error
% out.maxres and the wall time out.wall.
d = struct('L',256,'dx',1,'dt',1e-3,'k',sqrt(2),'lambda_theta',0.1,'theta_bg',3, ...
    'wrap','2pi','theta_weight','eta','rseed',2.5,'reject_covered',true,'rgrain',false, ...
    'rng_seed',0,'nuc_tmax',15,'snap_times',[],'diag_every',1000,'seeds',zeros(0,5), ...
    'survive_dt',5,'verbose',false,'snap_every',0,'stop_area',0,'t_min',0,'seed_on_bare',true,'progfile','');
fn = fieldnames(d); for q = 1:numel(fn), if ~isfield(p,fn{q}), p.(fn{q}) = d.(fn{q}); end, end

dx = p.dx; dt = p.dt; N = round(p.L/dx); nstep = round(p.tend/dt);
k = p.k; g = p.g; n = p.n; lambda = p.lambda; conc_eq = p.conc_eq;
tau_s = p.tau_s; tau_phi = p.tau_phi; D = p.D; F = p.F; lt = p.lambda_theta;
thetaSmallValue = 1e-3; gamma_param = 1e3; beta = 1e5; mu = 1e3;
tauT = 0.01; epsx = 0.008; s = 0.010; eps_gamma = 1/(gamma_param*10);

switch p.wrap
    case 'raw', wrapd = @(x) x;                          wrapv = @(x) mod(x + pi, 2*pi) - pi;
    case '2pi', wrapd = @(x) mod(x + pi, 2*pi) - pi;     wrapv = wrapd;
    case 'sym', per = 2*pi/n; wrapd = @(x) mod(x + per/2, per) - per/2; wrapv = wrapd;
    case 'sym3', per = 2*pi/3; wrapd = @(x) mod(x + per/2, per) - per/2; wrapv = wrapd;
end

rng(p.rng_seed);
phi = -ones(N); theta = p.theta_bg*ones(N);
conc = ones(N)*F*tau_s; conc_old = conc;
[J, I] = meshgrid(1:N, 1:N);
Xc = (I - N/2)*dx; Yc = (J - N/2)*dx;
kk = 0:N-1; l1 = (2*cos(2*pi*kk/N) - 2)/dx^2; Lhat = l1' + l1;
Ahat = 1/dt + 0.5/tau_s - 0.5*D*Lhat;
lap = @(f) (circshift(f,[0,1])+circshift(f,[0,-1])+circshift(f,[1,0])+circshift(f,[-1,0]) - 4*f)/dx^2;
Dc  = @(f) (circshift(f,[0,-1]) - circshift(f,[0,1]))/(2*dx);
Dr  = @(f) (circshift(f,[-1,0]) - circshift(f,[1,0]))/(2*dx);

NUC = zeros(0,7);   % [t, i, j, theta, c_at_site, placed(1)/rejected(0), survived(1/0/NaN)]
SEEDLOG = zeros(0,3); % [t, seed index, placed]
snaps = struct('t',{},'phi',{},'conc',{},'theta',{});
phis = struct('t',{},'phi',{});
At = zeros(0,3);
maxres = 0; srccum = 0;
nse = 0; if p.snap_every > 0, nse = round(p.snap_every/dt); end
tic;
for istep = 1:nstep
    t_prev = (istep-1)*dt; if istep==1, t_prev = -inf; end
    t_now = istep*dt;
    ia = find(p.seeds(:,1) > t_prev & p.seeds(:,1) <= t_now);
    for a = ia'
        m = (Xc - p.seeds(a,2)).^2 + (Yc - p.seeds(a,3)).^2 < p.seeds(a,5)^2;
        [~, ic] = min(((Xc(:) - p.seeds(a,2)).^2 + (Yc(:) - p.seeds(a,3)).^2));
        ok = ~(p.seed_on_bare && phi(ic) >= 0);
        if ok, phi(m) = 1; theta(m) = p.seeds(a,4); end
        SEEDLOG(end+1,:) = [t_now, a, ok]; %#ok<AGROW>
    end
    P_N = randi([0, 999]);
    if p.rgrain && P_N < 1 && t_now < p.nuc_tmax
        Nxi = randi(N); Nyi = randi(N); rt = deg2rad(randi(360));
        ok = ~(p.reject_covered && phi(Nxi,Nyi) >= 0);
        if ok
            if p.rseed <= 1
                phi(Nxi,Nyi) = 1; theta(Nxi,Nyi) = rt;
            else
                di = mod(I - Nxi + N/2, N) - N/2; dj = mod(J - Nyi + N/2, N) - N/2;
                m = (di.^2 + dj.^2)*dx^2 < p.rseed^2; phi(m) = 1; theta(m) = rt;
            end
        end
        NUC(end+1,:) = [t_now, Nxi, Nyi, rt, conc(Nxi,Nyi), ok, NaN]; %#ok<AGROW>
    end
    if ~isempty(NUC)
        chk = find(isnan(NUC(:,7)) & NUC(:,6)==1 & t_now >= NUC(:,1) + p.survive_dt);
        for c = chk', NUC(c,7) = phi(NUC(c,2), NUC(c,3)) > 0; end
    end
    phiold = phi; concold = conc; th = theta;

    % theta equation
    tip1 = circshift(th,[0,-1]); tim1 = circshift(th,[0,1]);
    tjp1 = circshift(th,[-1,0]); tjm1 = circshift(th,[1,0]);
    phi_ip1 = circshift(phiold,[0,-1]); phi_jp1 = circshift(phiold,[-1,0]);
    gx = wrapd(tip1 - tim1)/(2*dx); gy = wrapd(tjp1 - tjm1)/(2*dx);
    gmag = max(sqrt(gx.^2 + gy.^2), thetaSmallValue);
    expo = min(epsx*beta*gmag, 100);
    pFunc = 1 + exp(-expo).*(mu/epsx - 1);
    gphi = (phiold + 1)/2;
    phi_mod_sq = (gphi + (gphi < thetaSmallValue)*thetaSmallValue).^2;
    dxp = wrapd(tip1 - th); dyp = wrapd(tjp1 - th);
    gfx = max(abs(dxp)/dx, eps_gamma);
    IGx = (gfx > 1/gamma_param).*(1./gfx - gamma_param) + gamma_param;
    Dip = ((gphi + (phi_ip1+1)/2)/2).^2 .* (s*IGx + epsx^2);
    gfy = max(abs(dyp)/dx, eps_gamma);
    IGy = (gfy > 1/gamma_param).*(1./gfy - gamma_param) + gamma_param;
    Djp = ((gphi + (phi_jp1+1)/2)/2).^2 .* (s*IGy + epsx^2);
    fx = Dip.*dxp; fy = Djp.*dyp;
    div = (fx - circshift(fx,[0,1]) + fy - circshift(fy,[1,0]))/dx^2;
    theta = wrapv(th + dt*div./(tauT*phi_mod_sq.*pFunc));

    % phase-field equation
    ty = wrapd(circshift(theta,[0,-1]) - circshift(theta,[0,1]))/(2*dx);
    tx = wrapd(circshift(theta,[-1,0]) - circshift(theta,[1,0]))/(2*dx);
    gtm = sqrt(tx.^2 + ty.^2);
    phidy = Dc(phi); phidx = Dr(phi);
    psi = theta + atan2(phidy, phidx);
    kap = k^2*(1 + g*cos(n*psi)); kapd = -k^2*g*n*sin(n*psi); kap2 = kap.^2;
    term1 = Dc(kap.*kapd.*phidx); term2 = Dr(-kap.*kapd.*phidy);
    gk2 = Dr(kap2).*phidx + Dc(kap2).*phidy;
    dw = sin(pi*phiold) + lambda*(concold - conc_eq).*(1 + cos(pi*phiold));
    if strcmp(p.theta_weight,'eta'), wfac = (1+phiold)/2; else, wfac = phiold; end
    tterm = -2*s*wfac.*gtm - epsx^2*wfac.*gtm.^2;
    phi = phi + (dt/tau_phi)*(term1 + term2 + kap2.*lap(phi) + gk2 + dw + lt*tterm);

    % concentration
    ct = 1.5*concold - 0.5*conc_old;
    rhs = (1/dt - 0.5/tau_s)*concold + 0.5*D*lap(concold) + (1/tau_s)*phi.*ct ...
          + F*(1 - phi) - 0.5*(phi - phiold)/dt;
    cnew = real(ifft2(fft2(rhs)./Ahat));
    srcn = dt*sum(sum(F*(1-phi) - (concold+cnew)/(2*tau_s) + phi.*ct/tau_s))*dx^2;
    maxres = max(maxres, abs(sum(sum((cnew - concold) + (phi - phiold)/2))*dx^2 - srcn));
    srccum = srccum + srcn;
    conc_old = concold; conc = cnew;

    for ss = 1:numel(p.snap_times)
        if abs(t_now - p.snap_times(ss)) < dt/2
            snaps(end+1) = struct('t',t_now,'phi',single(phi),'conc',single(conc),'theta',single(theta)); %#ok<AGROW>
        end
    end
    if nse > 0 && mod(istep, nse) == 0
        phis(end+1) = struct('t',t_now,'phi',single(phi)); %#ok<AGROW>
    end
    if mod(istep, p.diag_every) == 0
        A = sum(phi(:) > 0)*dx^2; At(end+1,:) = [t_now, A, max(phi(:))]; %#ok<AGROW>
        if p.stop_area > 0 && A >= p.stop_area && t_now >= p.t_min, break; end
    end
    if p.verbose && mod(istep,10000)==0, fprintf('t=%g\n', t_now); end
    if ~isempty(p.progfile) && mod(istep,10000)==0
        fp = fopen(p.progfile,'w'); fprintf(fp,'%g %g %g %.1f\n', t_now, p.tend, sum(phi(:)>0)*dx^2, toc); fclose(fp);
    end
end
out.wall = toc; out.p = p; out.NUC = NUC; out.SEEDLOG = SEEDLOG; out.maxres = maxres;
out.snaps = snaps; out.phis = phis; out.At = At; out.tfinal = t_now;
out.phi = single(phi); out.conc = single(conc); out.theta = single(theta);
out.maxphi = max(phi(:)); out.minphi = min(phi(:));
end
