%% Phase-Field Modeling of Multi-Flake Growth and Grain-Boundary Formation in CVD-Grown 2D Materials
% KWC_EulerCN_standalone_v2_Fig6.m  -  KWC_EulerCN_standalone_v2.m set up for Fig. 6(a-c)
% (two flakes, 40 and 110 deg, nucleated together; for Fig. 6(d-f) swap the commented
% rows of nucleation_schedule)
%
% Changes relative to the original live script (all marked "% FIX" in the code):
%   1a  Orientation differences and gradients are taken the short way round the
%       circle (modulo 2*pi), as FiPy's ModularVariable does. Previously raw
%       differences were used, so neighbouring values on either side of +/-pi
%       produced a ~2*pi jump (spurious grain boundaries, and drift of the
%       orientation of flakes assigned near 180 deg).
%   1b  The orientation term of the phase-field equation (Eq. 1) is weighted by
%       eta = (1+phi)/2 instead of phi, so that it vanishes on the bare substrate,
%       as in FiPy where phi = 0 in the liquid.
%   1d  Seeds are discs of radius seed = 2.5 grid cells (periodic distance for
%       random nuclei) instead of a single cell.
%   1e  Random nuclei are placed only if the nucleation site is on bare
%       substrate (phi < 0). The random-number draws are unchanged.
%   Not adopted (stated as a limitation in the paper): reduction of theta modulo
%   2*pi/n. To test it, set sym_reduce = true below.
%
% Everything else (equations, parameters, discretisation, CN solver, plotting)
% is identical to the original.

clear;
clc;

hFig = figure('Color','white');
gif_filename = 'simulation_output.gif';
if exist(gif_filename, 'file')
    delete(gif_filename); % Delete old file so it doesn't append to a previous run
end

%% Switches
% Nucleation switch
%if both are false, nucleation takes place using the nucleus_kwc function
rgrain              = false;    % random nucleation
rsched              = true;     % scheduled nucleation for bowtie simulations (Fig. 6)

% Visualization switches
do_plot             = true;     % live visualisation
do_gif              = false;    % create gif
do_vtk              = false;    % create vtk
show_res            = false;    % show Step-to-step residual information

% FIX 1a: orientation differences on the circle
sym_reduce          = false;    % false = modulo 2*pi (adopted, as FiPy)
                                % true  = modulo 2*pi/n (Warren et al. 2003; not adopted)

%% Grid and Time
% Get initial wall time:
time0  = clock();

% Simulation cell parameters:
Nx   = 256;
Ny   = 256;
NxNy = Nx * Ny;
dx   = 1;
dy   = dx;

% Time integration parameters:
nstep  = 100000;   % Fig. 6: t = 100
nprint =    500;
dtime  = 1.e-3;

%% Phase field and Concentration equation parameters
% tau_phi dphi/dt = d/dy(kappa kappa' dphi/dx) - d/dx(kappa kappa' dphi/dy) + div(kappa^2 grad phi)
%                   + sin(pi phi) + lambda (c - c_eq)[1 + cos(pi phi)]
%                   + lambda_theta [-2 s eta|grad theta| - eps^2 eta |grad theta|^2],  eta = (1+phi)/2   (FIX 1b)
% dc/dt = D lap(c) + (F - c/tau_s)(1 - phi) - (1/2) dphi/dt
k        = sqrt(2);
lambda_theta = 0.1;
%============================
g        = 0.1;
n        = 3.0;        % symmetry order
lambda   = 18.0;
conc_eq  = 0.1;
tau_s    = 200.0;

%====================
tau_phi = 5.0;
D       = 2.0;
F       = 0.002;
%====================

%============================
% % Figure 5 Settings for 6-fold symmetry (Uncomment to use; also turn on rgrain switch)
%============================
% g = 0.04;
% n = 6.0; %symmetry
% lambda = 200.0;
% conc_eq = 0.03;
% tau_s = 18.0;
% lambda_theta = 0.1;
% %====================
% tau_phi = 100.0;    %compact fig. 5(a-c)
% % tau_phi = 5.0;      %dendrite fig. 5(d-f)
% D = 10.0;
% F = 0.007;
% %====================

% FIX 1a: shortest signed difference between two orientations
if sym_reduce
    per = 2*pi/n;
else
    per = 2*pi;
end
wrapd = @(x) mod(x + per/2, per) - per/2;

%% Theta equation parameters
% P(eps|grad theta|) tau_theta eta^2 dtheta/dt = div[ eta^2 (s/|grad theta| + eps^2) grad theta ]
% P(w) = 1 - exp(-beta w) + (mu/eps) exp(-beta w)
% Interested readers are referred to the implementation of this equation in FiPy:
% FiPy > examples > examples.phase.polyxtalCoupled
% Guyer, Jonathan E., Daniel Wheeler, and James A. Warren. "FiPy: Partial differential
% equations with Python." Computing in Science & Engineering 11.3 (2009): 6-15.
thetaSmallValue     = 1e-3;

gamma_param         = 1e3;
beta                = 1e5;
mu                  = 1e3;

thetaTransientCoeff = 0.01;

epsilonex           = 0.008;
s                   = 0.010;

%% Nuclei and Array Initialization
seed    = 2.5; % FIX 1d: seed radius in grid cells (original: 1, i.e. a single cell)
ngrain  = 1; % no of grains for pre determined positions and orientations,
             % you can change the values in nucleus_kwc

if rgrain == true
    ngrain = 0; % for random grain or bowtie mode so that the simulation starts
                % with no flakes
elseif rsched == true
    ngrain = 0;
    tau_phi = 5.0; D= 10.0; F= 0.002;
    g      = 0.3;
    % Time based nucleation scheduling to capture the effects of bowtie
    nucleation_schedule = [
    %time, Ni(y), Nj(x), theta_0
    %Fig 6(a-c)
    0, 5,  -5,  deg2rad(40);
    0, -5,  5,  deg2rad(110);
    %Fig 6(d-f)
    % 0, -36/2,  -30/2,  deg2rad(29);
    % 5, 30/2,  36/2,  deg2rad(105);
    ];
end

% Array initialization
[phi, theta] = nucleus_kwc(Nx, Ny, seed, ngrain);
conc         = ones(Nx, Ny) * F * tau_s;
conc_old     = conc;   % history for CN scheme for the C field

%--- Sparse periodic Laplacian:
laplacian_mat = build_laplacian(Nx, Ny, dx, dy);

fprintf('Building CN matrix and computing sparse LU factorisation...\n');

A_cn = (1/dtime + 0.5/tau_s) * speye(NxNy) - 0.5 * D * laplacian_mat;

[L_cn, U_cn, P_cn, Q_cn] = lu(A_cn);
fprintf('LU factorisation complete.\n\n');

%--- Stability information (informational only) ---
r_explicit = D * dtime / (dx^2);
fprintf('=== Stability information ===\n');
fprintf('  Explicit Fourier number r = D*dt/dx^2 = %.4f\n', r_explicit);
if r_explicit > 0.25
    fprintf('  WARNING: r > 0.25 - explicit conc. would be UNSTABLE.\n');
    fprintf('           CN removes this constraint.\n');
else
    fprintf('  r <= 0.25: explicit stable, but CN used for robustness.\n');
end
fprintf('\n');

%% TIME LOOP
for istep = 1:nstep

    % Random nucleation ---
    P_N = randi([0, 999]);
    if rgrain == true && P_N < 1 && istep < 15000
    % Seeds are allowed to grow only until step 15000
    % so that that the simulation space is not over-saturated
    % and the initilialized grains get the time to fully grow
        Nxi        = randi(Nx);
        Nyi        = randi(Ny);
        rand_theta = deg2rad(randi(360));
        if phi(Nxi, Nyi) < 0                       % FIX 1e: nucleate only on bare substrate
            for i = 1:Nx
                for j = 1:Ny
                    di = mod(i - Nxi + Nx/2, Nx) - Nx/2;   % FIX 1d: periodic distance
                    dj = mod(j - Nyi + Ny/2, Ny) - Ny/2;
                    if di^2 + dj^2 < seed^2                % FIX 1d: disc of radius seed
                        phi(i,j)   = 1.0;
                        theta(i,j) = rand_theta;
                    end
                end
            end
        end
    end

    % Scheduled nucleation for bowtie
    if rsched == true
        if istep == 1
            t_prev = -inf;              % include seeds scheduled at t = 0
        else
            t_prev = (istep - 1) * dtime;
        end
        t_now = istep * dtime;

        active_seeds = nucleation_schedule( ...
            nucleation_schedule(:,1) > t_prev & ...
            nucleation_schedule(:,1) <= t_now, :);

        ic = (Nx/2);               % domain centre, in cells
        jc = (Ny/2);

        for iseed = 1:size(active_seeds, 1)
            i_s    = active_seeds(iseed, 2);   % x offset from centre
            j_s    = active_seeds(iseed, 3);   % y offset from centre
            orient = active_seeds(iseed, 4);

            for i = 1:Nx
                for j = 1:Ny
                    i_centered = i - ic;
                    j_centered = j - jc;
                    if (i_centered - i_s)^2 + (j_centered - j_s)^2 < seed^2
                        phi(i,j)   = 1.0;
                        theta(i,j) = orient;
                    end
                end
            end
        end
    end
    % Store fields at time n
    phiold    = phi;
    concold   = conc;
    theta_old = theta;

    %% THETA EQUATION
    theta_ip1 = circshift(theta_old, [0, -1]);
    theta_im1 = circshift(theta_old, [0,  1]);
    theta_jp1 = circshift(theta_old, [-1, 0]);
    theta_jm1 = circshift(theta_old, [ 1, 0]);

    phi_ip1 = circshift(phiold, [0, -1]);
    phi_jp1 = circshift(phiold, [-1, 0]);

    % Calculate |grad theta| at cell centers for PFunc
    grad_theta_x       = wrapd(theta_ip1 - theta_im1) / (2*dx);   % FIX 1a
    grad_theta_y       = wrapd(theta_jp1 - theta_jm1) / (2*dy);   % FIX 1a
    grad_theta_mag_old = sqrt(grad_theta_x.^2 + grad_theta_y.^2);

    % Calculate P function
    grad_theta_mag_old = max(grad_theta_mag_old, thetaSmallValue);
    expo   = min(epsilonex * beta * grad_theta_mag_old, 100);
    pFunc  = 1.0 + exp(-expo) .* (mu / epsilonex - 1.0);

    % Calculate phi modulation term
    g_phi      = (phiold + 1.0) / 2.0;
    phi_mod_sq = (g_phi + (g_phi < thetaSmallValue) * thetaSmallValue).^2;

    % Face differences of theta, taken on the circle                % FIX 1a
    dth_ip = wrapd(theta_ip1 - theta_old);     % across face i+1/2
    dth_im = wrapd(theta_old - theta_im1);     % across face i-1/2
    dth_jp = wrapd(theta_jp1 - theta_old);     % across face j+1/2
    dth_jm = wrapd(theta_old - theta_jm1);     % across face j-1/2

    % Calculate Diffusion Coefficients D at the 4 cell faces
    eps_gamma       = 1.0 / (gamma_param * 10.0);

    % Faces i+1/2 (right face of each cell)
    grad_mag_face_x = max(abs(dth_ip) / dx, eps_gamma);              % FIX 1a
    IGamma_x        = (grad_mag_face_x > 1/gamma_param) .* (1./grad_mag_face_x - gamma_param) + gamma_param;
    g_phi_face_x    = ((phiold + 1.0)/2.0 + (phi_ip1 + 1.0)/2.0) / 2.0;
    D_ip1_2         = g_phi_face_x.^2 .* (s * IGamma_x + epsilonex^2);
    D_im1_2         = circshift(D_ip1_2, [0, 1]);

    % Faces j+1/2 (top face of each cell)
    grad_mag_face_y = max(abs(dth_jp) / dy, eps_gamma);              % FIX 1a
    IGamma_y        = (grad_mag_face_y > 1/gamma_param) .* (1./grad_mag_face_y - gamma_param) + gamma_param;
    g_phi_face_y    = ((phiold + 1.0)/2.0 + (phi_jp1 + 1.0)/2.0) / 2.0;
    D_jp1_2         = g_phi_face_y.^2 .* (s * IGamma_y + epsilonex^2);
    D_jm1_2         = circshift(D_jp1_2, [1, 0]);

    % Calculate the divergence term: div(D * grad(theta))
    flux_x     = (D_ip1_2 .* dth_ip - D_im1_2 .* dth_im) / dx^2;    % FIX 1a
    flux_y     = (D_jp1_2 .* dth_jp - D_jm1_2 .* dth_jm) / dy^2;    % FIX 1a
    divergence = flux_x + flux_y;

    % Calculate the time derivative d(theta)/dt
    d_theta_dt = (1.0 ./ (thetaTransientCoeff * phi_mod_sq .* pFunc)) .* divergence;

    % Time integration
    theta = theta_old + dtime * d_theta_dt;
    % Wrap theta back onto the circle
    theta = wrapd(theta);   % identical to mod(theta + pi, 2*pi) - pi when sym_reduce = false

    %% PHASE FIELD EQUATION

    % Calculate derivatives
    % FIX 1a: centred differences of theta taken on the circle
    % (replaces [theta_y, theta_x] = gradient_mat(theta, Nx, Ny, dx, dy))
    theta_y = wrapd(circshift(theta, [0, -1]) - circshift(theta, [0, 1])) / (2*dx);
    theta_x = wrapd(circshift(theta, [-1, 0]) - circshift(theta, [1, 0])) / (2*dy);
    gradthetamag       = sqrt(theta_x.^2 + theta_y.^2);
    gradthetamagsq     = gradthetamag.^2;

    % Calculate the laplacians and epsilon:
    phi2     = reshape(phi', NxNy, 1);
    lap_phi2 = laplacian_mat * phi2;
    lap_phi  = vec2matx(lap_phi2, Nx);

    % Gradients of phi
    [phidy, phidx] = gradient_mat(phi, Nx, Ny, dx, dy);

    % Calculate angle
    psi           = theta + atan2(phidy, phidx);

    % Kappa and its derivative
    kappa       = (k^2) * (1 + g * cos(n * psi));
    kappa_deriv = -k*k * g * n * sin(n * psi);
    kappa2 = kappa.^2;

    % Kobayashi first term
    dummyx = kappa .* kappa_deriv .* phidx;
    [term1, ~] = gradient_mat(dummyx, Nx, Ny, dx, dy);

    % Kobayashi second term
    dummyy = -kappa .* kappa_deriv .* phidy;
    [~, term2] = gradient_mat(dummyy, Nx, Ny, dx, dy);

    % Second divergence term
    [dkappa2_dy,dkappa2_dx]=gradient_mat(kappa2,Nx,Ny,dx,dy);
    grad_kappa2_dot_grad_phi = dkappa2_dx.*phidx + dkappa2_dy.*phidy;

    % Double well
    dw = sin(pi * phiold) + lambda * (concold - conc_eq) .* (1 + cos(pi * phiold));
    % Theta term                                                      % FIX 1b
    eta_old    = (phiold + 1.0) / 2.0;      % 0 on the substrate, 1 inside the flake
    theta_term = -2*s * eta_old .* gradthetamag - epsilonex^2 * eta_old .* gradthetamagsq;

    phi = phi + (dtime / tau_phi) * (term1 + term2 + kappa.^2 .* lap_phi + ...
                 grad_kappa2_dot_grad_phi + dw + lambda_theta * theta_term);

    %% CONCENTRATION (Crank Nicolson)
    c_tilde = 1.5 * concold - 0.5 * conc_old;

    % Laplacian of c^n (explicit half of diffusion)
    tempx     = reshape(concold', NxNy, 1);
    lap_c_n   = vec2matx(laplacian_mat * tempx, Nx);

    % Build RHS
    rhs_conc = (1/dtime - 0.5/tau_s) * concold ...      % time derivative + explicit -c/tau_s half
             + 0.5 * D * lap_c_n ...                    % explicit half of diffusion
             + (1/tau_s) * phi .* c_tilde ...           % nonlinear phi*c via AB2
             + F * (1 - phi) ...                        % deposition source (phi^{n+1})
             - 0.5 * (phi - phiold) / dtime;            % phase-field coupling

    % Solve A_cn * c^{n+1} = rhs
    rhs_vec      = reshape(rhs_conc', NxNy, 1);
    conc_new_vec = Q_cn * (U_cn \ (L_cn \ (P_cn * rhs_vec)));
    conc_new     = vec2matx(conc_new_vec, Nx);

    % Advance AB2 history BEFORE overwriting conc
    conc_old = concold;   % c^{n-1} <- c^n
    conc     = conc_new;  % c^n     <- c^{n+1}

    %% VISUALIZATION
    if mod(istep, nprint) == 0
        if show_res == true
        % Step-to-step residual
        err_phi  = sqrt(dx*dy) * norm(phi(:)  - phiold(:));
        err_conc = sqrt(dx*dy) * norm(conc(:) - conc_new_vec);   % uses c^{n+1} vs c^n
        fprintf('Step %d/%d | res_phi = %.4e | res_c = %.4e\n', ...
                istep, nstep, err_phi, err_conc);
        else
        fprintf('Step %d/%d\n', istep, nstep);
        end

        if do_plot
            set(gcf, 'Color', 'white');
            set(gcf, 'Visible', 'on');
            set(0, 'CurrentFigure', hFig);

            subplot(1,3,1);
            imagesc([0 Nx], [0 Ny], phi);
            clim([-1 1]);
            axis equal tight; colorbar;
            set(gca, 'Layer', 'top');
            title(['\phi at t = ', num2str(istep * dtime)]);
            xlabel('x'); ylabel('y');

            subplot(1,3,2);
            imagesc([0 Nx], [0 Ny], conc);
            axis equal tight; colorbar;
            set(gca, 'Layer', 'top');
            title(['\Delta C at t = ', num2str(istep * dtime)]);
            xlabel('x'); ylabel('y');

            subplot(1,3,3);
            imagesc([0 Nx], [0 Ny], theta);
            axis equal tight; colorbar;
            clim([-per/2 per/2]);
            set(gca, 'Layer', 'top');
            title(['\theta at t = ', num2str(istep * dtime)]);
            xlabel('x'); ylabel('y');

            sgtitle(['Step: ' num2str(istep) ' / ' num2str(nstep)]);
            drawnow;
        end
        % --- Append current figure to the GIF ---
        if do_gif
            exportgraphics(gcf, gif_filename, 'Append',true);
        end

        if do_vtk
            write_vtk_grid_values(Nx,Ny,dx,dy,istep,phi,conc,theta);
        end

    end

end % istep

compute_time = etime(clock(), time0);
fprintf('Compute Time: %10d s\n', compute_time);


%% LOCAL FUNCTIONS
function [phi, theta] = nucleus_kwc(Nx, Ny, seed, ngrain)
    format long;
    phi   = -1 * ones(Nx, Ny);
    theta =  3 * ones(Nx, Ny);

    orientations = [0, deg2rad(105), -pi/3, -2*pi/3, pi/2];
    Ni_centered  = [0,    30/3,  Nx/4,  -Nx/4, Nx/4];
    Nj_centered  = [0,    37/3, -Ny/4,   Ny/4, Ny/4];

    for igrain = 1:ngrain
        for i = 1:Nx
            for j = 1:Ny
                i_c = i - Nx/2;
                j_c = j - Ny/2;
                if (i_c - Ni_centered(igrain))^2 + (j_c - Nj_centered(igrain))^2 < seed^2   % FIX 1d
                    phi(i,j)   = 1.0;
                    theta(i,j) = orientations(igrain);
                end
            end
        end
    end
end


function [grad] = build_laplacian(nx, ny, dx, dy)
    format long;
    nxny   = nx * ny;
    r      = zeros(1, nx);
    r(1:2) = [2, -1];
    T      = toeplitz(r);
    E      = speye(nx);
    grad   = -(kron(T, E) + kron(E, T));

    for i = 1:nx
        ii = (i-1)*nx + 1;
        jj = ii + nx - 1;
        grad(ii, jj) = 1.0;
        grad(jj, ii) = 1.0;
        kk = nxny - nx + i;
        grad(i,  kk) = 1.0;
        grad(kk,  i) = 1.0;
    end
    grad = grad / (dx * dy);
end


function [matdx, matdy] = gradient_mat(matx, Nx, Ny, dx, dy)
    format long;
    [matdx, matdy] = gradient(matx);

    matdx(1:Nx, 1)  = 0.5*(matx(1:Nx, 2)   - matx(1:Nx, Nx));
    matdx(1:Nx, Nx) = 0.5*(matx(1:Nx, 1)   - matx(1:Nx, Nx-1));
    matdy(1,  1:Ny) = 0.5*(matx(2,    1:Ny) - matx(Ny,   1:Ny));
    matdy(Ny, 1:Ny) = 0.5*(matx(1,    1:Ny) - matx(Ny-1, 1:Ny));

    matdx = matdx / dx;
    matdy = matdy / dy;
end


function [M] = vec2matx(V, N)
    format long;
    V = V.';
    V = V(:);
    R = ceil(length(V) / N);
    M = reshape(V, N, R).';
end


function [  ]= write_vtk_grid_values(nx,ny,dx,dy,istep,data1,data2,data3)

% open output file

fname=sprintf('time_%d.vtk',istep);
out =fopen(fname,'w');

nz=1;

npoin =nx*ny*nz;

% start writing ASCII VTK file:

% header of VTK file

fprintf(out,'# vtk DataFile Version 2.0\n');
fprintf(out,'time_10.vtk\n');
fprintf(out,'ASCII\n');
fprintf(out,'DATASET STRUCTURED_GRID\n');

%--- coords of grid points:

fprintf(out,'DIMENSIONS %5d  %5d  %5d\n',nx,ny,nz);
fprintf(out,'POINTS %7d   float64\n',npoin);

for i = 1:nx
for j = 1:ny

x =(i-1)*dx;
y =(j-1)*dy;
z = 0.0;

fprintf(out, '%14.6E   %14.6E  %14.6E\n',x,y,z);

end
end

%--- write grid point values:

fprintf(out,'POINT_DATA %5d\n',npoin);

fprintf(out,'SCALARS Phi  float64  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

for i = 1:nx
for j = 1:ny

fprintf(out,'%14.6E\n',data1(i,j));
end
end

fprintf(out,'SCALARS Conc  float64  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

for i = 1:nx
for j = 1:ny
fprintf(out,'%14.6E\n',data2(i,j));
end
end

fprintf(out,'SCALARS Theta  float64  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

for i = 1:nx
for j = 1:ny
fprintf(out,'%14.6E\n',data3(i,j));
end
end

fclose(out);

end
