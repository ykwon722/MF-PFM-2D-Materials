function p = case_params(name)
% CASE_PARAMS  Parameters of a named simulation case, for mfpfm_run.
% All cases use the defaults of mfpfm_run (orientation differences modulo 2*pi, orientation
% term weighted by (1+phi)/2, seed radius 2.5, random nuclei only on bare substrate).
% Names are <family>_<options>; the families are listed below and in case_list.m.
n3 = struct('n',3,'lambda',18,'conc_eq',0.1,'tau_s',200,'F',0.002);
n6 = struct('n',6,'lambda',200,'conc_eq',0.03,'tau_s',18,'F',0.007,'D',10,'g',0.04);
star = n3; star.tau_phi=2; star.D=5;  star.g=0.1; star.tend=185; star.snap_times=[60 120 185]; star.snap_every=5;
cmp  = n3; cmp.tau_phi=5;  cmp.D=10;  cmp.g=0.3; cmp.tend=100; cmp.snap_times=[50 100];     cmp.snap_every=5;
tok = regexp(name,'_','split');
switch tok{1}
  % ---------------- convergence (star = Fig 2a; compact = Fig 6 regime single flake) -----
  case 'CS'
    p = star; p.seeds = [0 0 0 0 2.5];
    switch tok{2}
      case 'base'
      case 'dt5',   p.dt = 5e-4; p.diag_every = 2000;
      case 'dt1',   p.dt = 1e-4; p.diag_every = 10000;
      case 'r30',   p.seeds(4) = deg2rad(30);
      case 'r45',   p.seeds(4) = deg2rad(45);
      case 'dx05',  p.dx = 0.5;
      case 'dx025', p.dx = 0.25;
      case 'L512',  p.L = 512;
    end
  case 'CC'
    p = cmp; p.seeds = [0 0 0 0 2.5];
    switch tok{2}
      case 'base'
      case 'dt5',   p.dt = 5e-4; p.diag_every = 2000;
      case 'dt1',   p.dt = 1e-4; p.diag_every = 10000;
      case 'dx05',  p.dx = 0.5;
    end
  case 'C6'   % Fig 6(a-c) configuration; 'rot' = rigid 45 deg rotation of the whole set-up
    p = cmp; p.snap_every = 0;
    if strcmp(tok{2},'base')
      p.seeds = [0 5 -5 deg2rad(40) 2.5; 0 -5 5 deg2rad(110) 2.5];
    else
      p.seeds = [0 5*sqrt(2) 0 deg2rad(40-45) 2.5; 0 -5*sqrt(2) 0 deg2rad(110-45) 2.5];
      p.theta_bg = 3 - pi/4;
    end
  case 'N6C'  % sixfold compact single flake (Fig 5a-c parameters), grid convergence
    p = n6; p.tau_phi = 100; p.tend = 245; p.seeds = [0 0 0 0 2.5]; p.snap_times = [100 150 200 245]; p.snap_every = 5;
    if strcmp(tok{2},'dx05'), p.dx = 0.5; end
  case 'N6B'  % sixfold branched single flake (Fig 5d-f parameters)
    p = n6; p.tau_phi = 5; p.tend = 125; p.seeds = [0 0 0 0 2.5]; p.snap_times = [50 75 100 125]; p.snap_every = 5;
    if strcmp(tok{2},'dx05'), p.dx = 0.5; end
  % ---------------- interface-width test: W -> f*W (kappa = k^2 ~ W), dx/W, lambda/W, tau/W^2 fixed
  case {'WS','WC'}
    if strcmp(tok{1},'WS'), p = star; else, p = cmp; end
    f = str2double(tok{2})/100;
    p.k = sqrt(2*f); p.dx = f; p.lambda = p.lambda*f; p.tau_phi = p.tau_phi*f^2;
    p.seeds = [0 0 0 0 2.5*f];
  % ---------------- equal-area single-parameter sweeps (Figs 1-3, S1, S2) -----------------
  case 'SW'   % SW_D<D>_T<tau>_F<1000F>
    p = n3; p.g = 0.1;
    p.D = str2double(tok{2}(2:end)); p.tau_phi = str2double(tok{3}(2:end)); p.F = str2double(tok{4}(2:end))/1000;
    p.seeds = [0 0 0 0 2.5]; p.snap_every = 5; p.stop_area = 7000; p.tend = 1000;
    pt = paper_time(p.D, p.tau_phi, p.F); p.snap_times = pt; if ~isempty(pt), p.t_min = max(pt); end
  % ---------------- random multi-flake realisations ---------------------------------------
  case 'R'    % R_<cfg>_s<seed>; Figs. 4 and 5 use seed 1 (one shared nucleation sequence), ESI Table S1 seeds 1-5
    s = str2double(tok{3}(2:end));
    switch tok{2}
      case 'F4C', p = n3; p.tau_phi=2; p.D=10; p.g=0.3; p.tend=100; p.snap_times=[17 21 40 60 80 100];
      case 'F4S', p = n3; p.tau_phi=2; p.D=2;  p.g=0.1; p.tend=100; p.snap_times=[17 21 40 60 80 100];
      case 'F4M', p = n3; p.tau_phi=2; p.D=10; p.g=0.1; p.tend=100; p.snap_times=[17 21 40 60 80 100];
      case 'F5C', p = n6; p.tau_phi=100; p.tend=245; p.snap_times=[100 150 200 245];
      case 'F5B', p = n6; p.tau_phi=5;   p.tend=125; p.snap_times=[50 75 100 125];
    end
    p.rgrain = true; p.rng_seed = s;
  % ---------------- symmetric bicrystals, misorientation a (deg), tip-to-tip family ---------
  case 'BX'
    a = str2double(tok{2});
    p = cmp; p.L = 128; p.tend = 60; p.snap_times = [30 60]; p.snap_every = 0;
    p.seeds = [0 0 -10 deg2rad(-a/2) 2.5; 0 0 10 deg2rad(a/2) 2.5];
  % ---------------- bicrystal outcome map: OM_a<mis>_t<delay>_o<offset> --------------------
  case 'OM'
    a = str2double(tok{2}(2:end)); dl = str2double(tok{3}(2:end)); of = str2double(tok{4}(2:end));
    p = cmp; p.L = 192; p.tend = 100; p.snap_times = [25 50 75 100]; p.snap_every = 0;
    p.seeds = [0 0 -15 deg2rad(-a/2) 2.5; dl of 15 deg2rad(a/2) 2.5];
  % ---------------- star grid convergence with uniform theta (theta_bg = 0 = flake orientation):
  % the explicit theta equation is unstable at dx = 0.25 with dt = 1e-3 (limit ~ dx^2), so the
  % phi-c convergence is tested with theta uniform, where the theta equation is inactive.
  case 'CSZ'
    p = star; p.seeds = [0 0 0 0 2.5]; p.theta_bg = 0;
    switch tok{2}, case 'dx1', case 'dx05', p.dx = 0.5; case 'dx025', p.dx = 0.25; end
  % ---------------- figure regeneration: single flakes stopped at equal area 6000 (Figs 1-3, S1, S2)
  case 'FG'   % FG_D<D>_T<tau>_F<1000F>
    p = n3; p.g = 0.1;
    p.D = str2double(tok{2}(2:end)); p.tau_phi = str2double(tok{3}(2:end)); p.F = str2double(tok{4}(2:end))/1000;
    p.seeds = [0 0 0 0 2.5]; p.stop_area = 6000; p.tend = 1000; p.diag_every = 100;
  case 'MG'   % MoS2 case: g from DFT edge-energy ratios, compact set (D = 10, tau = 5, F = 0.002), area 6000
    p = n3; p.D = 10; p.tau_phi = 5; p.g = str2double(tok{2})/100;
    p.seeds = [0 0 0 0 2.5]; p.stop_area = 6000; p.tend = 1000; p.diag_every = 100; p.snap_every = 5;
  case 'MT'   % MoS2 case, attachment-limited set (D = 5, tau = 10, F = 0.002), area 6000
    p = n3; p.D = 5; p.tau_phi = 10; p.g = str2double(tok{2})/100;
    p.seeds = [0 0 0 0 2.5]; p.stop_area = 6000; p.tend = 1000; p.diag_every = 100; p.snap_every = 5;
  case 'F6B'  % Fig 6(d-f): 29 deg at (-15,-18), t = 0; 105 deg at (18,15), t = 5
    p = cmp; p.snap_every = 0; p.snap_times = [50 100];
    p.seeds = [0 -15 -18 deg2rad(29) 2.5; 5 18 15 deg2rad(105) 2.5];
  % ---------------- W x 0.75 sweep: does the D / tau / F ordering survive at a thinner interface?
  case 'WSW'  % WSW_D<D>_T<tau>_F<1000F>
    f = 0.75; p = n3; p.g = 0.1;
    p.D = str2double(tok{2}(2:end)); tau0 = str2double(tok{3}(2:end)); p.F = str2double(tok{4}(2:end))/1000;
    p.k = sqrt(2*f); p.dx = f; p.lambda = p.lambda*f; p.tau_phi = tau0*f^2;
    p.seeds = [0 0 0 0 2.5*f]; p.snap_every = 5; p.stop_area = 6000; p.tend = 1000;
  % ---------------- bicrystal map v2: B placed at angle phi from A's tip direction (0 = tip-to-tip,
  % 60 = edge-to-edge), separation 20, delay dl, misorientation a
  case 'OB'   % OB_a<mis>_p<phi>_t<delay>
    a = str2double(tok{2}(2:end)); ph = str2double(tok{3}(2:end)); dl = str2double(tok{4}(2:end));
    p = cmp; p.L = 256; p.tend = 120; p.snap_times = [30 60 90 120]; p.snap_every = 0;
    ang = deg2rad(90 + ph); dd = 20;   % A's tip points along +j (angle 90 deg) for theta_A = -30 deg
    p.seeds = [0  -dd/2*cos(ang) -dd/2*sin(ang) deg2rad(-30)   2.5; ...
               dl  dd/2*cos(ang)  dd/2*sin(ang) deg2rad(-30+a) 2.5];
  % ---------------- symmetry test: SYM_<cfg>_<wrap>, as R_<cfg>_s1 (C6_base for F6) with the
  % orientation differences reduced modulo 2*pi/n ('sym', 'sym6') or 2*pi/3 ('sym3')
  case 'SYM'
    if strcmp(tok{2},'F6'), p = case_params('C6_base'); else, p = case_params(['R_' tok{2} '_s1']); end
    switch tok{3}, case {'sym','sym6'}, p.wrap = 'sym'; case 'sym3', p.wrap = 'sym3'; end
  otherwise
    error('unknown case %s', name);
end
end

function t = paper_time(D, tau, F)
% snapshot times used in the paper for Figs 1-3 (empty if the case is only in S1/S2)
T = [5 5 0.002 386; 10 5 0.002 260; 5 2 0.002 185; 5 5 0.002 295; 5 10 0.002 400; ...
     10 5 0.003 93; 10 5 0.006 40];
m = T(:,1)==D & T(:,2)==tau & abs(T(:,3)-F)<1e-9; t = T(m,4)';
end
