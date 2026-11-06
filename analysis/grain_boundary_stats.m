function T = grain_boundary_stats(groups, times)
% GRAIN_BOUNDARY_STATS  Symmetry test of the orientation field (Scope and limitations, ESI).
% Compares runs with orientation differences taken modulo 2*pi (used in the paper) with
% runs in which they are reduced modulo 2*pi/3 (lattice symmetry) or 2*pi/6.
% For the final state of each run:
%   GB area   cells inside flakes (phi > 0.3 on both sides) whose neighbour differs in
%             orientation by more than 5 deg, counted modulo 2*pi, 2*pi/3 and 2*pi/6
%   grains    connected regions (>= 30 cells) of phi > 0.3 without grain boundary
%   pairs     neighbouring grains, classed by the difference of their mean orientations
%             modulo 120 deg (60 deg for the 2*pi/6 runs): twin > 50 deg, low angle < 10 deg
% The first run of each group is the reference for the shape overlap (IoU). times(c) is the
% time analysed for group c (the stored snapshot closest to it); Inf means the final state.
% The paper quotes GB_cells_mod120 for all runs (GB_cells_mod60 for the 2*pi/6 runs).
cfg = paths_config();
if nargin < 1
    groups = {{'R_F5C_s1','SYM_F5C_sym3','SYM_F5C_sym6'}, {'R_F5B_s1','SYM_F5B_sym3','SYM_F5B_sym6'}, ...
              {'R_F4C_s1','SYM_F4C_sym'}, {'R_F4S_s1','SYM_F4S_sym'}, {'C6_base','SYM_F6_sym'}};
    times = [Inf Inf 60 60 100];   % Fig. 5 final state, Fig. 4 at t = 60, Fig. 6 at t = 100
end
if ~exist('times', 'var'), times = Inf(1, numel(groups)); end
thr = deg2rad(5);
w = @(x,per) abs(mod(x + per/2, per) - per/2);
T = table();
for c = 1:numel(groups)
  ref = [];
  for v = 1:numel(groups{c})
    name = groups{c}{v}; f = fullfile(cfg.data, [name '.mat']);
    if ~exist(f, 'file'), fprintf('missing %s\n', name); continue; end
    o = getfield(load(f), 'o'); %#ok<GFLD>
    if isfield(o, 'tfinal'), tf = o.tfinal; else, tf = o.p.tend; end
    if isinf(times(c)) || abs(tf - times(c)) < 1e-6
        phi = double(o.phi); th = double(o.theta); tt = tf;
    else
        [~, is] = min(abs([o.snaps.t] - times(c))); tt = o.snaps(is).t;
        phi = double(o.snaps(is).phi); th = double(o.snaps(is).theta);
    end
    in = phi > 0.3;
    G = struct('m360',false(size(phi)),'m120',false(size(phi)),'m60',false(size(phi)));
    for sh = {[1 0],[0 1],[-1 0],[0 -1]}
      thn = circshift(th, sh{1}); both = in & circshift(in, sh{1}); dth = thn - th;
      G.m360 = G.m360 | (both & w(dth,2*pi) > thr);
      G.m120 = G.m120 | (both & w(dth,2*pi/3) > thr);
      G.m60  = G.m60  | (both & w(dth,2*pi/6) > thr);
    end
    six = endsWith(name, 'sym6');
    if six, own = G.m60; per = 2*pi/6; else, own = G.m120; per = 2*pi/3; end
    CC = bwconncomp(in & ~imdilate(own, ones(3)), 4);
    A = cellfun(@numel, CC.PixelIdxList); CC.PixelIdxList = CC.PixelIdxList(A >= 30);
    ng = numel(CC.PixelIdxList); mth = zeros(ng,1); Lb = zeros(size(phi));
    for q = 1:ng
      t = th(CC.PixelIdxList{q});
      mth(q) = angle(mean(exp(1i*t*2*pi/per)))*per/(2*pi); Lb(CC.PixelIdxList{q}) = q;
    end
    ntw = 0; nlow = 0; ngen = 0;
    for q = 1:ng
      nb = unique(Lb(imdilate(Lb==q, ones(9)) & Lb>0 & Lb~=q));
      for r = nb(nb>q)'
        d = rad2deg(w(mth(q)-mth(r), per));
        if d > 50, ntw = ntw + 1; elseif d < 10, nlow = nlow + 1; else, ngen = ngen + 1; end
      end
    end
    eta = (1 + phi)/2;
    if isempty(ref), ref = in; iou = 1; else, iou = nnz(in & ref)/nnz(in | ref); end
    row = table(string(name), tt, mean(eta(:)), iou, nnz(G.m360), nnz(G.m120), nnz(G.m60), ng, ntw, ngen, nlow, ...
      'VariableNames', {'run','t','coverage','IoU_vs_reference','GB_cells_mod360','GB_cells_mod120','GB_cells_mod60', ...
                        'grains','pairs_twin','pairs_general','pairs_lowangle'});
    T = [T; row]; %#ok<AGROW>
  end
end
writetable(T, fullfile(cfg.tables, 'grain_boundary_symmetry.csv'));
disp(T)
end
