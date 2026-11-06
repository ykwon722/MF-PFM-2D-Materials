function analysis_5()
% sixfold compact grid, W x 0.75 sweep (complete), collision map v2, CSZ dx1 vs dx05
cfg = paths_config(); d = cfg.data; dg = fullfile(cfg.figures, 'diagnostics'); if ~exist(dg,'dir'), mkdir(dg); end
old = cd(d); back = onCleanup(@() cd(old));
fo = fopen(fullfile(cfg.tables,'analysis5.txt'),'w'); P = @(varargin) fprintf(fo, varargin{:});
L = @(n) getfield(load(fullfile(d,[n '.mat'])),'o');
wr = @(x) mod(x+pi,2*pi)-pi;
% ---- sixfold compact grid
A = L('N6C_dx1'); B = L('N6C_dx05');
for s = 1:numel(A.snaps)
  ma = shape_metrics(A.snaps(s).phi, 1); mb = shape_metrics(B.snaps(s).phi, 0.5);
  P('N6C t=%g: dx1 A %.0f sol %.3f circ %.3f R %.3f | dx0.5 A %.0f sol %.3f circ %.3f R %.3f\n', A.snaps(s).t, ma.area, ma.solidity, ma.circularity, ma.Rmax_over_Req, mb.area, mb.solidity, mb.circularity, mb.Rmax_over_Req);
end
% ---- CSZ
O = {L('CSZ_dx1'), L('CSZ_dx05')}; base = L('CS_base');
for s = 1:3
  m0 = shape_metrics(base.snaps(s).phi); m1 = shape_metrics(O{1}.snaps(s).phi); m2 = shape_metrics(O{2}.snaps(s).phi, 0.5);
  P('CSZ t=%g: CS_base A %.0f sol %.3f R %.3f | CSZ dx1 A %.0f sol %.3f R %.3f | CSZ dx0.5 A %.0f sol %.3f R %.3f\n', base.snaps(s).t, ...
    m0.area, m0.solidity, m0.Rmax_over_Req, m1.area, m1.solidity, m1.Rmax_over_Req, m2.area, m2.solidity, m2.Rmax_over_Req);
end
% ---- W x 0.75 sweep, all
S = load('SW_descriptors.mat'); SW = S.SW; w = dir('WSW_*.mat');
for i = 1:numel(w)
  nm = strrep(w(i).name, '.mat', ''); o = L(nm); M = [];
  for s = 1:numel(o.phis), if nnz(o.phis(s).phi>0) < 30, continue; end
    m = shape_metrics(o.phis(s).phi, o.p.dx); M(end+1,:) = [o.phis(s).t m.area m.solidity m.Rmax_over_Req]; end %#ok<AGROW>
  [~,u] = unique(M(:,2)); M = M(u,:); ref = SW.(strrep(nm,'WSW','SW'));
  for Aq = [2000 4000]
    v = interp1(M(:,2), M(:,3:4), Aq); vr = interp1(ref(:,2), ref(:,3:5), Aq);
    P('WSW %-16s A=%d: W0.75 sol %.3f R %.2f | W1 sol %.3f R %.2f\n', nm, Aq, v(1), v(2), vr(1), vr(3));
  end
  th = double(o.theta(o.phi > 0.5)); P('   theta in flake [%.4f %.4f]\n', min(th), max(th));
end
% ---- collision map v2
C = L('CC_base'); Af = @(t) interp1(C.At(:,1), C.At(:,2), t, 'linear', 'extrap');
gbm = @(th,ph) (abs(wr(circshift(th,[0 1])-th)) > 0.0873 | abs(wr(circshift(th,[1 0])-th)) > 0.0873) & imerode(ph>0, ones(3));
for a = [60 40]
  if a == 60, PH = [0 15 30 45 60]; DL = [0 5 10 15]; else, PH = [0 30 60]; DL = [0 10]; end
  f = figure('Visible','off','Color','w','Position',[50 50 260*numel(DL) 240*numel(PH)]); k = 0;
  for ph = PH, for dl = DL
    k = k + 1; nm = sprintf('OB_a%d_p%d_t%d', a, ph, dl);
    if ~exist([nm '.mat'],'file'), continue; end
    o = L(nm); th = double(o.theta); phv = double(o.phi);
    tA = deg2rad(-30); tB = deg2rad(-30+a);
    A1 = nnz(phv>0 & abs(wr(th - tA)) < deg2rad(min(a,60)/3)); A2 = nnz(phv>0 & abs(wr(th - tB)) < deg2rad(min(a,60)/3));
    g = gbm(th, phv); cc = bwconncomp(phv > 0); mm = shape_metrics(o.phi);
    exp_ratio = Af(120 - dl)/Af(120);
    P('OB %s: A1 %d A2 %d ratio %.3f (indep %.3f) | GB px %d | islands %d | union sol %.3f circ %.3f | placed %s\n', nm, A1, A2, A2/A1, exp_ratio, nnz(g), cc.NumObjects, mm.solidity, mm.circularity, mat2str(o.SEEDLOG(:,3)'));
    thp = th; thp(phv < 0) = NaN;
    subplot(numel(PH), numel(DL), k); imagesc(thp, 'AlphaData', ~isnan(thp)); axis image off; clim([-pi/2 pi/2]); hold on;
    [r,c] = find(g); plot(c, r, 'k.', 'MarkerSize', 2);
    N = size(phv,1); xlim(N/2 + [-90 90]); ylim(N/2 + [-90 90]);
    title(sprintf('\\phi_{pos}=%d^\\circ, \\Deltat=%d: A_2/A_1=%.2f (%.2f)', ph, dl, A2/A1, exp_ratio), 'FontSize', 7);
  end, end
  exportgraphics(f, fullfile(dg, sprintf('fig_OB_a%d.png', a)), 'Resolution', 90); close(f);
end
fclose(fo);
end
