function analysis_3()
% Third analysis: grid convergence (star dx 1/0.5/0.25, compact, sixfold), F5 realisations.
cfg = paths_config(); d = cfg.data; dg = fullfile(cfg.figures, 'diagnostics'); if ~exist(dg,'dir'), mkdir(dg); end
old = cd(d); back = onCleanup(@() cd(old));
fo = fopen(fullfile(cfg.tables,'analysis3.txt'),'w'); P = @(varargin) fprintf(fo, varargin{:});
L = @(n) getfield(load(fullfile(d,[n '.mat'])),'o');
wr = @(x) mod(x+pi,2*pi)-pi; wr60 = @(x) mod(x+pi/6,pi/3)-pi/6;
iou = @(a,b) nnz(a&b)/nnz(a|b);

% ---------- grid convergence at equal time and equal area ----------
sets = {{'CS_base','CS_dx05','CS_dx025'}, {'CC_base','CC_dx05'}, {'N6B_dx1','N6B_dx05'}};
if exist('N6C_dx05.mat','file'), sets{end+1} = {'N6C_dx1','N6C_dx05'}; else, sets{end+1} = {'N6C_dx1'}; end
for q = 1:numel(sets)
  O = cellfun(L, sets{q}, 'UniformOutput', false);
  ref = O{1};
  for s = 1:numel(ref.snaps)
    P('GRID t=%g:', ref.snaps(s).t);
    mr = shape_metrics(ref.snaps(s).phi, ref.p.dx);
    for k = 1:numel(O)
      j = find(abs([O{k}.snaps.t] - ref.snaps(s).t) < 1e-6, 1);
      m = shape_metrics(O{k}.snaps(j).phi, O{k}.p.dx);
      f = round(ref.p.dx/O{k}.p.dx); pb = imresize(double(O{k}.snaps(j).phi), 1/f, 'bilinear'); mb = shape_metrics(pb, ref.p.dx);
      P(' | %s A %.0f sol %.3f circ %.3f R %.3f IoU %.3f', sets{q}{k}, m.area, m.solidity, m.circularity, m.Rmax_over_Req, iou(mr.mask, mb.mask));
    end
    P('\n');
  end
  % equal area from phis
  for Aq = [2000 4000 6000]
    P('GRID A=%d:', Aq);
    for k = 1:numel(O)
      ph = O{k}.phis; if isempty(ph), continue; end
      M = zeros(numel(ph),4);
      for s = 1:numel(ph), if nnz(ph(s).phi>0) < 30, M(s,:) = NaN; continue; end
        m = shape_metrics(ph(s).phi, O{k}.p.dx); M(s,:) = [ph(s).t m.area m.solidity m.Rmax_over_Req]; end
      M = M(~isnan(M(:,1)),:); [~,u] = unique(M(:,2)); M = M(u,:);
      if Aq > max(M(:,2)), continue; end
      v = interp1(M(:,2), M(:,[1 3 4]), Aq);
      P(' | %s t %.1f sol %.3f R %.3f', sets{q}{k}, v(1), v(2), v(3));
    end
    P('\n');
  end
  for k = 1:numel(O), P('   %s: dx %.2f wall %.0f s maxres %.1e maxphi %.5f\n', sets{q}{k}, O{k}.p.dx, O{k}.wall, O{k}.maxres, O{k}.maxphi); end
end
% figure: star contours dx 1 / 0.5 / 0.25
f = figure('Visible','off','Color','w','Position',[50 50 1200 420]);
nm = {'CS_base','CS_dx05','CS_dx025'}; col = {'k','b','r'};
subplot(1,2,1); hold on;
for k = 1:3, o = L(nm{k}); N = size(o.phi,1); x = ((1:N)-N/2)*o.p.dx; contour(x, x, double(o.phi), [0 0], col{k}, 'LineWidth', 1.2); end
axis image; xlim([-110 110]); ylim([-110 110]); title('Star flake, t = 185: \Deltax = 1 (k), 0.5 (b), 0.25 (r)'); set(gca,'YDir','reverse');
subplot(1,2,2); hold on;
for k = 1:3, o = L(nm{k}); ph = o.phis; M = [];
  for s = 1:numel(ph), if nnz(ph(s).phi>0) < 30, continue; end, m = shape_metrics(ph(s).phi, o.p.dx); M(end+1,:) = [m.area m.solidity]; end %#ok<AGROW>
  plot(M(:,1), M(:,2), col{k}, 'LineWidth', 1.5); end
xlabel('flake area'); ylabel('solidity'); legend('\Deltax = 1','0.5','0.25'); grid on; title('Solidity vs area');
exportgraphics(f, fullfile(dg, 'fig_grid_star.png'),'Resolution',100); close(f);

% ---------- sixfold realisations ----------
gbm = @(th,ph,w) (abs(w(circshift(th,[0 1])-th)) > 0.0873 | abs(w(circshift(th,[1 0])-th)) > 0.0873) & imerode(ph>0, ones(3));
for cfg = {'F5C','F5B'}
  R = [];
  for s = 1:5
    o = L(sprintf('R_%s_s%d', cfg{1}, s));
    ph = double(o.phi); th = double(o.theta);
    [~, n] = plabel(ph > 0); g = gbm(th, ph, wr); gc = gbm(th, ph, wr60);
    R(end+1,:) = [mean(ph(:)>0), n, sum(o.NUC(:,6)==1), nnz(g), nnz(gc)]; %#ok<AGROW>
  end
  P('REAL %s final: coverage %.3f+-%.3f | islands %.1f+-%.1f | nuclei %.1f+-%.1f | GB px (raw >5deg) %.0f+-%.0f | crystallographic GB px (mod 60 >5deg) %.0f+-%.0f\n', ...
    cfg{1}, mean(R(:,1)), std(R(:,1)), mean(R(:,2)), std(R(:,2)), mean(R(:,3)), std(R(:,3)), mean(R(:,4)), std(R(:,4)), mean(R(:,5)), std(R(:,5)));
end
fclose(fo);
end

function [lab, n] = plabel(bw)
[l, n0] = bwlabel(bw, 4); par = 1:n0;
    function r = fnd(x), while par(x) ~= x, x = par(x); end, r = x; end
N = size(bw,1); M = size(bw,2);
for i = 1:N, if l(i,1) && l(i,M), a = fnd(l(i,1)); b = fnd(l(i,M)); if a ~= b, par(max(a,b)) = min(a,b); end, end, end
for j = 1:M, if l(1,j) && l(N,j), a = fnd(l(1,j)); b = fnd(l(N,j)); if a ~= b, par(max(a,b)) = min(a,b); end, end, end
root = arrayfun(@fnd, 1:n0); [~,~,id] = unique(root); lab = zeros(size(l));
lab(l>0) = id(l(l>0)); n = max([id(:); 0]);
end
