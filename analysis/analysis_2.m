function analysis_2()
% Second analysis: equal-area sweeps, interface-width test, outcome map, Fig 4 per-island stats.
cfg = paths_config(); d = cfg.data; dg = fullfile(cfg.figures, 'diagnostics'); if ~exist(dg,'dir'), mkdir(dg); end
old = cd(d); back = onCleanup(@() cd(old));
fo = fopen(fullfile(cfg.tables,'analysis2.txt'),'w'); P = @(varargin) fprintf(fo, varargin{:});
L = @(n) getfield(load(fullfile(d,[n '.mat'])),'o');
wr = @(x) mod(x+pi,2*pi)-pi;

% ================= 1. equal-area sweeps ==================
Aeq = [2000 4000 6000];
names = {}; for D = [2 5 10], for T = [2 5 10], names{end+1} = sprintf('SW_D%d_T%d_F2',D,T); end, end %#ok<AGROW>
for D = [2 5 10], for F = [3 6], names{end+1} = sprintf('SW_D%d_T5_F%d',D,F); end, end %#ok<AGROW>
SW = struct();
for i = 1:numel(names)
  o = L(names{i}); ns = numel(o.phis); M = nan(ns,5);
  for s = 1:ns
    if nnz(o.phis(s).phi > 0) < 30, continue; end
    m = shape_metrics(o.phis(s).phi); M(s,:) = [o.phis(s).t m.area m.solidity m.circularity m.Rmax_over_Req];
  end
  M = M(~isnan(M(:,1)),:); [~,u] = unique(M(:,2)); M = M(u,:);
  SW.(names{i}) = M;
  vals = interp1(M(:,2), M(:,3:5), Aeq, 'linear', NaN);
  P('SW %-14s tfinal %5.0f Afinal %5.0f | ', names{i}, o.tfinal, M(end,2));
  for a = 1:numel(Aeq), P('A=%d: t=%5.1f sol %.3f circ %.3f R %.2f | ', Aeq(a), interp1(M(:,2),M(:,1),Aeq(a)), vals(a,1), vals(a,2), vals(a,3)); end
  P('\n');
  for k = 1:numel(o.snaps)
    m = shape_metrics(o.snaps(k).phi);
    P('     paper snapshot t=%g: area %.0f sol %.3f circ %.3f R %.2f\n', o.snaps(k).t, m.area, m.solidity, m.circularity, m.Rmax_over_Req);
  end
end
save('SW_descriptors.mat','SW');
% descriptor-vs-area curves
f = figure('Visible','off','Color','w','Position',[50 50 1300 400]);
groups = {{'SW_D2_T5_F2','SW_D5_T5_F2','SW_D10_T5_F2'}, {'SW_D5_T2_F2','SW_D5_T5_F2','SW_D5_T10_F2'}, {'SW_D10_T5_F2','SW_D10_T5_F3','SW_D10_T5_F6'}};
ttl = {'D sweep (\tau=5, F=0.002)','\tau_\phi sweep (D=5, F=0.002)','F sweep (D=10, \tau=5)'};
for gI = 1:3
  subplot(1,3,gI); hold on; for k = 1:3, M = SW.(groups{gI}{k}); plot(M(:,2), M(:,3), '-', 'LineWidth', 1.5); end
  xlabel('flake area'); ylabel('solidity'); ylim([0.4 1]); xlim([0 8000]); grid on; title(ttl{gI});
  legend(strrep(groups{gI},'_','\_'),'Location','southwest','FontSize',7);
end
exportgraphics(f, fullfile(dg, 'fig_SW_solidity.png'),'Resolution',110); close(f);
% equal-area montages (S1 and S2 grids) at A = 6000
for gg = 1:2
  f = figure('Visible','off','Color','w','Position',[50 50 900 900]);
  for r = 1:3, for c = 1:3
    D = [2 5 10]; D = D(c);
    if gg == 1, T = [2 5 10]; nm = sprintf('SW_D%d_T%d_F2',D,T(r)); else, F = [2 3 6]; nm = sprintf('SW_D%d_T5_F%d',D,F(r)); end
    o = L(nm); A = arrayfun(@(s) nnz(s.phi>0), o.phis); [~,j] = min(abs(A - 6000));
    subplot(3,3,(r-1)*3+c); imagesc(o.phis(j).phi); axis image off; clim([-1 1]);
    m = shape_metrics(o.phis(j).phi);
    title(sprintf('%s: t=%g, A=%.0f, sol %.2f', strrep(nm(4:end),'_',' '), o.phis(j).t, m.area, m.solidity),'FontSize',8);
  end, end
  exportgraphics(f, fullfile(dg, sprintf('fig_S%d_equalarea.png',gg)),'Resolution',90); close(f);
end

% ================= 2. interface-width test ==================
pairs = {'CS_base','WS_75','WS_125'; 'CC_base','WC_75','WC_125'};
for q = 1:2
  for k = 1:3
    o = L(pairs{q,k}); dx = o.p.dx;
    for s = 1:numel(o.snaps)
      m = shape_metrics(o.snaps(s).phi, dx);
      P('W %-7s (dx=%.2f) t=%g: area %.0f sol %.3f circ %.3f R %.3f\n', pairs{q,k}, dx, o.snaps(s).t, m.area, m.solidity, m.circularity, m.Rmax_over_Req);
    end
  end
end
f = figure('Visible','off','Color','w','Position',[50 50 1100 700]);
for q = 1:2, for k = 1:3
  o = L(pairs{q,k}); dx = o.p.dx; N = size(o.phi,1); x = ((1:N)-N/2)*dx;
  subplot(2,3,(q-1)*3+k); imagesc(x, x, o.phi); axis image; clim([-1 1]); xlim([-110 110]); ylim([-110 110]);
  m = shape_metrics(o.phi, dx); title(sprintf('%s: A=%.0f, sol %.3f', strrep(pairs{q,k},'_','\_'), m.area, m.solidity));
end, end
exportgraphics(f, fullfile(dg, 'fig_Wtest.png'),'Resolution',90); close(f);

% ================= 3. outcome map ==================
gbm = @(th,ph) (abs(wr(circshift(th,[0 1])-th)) > 0.0873 | abs(wr(circshift(th,[1 0])-th)) > 0.0873) & imerode(ph>0, ones(3));
for a = [60 40]
  f = figure('Visible','off','Color','w','Position',[50 50 1300 1000]); k = 0;
  for of = [0 6 12], for dl = [0 5 10 15]
    k = k+1; nm = sprintf('OM_a%d_t%d_o%d', a, dl, of); o = L(nm);
    th = double(o.theta); ph = double(o.phi); th(ph<0) = NaN; g = gbm(double(o.theta), ph);
    subplot(3,4,k); imagesc(th,'AlphaData',~isnan(th)); axis image off; clim([-pi/2 pi/2]); hold on;
    [r,c] = find(g); plot(c,r,'k.','MarkerSize',2);
    A1 = nnz(ph>0 & abs(wr(double(o.theta)+deg2rad(a/2))) < deg2rad(a/4)); A2 = nnz(ph>0 & abs(wr(double(o.theta)-deg2rad(a/2))) < deg2rad(a/4));
    title(sprintf('\\Deltat=%d, offset=%d: A2/A1=%.2f', dl, of, A2/A1),'FontSize',8);
    P('OM %s: area %d, A1 %d, A2 %d, ratio %.3f, GB px %d, seeds placed %s\n', nm, nnz(ph>0), A1, A2, A2/A1, nnz(g), mat2str(o.SEEDLOG(:,3)'));
  end, end
  exportgraphics(f, fullfile(dg, sprintf('fig_OM_a%d.png',a)),'Resolution',90); close(f);
end

% ================= 4. Fig 4 realisations: periodic islands, per-island solidity ==================
cfgs = {'F4C','F4S','F4M'};
for c = 1:3
  for t = [60 100]
    R = [];
    for s = 1:5
      nm = sprintf('R_%s_s%d', cfgs{c}, s);
      if exist([nm '.mat'],'file'), o = L(nm);
      else, o = L(nm); end
      j = find(abs([o.snaps.t]-t) < 1e-6, 1); ph = double(o.snaps(j).phi); th = double(o.snaps(j).theta);
      [lab, n] = plabel(ph > 0); g = gbm(th, ph); sol = [];
      for q = 1:n
        mk = lab == q; if nnz(mk) < 50, continue; end
        [r0, c0] = find(mk, 1); mk = circshift(mk, [round(size(mk,1)/2)-r0, round(size(mk,2)/2)-c0]);
        st = regionprops(bwlabel(mk), 'Area', 'Solidity'); sol(end+1) = sum([st.Area].*[st.Solidity])/sum([st.Area]); %#ok<AGROW>
      end
      R(end+1,:) = [mean(ph(:)>0), n, sum(o.NUC(:,6)==1), nnz(g), median(sol)]; %#ok<AGROW>
    end
    P('REAL %s t=%d: coverage %.3f+-%.3f | islands %.1f+-%.1f | nuclei %.1f+-%.1f | GB px %.0f+-%.0f | median island solidity %.3f+-%.3f\n', ...
      cfgs{c}, t, mean(R(:,1)), std(R(:,1)), mean(R(:,2)), std(R(:,2)), mean(R(:,3)), std(R(:,3)), mean(R(:,4)), std(R(:,4)), mean(R(:,5)), std(R(:,5)));
  end
end
fclose(fo);
end

function [lab, n] = plabel(bw)
% connected components with periodic boundaries
[l, n0] = bwlabel(bw, 4); par = 1:n0;
    function r = fnd(x), while par(x) ~= x, x = par(x); end, r = x; end
N = size(bw,1); M = size(bw,2);
for i = 1:N, if l(i,1) && l(i,M), a = fnd(l(i,1)); b = fnd(l(i,M)); if a ~= b, par(max(a,b)) = min(a,b); end, end, end
for j = 1:M, if l(1,j) && l(N,j), a = fnd(l(1,j)); b = fnd(l(N,j)); if a ~= b, par(max(a,b)) = min(a,b); end, end, end
root = arrayfun(@fnd, 1:n0); [~,~,id] = unique(root); lab = zeros(size(l));
lab(l>0) = id(l(l>0)); n = max([id(:); 0]);
end
