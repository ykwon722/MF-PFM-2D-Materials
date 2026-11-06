function make_esi_figures()
% ESI figures for the revision: convergence, interface width, descriptor curves, bicrystals.
cfg = paths_config(); R = cfg.data; out = cfg.figures;
Lr = @(n) getfield(load(fullfile(R,[n '.mat'])),'o');
wr = @(x) mod(x+pi,2*pi)-pi;
% ---------- ESI convergence: contour overlays ----------
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 11.5]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
sets = {{'CS_base','CS_dt5','CS_dt1'}, {'CS_base','CS_dx05'}, {'CC_base','CC_dx05'}, ...
        {'CS_base','CS_r30','CS_r45'}, {'N6B_dx1','N6B_dx05'}, {'N6C_dx1','N6C_dx05'}};
ttl = {'(a) star: \Deltat = 10^{-3}, 5\times10^{-4}, 10^{-4}', '(b) star: \Deltax = 1, 0.5', '(c) compact n=3: \Deltax = 1, 0.5', ...
       '(d) star: seed at 0^\circ, 30^\circ, 45^\circ (rotated back)', '(e) branched n=6: \Deltax = 1, 0.5', '(f) compact n=6: \Deltax = 1, 0.5'};
col = {[0 0 0], [0.85 0.2 0.1], [0.1 0.4 0.85]}; lsty = {'-','--',':'};
for k = 1:6
  nexttile; hold on; S = sets{k};
  for j = 1:numel(S)
    o = Lr(S{j}); ph = double(o.phi);
    if k == 4 && j > 1, ang = str2double(S{j}(5:end)); ph = imrotate(ph + 1, ang, 'bilinear', 'crop') - 1; end
    N = size(ph,1); x = ((1:N) - N/2)*o.p.dx;
    contour(x, x, ph, [0 0], 'LineColor', col{j}, 'LineStyle', lsty{j}, 'LineWidth', 1.1);
  end
  axis image; set(gca,'YDir','reverse','FontSize',7); box on;
  lim = 110; if k == 6, lim = 45; elseif k == 3, lim = 75; end
  xlim([-lim lim]); ylim([-lim lim]); title(ttl{k}, 'FontWeight','normal','FontSize',8);
end
save_fig(f, out, 'FigS_convergence');

% ---------- ESI interface width ----------
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 11]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
P = {'WS_75','CS_base','WS_125'; 'WC_75','CC_base','WC_125'}; fw = [0.75 1 1.25]; lab = 'abcdef'; k = 0;
for r = 1:2, for c = 1:3
  k = k+1; o = Lr(P{r,c}); N = size(o.phi,1); x = ((1:N)-N/2)*o.p.dx; m = shape_metrics(o.phi, o.p.dx);
  nexttile; imagesc(x, x, o.phi); clim([-1 1]); axis image; colormap(gca, parula); set(gca,'FontSize',7);
  xlim([-110 110]); ylim([-110 110]);
  title(sprintf('(%s) W \\times %.2f: area %.0f, solidity %.2f', lab(k), fw(c), m.area, m.solidity), 'FontWeight','normal','FontSize',8);
end, end
save_fig(f, out, 'FigS_interface_width');

% ---------- ESI descriptor curves (W x 1 and W x 0.75) ----------
S = load(fullfile(R,'SW_descriptors.mat')); SW = S.SW;
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 6.2]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
G = {{'D2_T5_F2','D5_T5_F2','D10_T5_F2'}, {'D5_T2_F2','D5_T5_F2','D5_T10_F2'}, {'D10_T5_F2','D10_T5_F3','D10_T5_F6'}};
leg = {{'D = 2','D = 5','D = 10'}, {'\tau_\phi = 2','\tau_\phi = 5','\tau_\phi = 10'}, {'F = 0.002','F = 0.003','F = 0.006'}};
ttl = {'(a) D (\tau_\phi = 5, F = 0.002)', '(b) \tau_\phi (D = 5, F = 0.002)', '(c) F (D = 10, \tau_\phi = 5)'};
cc = lines(3);
for g = 1:3
  nexttile; hold on;
  for j = 1:3
    M = SW.(['SW_' G{g}{j}]); plot(M(:,2), M(:,3), '-', 'Color', cc(j,:), 'LineWidth', 1.4);
    o = Lr(['WSW_' G{g}{j}]); Mw = [];
    for s = 1:numel(o.phis), if nnz(o.phis(s).phi>0) < 30, continue; end
      m = shape_metrics(o.phis(s).phi, o.p.dx); Mw(end+1,:) = [m.area m.solidity]; end %#ok<AGROW>
    plot(Mw(:,1), Mw(:,2), '--', 'Color', cc(j,:), 'LineWidth', 1.1, 'HandleVisibility', 'off');
  end
  xlim([0 7000]); ylim([0.3 1]); grid on; box on; set(gca,'FontSize',7);
  xlabel('flake area'); ylabel('solidity'); title(ttl{g},'FontWeight','normal','FontSize',8);
  legend(leg{g}, 'Location','southwest','FontSize',6); text(200, 0.33, 'solid: W\times1, dashed: W\times0.75', 'FontSize', 6);
end
save_fig(f, out, 'FigS_descriptors');

% ---------- ESI symmetric bicrystals ----------
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 8.8]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact'); A = [0 20 40 60];
gbm = @(th,ph) (abs(wr(circshift(th,[0 1])-th)) > 0.0873 | abs(wr(circshift(th,[1 0])-th)) > 0.0873) & imerode(ph>0, ones(3));
for k = 1:4
  o = Lr(sprintf('BX_%d', A(k))); nexttile; imagesc(o.phi); clim([-1 1]); axis image off; colormap(gca, parula);
  title(sprintf('(%s) \\Delta\\theta = %d^\\circ', char('a'+k-1), A(k)), 'FontWeight','normal','FontSize',8);
end
for k = 1:4
  o = Lr(sprintf('BX_%d', A(k))); th = double(o.theta); ph = double(o.phi);
  nexttile; h = imagesc(th); set(h,'AlphaData',double(ph>=0)); set(gca,'Color',[1 1 1]); clim([-pi/2 pi/2]); axis image off; colormap(gca, parula); hold on;
  [r,c] = find(gbm(th,ph)); plot(c, r, 'k.', 'MarkerSize', 3);
  title(sprintf('(%s) orientation', char('e'+k-1)), 'FontWeight','normal','FontSize',8);
end
save_fig(f, out, 'FigS_bicrystals');
end

function save_fig(f, out, name)
exportgraphics(f, fullfile(out,[name '.png']), 'Resolution', 300);
exportgraphics(f, fullfile(out,[name '.pdf']), 'ContentType', 'vector');
close(f);
end
