function make_mos2_figure()
% Two-row MoS2 figure: row 1 attachment-limited (D = 5, tau = 10), row 2 diffusion-influenced (D = 10, tau = 5)
cfg = paths_config(); R = cfg.data; out = cfg.figures;
G = [35 20 10 3];
hdr = {{'g = 0.35','S-rich, \Delta\mu_S \approx 0'}, {'g = 0.20','\Delta\mu_S \approx -0.3 (C), -0.25 to -0.75 (S)'}, ...
       {'g = 0.10','\Delta\mu_S \approx -0.6 (C), -1.3 (S)'}, {'g = 0.03','\Delta\mu_S \approx -0.95 (C)'}};
rows = {'MT','MG'}; rl = {'D = 5, \tau_\phi = 10','D = 10, \tau_\phi = 5'};
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17.5 9.6]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
fi = fopen(fullfile(cfg.tables,'mos2_info.txt'),'w'); k = 0;
for r = 1:2, for c = 1:4
  k = k+1; o = getfield(load(fullfile(R,sprintf('%s_%02d.mat',rows{r},G(c)))),'o'); m = shape_metrics(o.phi);
  nexttile; imagesc(o.phi(28:228,28:228)); clim([-1 1]); axis image; set(gca,'XTick',[],'YTick',[]); colormap(gca, parula);
  if r == 1, L = hdr{c}; L{1} = sprintf('(%s) %s', char('a'+k-1), L{1}); title(L,'FontWeight','normal','FontSize',7);
  else, title(sprintf('(%s)', char('a'+k-1)),'FontWeight','normal','FontSize',7); end
  if c == 1, ylabel(rl{r},'FontSize',8); end
  text(196, 194, sprintf('t = %.0f, solidity %.2f', o.tfinal, m.solidity), 'Color','w','FontSize',6,'HorizontalAlignment','right');
  fprintf(fi,'%s g=%.2f: t=%.1f area %.0f solidity %.3f circ %.3f R/Req %.2f\n', rl{r}, G(c)/100, o.tfinal, m.area, m.solidity, m.circularity, m.Rmax_over_Req);
end, end
fclose(fi);
exportgraphics(f, fullfile(out,'FigS_MoS2.png'), 'Resolution', 300);
exportgraphics(f, fullfile(out,'FigS_MoS2.pdf'), 'ContentType', 'image', 'Resolution', 300);
close(f);
end
