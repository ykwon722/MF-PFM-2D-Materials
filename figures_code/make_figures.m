function make_figures(which)
% Paper figures from the saved runs in results/data.
% Output: results/figures/FigX.png (300 dpi) + FigX.pdf; FigX_masked.* = orientation panels with
% the substrate (phi<0) masked in grey. Snapshot info (time, flake area, coverage) -> figure_info.txt
if nargin < 1, which = {'F1','F2','F3','F4','F5','F6','S1','S2','F7'}; end
cfg = paths_config(); R = cfg.data; out = cfg.figures;
fi = fopen(fullfile(cfg.tables,'figure_info.txt'),'a'); fprintf(fi, '\n==== %s ====\n', datestr(now));
Lr = @(n) getfield(load(fullfile(R,[n '.mat'])),'o');
lab = 'abcdefghijklmnopqrstuvwxyz';
for w = which
  switch w{1}
    case 'F1'
      C = {Lr('FG_D5_T5_F2'), Lr('FG_D10_T5_F2')};
      F = fieldgrid(C, {'final','final'}); info(fi,'Fig1',{'D=5','D=10'},F);
      draw(out,'Fig1',F,{'phi','conc','theta'},[0 0.4],[]);
    case 'F2'
      C = {Lr('FG_D5_T2_F2'), Lr('FG_D5_T5_F2'), Lr('FG_D5_T10_F2')};
      F = fieldgrid(C, {'final','final','final'}); info(fi,'Fig2',{'tau=2','tau=5','tau=10'},F);
      drawrow(out,'Fig2',F);
    case 'F3'
      C = {Lr('FG_D10_T5_F2'), Lr('FG_D10_T5_F3'), Lr('FG_D10_T5_F6')};
      F = fieldgrid(C, {'final','final','final'}); info(fi,'Fig3',{'F=0.002','F=0.003','F=0.006'},F);
      drawrow(out,'Fig3',F);
    case 'F4'   % controlled pairs, same nucleation sequence (seed 1), t = 60
      C = {Lr('R_F4S_s1'), Lr('R_F4M_s1'), Lr('R_F4C_s1')};
      F = fieldgrid(C, {60,60,60}); info(fi,'Fig4',{'D=2 g=0.1','D=10 g=0.1','D=10 g=0.3'},F);
      draw(out,'Fig4',F,{'phi','conc','theta'},[0 0.4],[]);
      draw(out,'Fig4_masked',F,{'phi','conc','thetam'},[0 0.4],[]);
    case 'F5'
      C = {Lr('R_F5C_s1'), Lr('R_F5B_s1')};
      F = fieldgrid(C, {'final','final'}); info(fi,'Fig5',{'tau=100','tau=5'},F);
      draw(out,'Fig5',F,{'phi','conc','theta'},[],[]);
      draw(out,'Fig5_masked',F,{'phi','conc','thetam'},[],[]);
    case 'F6'
      C = {Lr('C6_base'), Lr('F6B')};
      F = fieldgrid(C, {100,100}); info(fi,'Fig6',{'symmetric','asymmetric'},F);
      crop = 128 + (-100:100);
      draw(out,'Fig6',F,{'phi','conc','theta'},[0.05 0.4],crop);
      draw(out,'Fig6_masked',F,{'phi','conc','thetam'},[0.05 0.4],crop);
    case {'S1','S2'}
      Dv = [2 5 10];
      if strcmp(w{1},'S1'), rv = [2 5 10]; nm = @(D,r) sprintf('FG_D%d_T%d_F2',D,r); rl = @(r) sprintf('\\tau_\\phi = %d',r);
      else,                 rv = [2 3 6];  nm = @(D,r) sprintf('FG_D%d_T5_F%d',D,r); rl = @(r) sprintf('F = 0.00%d',r); end
      f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 15 15.5]);
      tl = tiledlayout(3,3,'TileSpacing','none','Padding','compact'); k = 0;
      for r = 1:3, for c = 1:3
        k = k+1; o = Lr(nm(Dv(c),rv(r))); nexttile; imagesc(o.phi); clim([-1 1]); axis image; set(gca,'XTick',[],'YTick',[]);
        m = shape_metrics(o.phi);
        text(8,16,sprintf('(%s)',lab(k)),'Color','w','FontSize',9);
        text(250,248,sprintf('t = %.0f',o.tfinal),'Color','w','FontSize',7,'HorizontalAlignment','right');
        if r == 1, title(sprintf('D = %d',Dv(c)),'FontWeight','normal','FontSize',10); end
        if c == 1, ylabel(rl(rv(r)),'FontSize',10); end
        fprintf(fi,'%s (%s): %s t=%.1f area %.0f solidity %.3f\n', w{1}, lab(k), nm(Dv(c),rv(r)), o.tfinal, m.area, m.solidity);
      end, end
      colormap(f, parula); save_fig(f, out, w{1});
    case 'F7'
      collision_map(out, Lr, 60, [0 15 30 45 60], [0 5 10 15], 'Fig7');
      collision_map(out, Lr, 40, [0 30 60], [0 10], 'FigS_collision40');
  end
end
fclose(fi);
end

% ------------------------------------------------------------------------------------------
function F = fieldgrid(C, when)
F = struct('phi',{},'conc',{},'theta',{},'t',{});
for k = 1:numel(C)
  o = C{k};
  if ischar(when{k}), F(k).phi = double(o.phi); F(k).conc = double(o.conc); F(k).theta = double(o.theta);
      if isfield(o,'tfinal'), F(k).t = o.tfinal; else, F(k).t = o.p.tend; end
  else
      j = find(abs([o.snaps.t] - when{k}) < 1e-6, 1);
      F(k).phi = double(o.snaps(j).phi); F(k).conc = double(o.snaps(j).conc); F(k).theta = double(o.snaps(j).theta); F(k).t = o.snaps(j).t;
  end
end
end

function info(fi, name, lbl, F)
for k = 1:numel(F)
  ph = F(k).phi; cov = mean(ph(:) > 0);
  if max(bwlabel(ph>0),[],'all') == 1
    m = shape_metrics(ph); fprintf(fi,'%s row %d (%s): t=%.1f, flake area %.0f, solidity %.3f, R/Req %.2f\n', name, k, lbl{k}, F(k).t, m.area, m.solidity, m.Rmax_over_Req);
  else
    fprintf(fi,'%s row %d (%s): t=%.1f, coverage %.3f\n', name, k, lbl{k}, F(k).t, cov);
  end
end
end

function draw(out, name, F, cols, cl, crop)
lab = 'abcdefghijklmnopqrstuvwxyz'; nr = numel(F); nc = numel(cols);
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 5.6*nr]);
tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact'); k = 0;
for r = 1:nr, for c = 1:nc
  k = k+1; ax = nexttile;
  switch cols{c}
    case 'phi',   A = F(r).phi;   lim = [-1 1]; tk = -1:0.5:1;
    case 'conc',  A = F(r).conc;  if isempty(cl), lim = [min(A(:)) max(A(:))]; else, lim = cl; end; tk = [];
    case {'theta','thetam'}, A = F(r).theta; lim = [-pi pi]; tk = [-3 -1.5 0 1.5 3];
  end
  if ~isempty(crop), A = A(crop,crop); P = F(r).phi(crop,crop); else, P = F(r).phi; end
  if strcmp(cols{c},'thetam')
    h = imagesc(A); set(h,'AlphaData', double(P >= 0)); set(ax,'Color',[0.82 0.82 0.82]);
  else
    imagesc(A);
  end
  clim(lim); axis image; set(ax,'XTick',[],'YTick',[]); colormap(ax, parula);
  cb = colorbar; if ~isempty(tk), cb.Ticks = tk; end; cb.FontSize = 7;
  title(sprintf('(%s)',lab(k)),'FontWeight','normal','FontSize',11,'Units','normalized','Position',[-0.02 1.0],'HorizontalAlignment','left');
end, end
save_fig(f, out, name);
end

function drawrow(out, name, F)
lab = 'abc';
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 5.6]);
tiledlayout(1,numel(F),'TileSpacing','compact','Padding','compact');
for k = 1:numel(F)
  ax = nexttile; imagesc(F(k).phi); clim([-1 1]); axis image; set(ax,'XTick',[],'YTick',[]); colormap(ax, parula);
  cb = colorbar; cb.Ticks = -1:0.5:1; cb.FontSize = 7;
  title(sprintf('(%s)',lab(k)),'FontWeight','normal','FontSize',11,'Units','normalized','Position',[-0.02 1.0],'HorizontalAlignment','left');
end
save_fig(f, out, name);
end

function collision_map(out, Lr, a, PH, DL, name)
wr = @(x) mod(x+pi,2*pi)-pi;
gbm = @(th,ph) (abs(wr(circshift(th,[0 1])-th)) > 0.0873 | abs(wr(circshift(th,[1 0])-th)) > 0.0873) & imerode(ph>0, ones(3));
C = Lr('CC_base'); Af = @(t) interp1(C.At(:,1), C.At(:,2), t, 'linear', 'extrap');
f = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 3.6*numel(DL)+1.5 3.6*numel(PH)+1]);
tiledlayout(numel(PH),numel(DL),'TileSpacing','tight','Padding','compact');
for ip = 1:numel(PH), for id = 1:numel(DL)
  o = Lr(sprintf('OB_a%d_p%d_t%d', a, PH(ip), DL(id))); th = double(o.theta); ph = double(o.phi);
  tA = deg2rad(-30); tB = deg2rad(-30+a); w3 = deg2rad(min(a,60)/3);
  A1 = nnz(ph>0 & abs(wr(th-tA)) < w3); A2 = nnz(ph>0 & abs(wr(th-tB)) < w3);
  N = size(ph,1); cr = round(N/2) + (-85:85);
  ax = nexttile; h = imagesc(th(cr,cr)); set(h,'AlphaData',double(ph(cr,cr) >= 0)); set(ax,'Color',[1 1 1]);
  clim([-pi/2 pi/2]); colormap(ax, parula); axis image; set(ax,'XTick',[],'YTick',[]); hold on;
  g = gbm(th, ph); g = g(cr,cr); [r,c] = find(g); plot(c, r, 'k.', 'MarkerSize', 3);
  text(4, 165, sprintf('%.2f (%.2f)', A2/A1, Af(120-DL(id))/Af(120)), 'FontSize', 7);
  if ip == 1, title(sprintf('\\Deltat = %d', DL(id)),'FontWeight','normal','FontSize',9); end
  if id == 1, ylabel(sprintf('%d^\\circ', PH(ip)),'FontSize',9); end
end, end
save_fig(f, out, name);
end

function save_fig(f, out, name)
exportgraphics(f, fullfile(out,[name '.png']), 'Resolution', 300);
exportgraphics(f, fullfile(out,[name '.pdf']), 'ContentType', 'image', 'Resolution', 300);
close(f);
end
