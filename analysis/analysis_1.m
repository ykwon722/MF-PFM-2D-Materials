function analysis_1()
% First analysis of finished revision runs: convergence (dt, rotation, Fig 6 rotation),
% bicrystals, Fig 4 realisations. Writes results/tables/analysis1.txt and results/figures/diagnostics/fig_*.png
cfg = paths_config(); d = cfg.data; dg = fullfile(cfg.figures, 'diagnostics'); if ~exist(dg,'dir'), mkdir(dg); end
old = cd(d); back = onCleanup(@() cd(old));
fo = fopen(fullfile(cfg.tables,'analysis1.txt'),'w'); P = @(varargin) fprintf(fo, varargin{:});
L = @(n) getfield(load(fullfile(d,[n '.mat'])),'o');
iou = @(a,b) nnz(a&b)/nnz(a|b);
wr = @(x) mod(x+pi,2*pi)-pi;
gbmask = @(th,ph) (abs(wr(circshift(th,[0 1])-th)) > 0.0873 | abs(wr(circshift(th,[1 0])-th)) > 0.0873) ...
                  & imerode(ph>0, ones(3));   % >5 deg jump inside flakes

% ---- time step: CS_base vs CS_dt5, CC_base vs CC_dt5
pairs = {'CS_base','CS_dt5'; 'CC_base','CC_dt5'};
for q = 1:2
  A = L(pairs{q,1}); B = L(pairs{q,2});
  for s = 1:numel(A.snaps)
    ma = shape_metrics(A.snaps(s).phi); mb = shape_metrics(B.snaps(s).phi);
    P('DT %s vs %s t=%g: area %.0f vs %.0f (%.2f%%) | sol %.3f vs %.3f | circ %.3f vs %.3f | Rmax/Req %.3f vs %.3f | IoU %.4f\n', ...
      pairs{q,1}, pairs{q,2}, A.snaps(s).t, ma.area, mb.area, 100*(mb.area/ma.area-1), ma.solidity, mb.solidity, ...
      ma.circularity, mb.circularity, ma.Rmax_over_Req, mb.Rmax_over_Req, iou(ma.mask,mb.mask));
  end
  P('   maxres %.2e / %.2e, maxphi %.6f / %.6f, minphi %.4f / %.4f, wall %.0f / %.0f s\n', A.maxres, B.maxres, A.maxphi, B.maxphi, A.minphi, B.minphi, A.wall, B.wall);
end

% ---- rotation: CS_r30, CS_r45 rotated back onto CS_base
A = L('CS_base'); mref = shape_metrics(A.phi);
for nm = {'CS_r30','CS_r45'}
  C = L(nm{1}); ang = str2double(nm{1}(5:end)); best = -1;
  for sgn = [-1 1], for dd = -2:0.5:2
    m = shape_metrics(C.phi,1,sgn*ang+dd); v = iou(mref.mask, m.mask); if v > best, best = v; ba = sgn*ang+dd; end
  end, end
  mc = shape_metrics(C.phi);
  P('ROT %s: area %.0f vs %.0f | sol %.3f vs %.3f | circ %.3f vs %.3f | Rmax/Req %.3f vs %.3f | IoU %.4f at %.1f deg\n', ...
    nm{1}, mc.area, mref.area, mc.solidity, mref.solidity, mc.circularity, mref.circularity, mc.Rmax_over_Req, mref.Rmax_over_Req, best, ba);
end

% ---- Fig 6 configuration under 45 deg rigid rotation
A = L('C6_base'); B = L('C6_rot'); mA = shape_metrics(A.phi);
best = -1; for dd = [-46:0.5:-44 44:0.5:46], m = shape_metrics(B.phi,1,dd); v = iou(mA.mask, m.mask); if v > best, best = v; bd = dd; end, end
thB = imrotate(double(B.theta) + pi/4, bd, 'nearest', 'crop'); phB = imrotate(double(B.phi), bd, 'bilinear', 'crop');
g1 = gbmask(double(A.theta), double(A.phi)); g2 = gbmask(thB, phB);
d12 = bwdist(g2); d21 = bwdist(g1);
P('FIG6ROT: phi IoU %.4f at %.1f deg | GB px %d vs %d | mean GB distance %.2f / %.2f cells | area %.0f vs %.0f\n', ...
  best, bd, nnz(g1), nnz(g2), mean(d12(g1)), mean(d21(g2)), sum(A.phi(:)>0), sum(B.phi(:)>0));
f = figure('Visible','off','Color','w','Position',[50 50 1100 360]);
th = double(A.theta); th(A.phi<0) = NaN; subplot(1,3,1); imagesc(th,'AlphaData',~isnan(th)); axis image off; clim([-pi pi]); title('C6\_base \theta (substrate masked)');
th = double(B.theta); th(B.phi<0) = NaN; subplot(1,3,2); imagesc(th,'AlphaData',~isnan(th)); axis image off; clim([-pi pi]); title('C6\_rot \theta');
rgb = ones([size(g1) 3]); rgb(:,:,2) = rgb(:,:,2)-0.8*g1; rgb(:,:,3) = rgb(:,:,3)-0.8*g1; rgb(:,:,1) = rgb(:,:,1)-0.8*g2; rgb(:,:,2) = rgb(:,:,2)-0.8*g2; rgb = max(rgb,0);
subplot(1,3,3); image(rgb); axis image off; hold on; contour(double(A.phi),[0 0],'k'); title('GB: base (red), rotated back (blue)');
exportgraphics(f, fullfile(dg, 'fig_C6rot.png'),'Resolution',110); close(f);

% ---- bicrystals
f = figure('Visible','off','Color','w','Position',[50 50 1200 600]); k = 0;
for a = [0 20 40 60]
  B = L(sprintf('BX_%d',a)); k = k+1; g = gbmask(double(B.theta), double(B.phi));
  m = shape_metrics(B.phi);
  P('BX %d: area %.0f | sol %.3f | circ %.3f | GB px %d | max |dtheta| across GB %.1f deg\n', a, m.area, m.solidity, m.circularity, nnz(g), ...
    rad2deg(max(abs(wr(B.theta(B.phi>0.9) - median(B.theta(B.phi>0.9)))))));
  subplot(2,4,k); imagesc(B.phi); axis image off; clim([-1 1]); title(sprintf('\\Delta\\theta = %d^\\circ: \\phi', a));
  th = double(B.theta); th(B.phi<0) = NaN; subplot(2,4,k+4); imagesc(th,'AlphaData',~isnan(th)); axis image off; clim([-pi/2 pi/2]); hold on; [r,c] = find(g); plot(c,r,'k.','MarkerSize',3); title('\theta, GB in black');
end
exportgraphics(f, fullfile(dg, 'fig_BX.png'),'Resolution',110); close(f);

% ---- Fig 4 realisations
cfgs = {'F4C','F4S','F4M'}; T = [60 80 100];
for c = 1:3
  S = zeros(5,4,numel(T));  % [coverage, islands, nuclei, GBpx]
  for s = 1:5
    nm = sprintf('R_%s_s%d', cfgs{c}, s);
    if exist([nm '.mat'],'file'), o = L(nm);
    elseif s == 1, o = L(nm);
    else, continue; end
    for it = 1:numel(T)
      j = find(abs([o.snaps.t]-T(it)) < 1e-6, 1); ph = double(o.snaps(j).phi); th = double(o.snaps(j).theta);
      cc = bwconncomp(ph>0); g = gbmask(th, ph);
      S(s,:,it) = [mean(ph(:)>0), cc.NumObjects, sum(o.NUC(:,6)==1), nnz(g)];
    end
  end
  for it = 1:numel(T)
    X = S(:,:,it); X = X(any(X,2),:);
    P('REAL %s t=%g (n=%d): coverage %.3f +- %.3f | islands %.1f +- %.1f | nuclei %.1f +- %.1f | GB px %.0f +- %.0f\n', cfgs{c}, T(it), size(X,1), ...
      mean(X(:,1)), std(X(:,1)), mean(X(:,2)), std(X(:,2)), mean(X(:,3)), std(X(:,3)), mean(X(:,4)), std(X(:,4)));
  end
end
% montage of t=60 for all realisations
f = figure('Visible','off','Color','w','Position',[50 50 1500 900]);
for c = 1:3, for s = 1:5
  nm = sprintf('R_%s_s%d', cfgs{c}, s);
  if exist([nm '.mat'],'file'), o = L(nm); elseif s==1, o = L(nm); else, continue; end
  j = find(abs([o.snaps.t]-60) < 1e-6, 1); th = double(o.snaps(j).theta); th(o.snaps(j).phi<0) = NaN;
  subplot(3,5,(c-1)*5+s); imagesc(th,'AlphaData',~isnan(th)); axis image off; clim([-pi pi]); title(sprintf('%s seed %d, t=60', cfgs{c}, s));
end, end
exportgraphics(f, fullfile(dg, 'fig_R_F4_t60.png'),'Resolution',90); close(f);
fclose(fo);
end
