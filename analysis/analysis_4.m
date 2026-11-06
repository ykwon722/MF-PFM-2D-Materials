function analysis_4()
% dt = 1e-4, 512 box, W x 0.75 sweep (partial)
cfg = paths_config(); d = cfg.data; dg = fullfile(cfg.figures, 'diagnostics'); if ~exist(dg,'dir'), mkdir(dg); end
old = cd(d); back = onCleanup(@() cd(old));
fo = fopen(fullfile(cfg.tables,'analysis4.txt'),'w'); P = @(varargin) fprintf(fo, varargin{:});
L = @(n) getfield(load(fullfile(d,[n '.mat'])),'o');
iou = @(a,b) nnz(a&b)/nnz(a|b);
for pr = {{'CS_base','CS_dt5','CS_dt1'}, {'CC_base','CC_dt5','CC_dt1'}}
  O = cellfun(L, pr{1}, 'UniformOutput', false);
  for s = 1:numel(O{1}.snaps)
    P('DT t=%g', O{1}.snaps(s).t); mr = shape_metrics(O{1}.snaps(s).phi);
    for k = 1:3
      j = find(abs([O{k}.snaps.t]-O{1}.snaps(s).t) < 1e-6, 1); m = shape_metrics(O{k}.snaps(j).phi);
      P(' | %s (dt %.0e) A %.0f sol %.4f R %.3f IoU %.4f', pr{1}{k}, O{k}.p.dt, m.area, m.solidity, m.Rmax_over_Req, iou(mr.mask, m.mask));
    end
    P('\n');
  end
  P('   wall: %s s\n', mat2str(round(cellfun(@(o) o.wall, O))));
end
A = L('CS_base'); B = L('CS_L512');
for s = 1:numel(A.snaps)
  ma = shape_metrics(A.snaps(s).phi); mb = shape_metrics(B.snaps(s).phi);
  P('BOX t=%g: 256 A %.0f sol %.3f R %.3f | 512 A %.0f sol %.3f R %.3f\n', A.snaps(s).t, ma.area, ma.solidity, ma.Rmax_over_Req, mb.area, mb.solidity, mb.Rmax_over_Req);
end
P('   wall CS_base %.0f s (2 threads, machine fully loaded), CS_L512 %.0f s\n', A.wall, B.wall);
% W x 0.75 sweep vs W x 1 at equal area
S = load('SW_descriptors.mat'); SW = S.SW;
w = dir('WSW_*.mat');
for i = 1:numel(w)
  nm = strrep(w(i).name, '.mat', ''); o = L(nm); dx = o.p.dx; M = [];
  for s = 1:numel(o.phis), if nnz(o.phis(s).phi>0) < 30, continue; end
    m = shape_metrics(o.phis(s).phi, dx); M(end+1,:) = [o.phis(s).t m.area m.solidity m.Rmax_over_Req]; end %#ok<AGROW>
  [~,u] = unique(M(:,2)); M = M(u,:); ref = SW.(strrep(nm,'WSW','SW'));
  th = double(o.theta(o.phi > 0.5));
  for Aq = [2000 4000 6000]
    if Aq > max(M(:,2)), continue; end
    v = interp1(M(:,2), M(:,3:4), Aq); vr = interp1(ref(:,2), ref(:,3:5), Aq);
    P('WSW %-16s A=%d: W0.75 sol %.3f R %.2f | W1 sol %.3f R %.2f\n', nm, Aq, v(1), v(2), vr(1), vr(3));
  end
  P('   theta in flake [%.4f %.4f], tfinal %.0f\n', min(th), max(th), o.tfinal);
end
fclose(fo);
end
