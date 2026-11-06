function R = nucleus_survival()
% NUCLEUS_SURVIVAL  ESI Table S2. A single seed of radius r is placed on a uniform substrate
% at the far-field concentration c_inf and grown for t = 10 on a 64 x 64 domain. The seed
% counts as grown if max(phi) > 0.5 at the end.
cfg = paths_config();
R = struct('cfg',{},'r',{},'cinf',{},'maxphi',{},'area',{},'grew',{});
cases = {
  'n3', struct('n',3,'lambda',18,'conc_eq',0.1,'tau_s',200,'tau_phi',2,'D',10,'g',0.3), [0.4 0.2 0.15 0.12 0.11];
  'n6', struct('n',6,'lambda',200,'conc_eq',0.03,'tau_s',18,'tau_phi',5,'D',10,'g',0.04), [0.126 0.06 0.045 0.035 0.032]};
for c = 1:size(cases,1)
    for cinf = cases{c,3}
        for r = [1 2.5 4]
            p = cases{c,2}; p.F = cinf/p.tau_s; p.L = 64; p.tend = 10;
            p.seeds = [0 0 0 0 r];
            o = mfpfm_run(p); eta = (1 + double(o.phi))/2;
            R(end+1) = struct('cfg',cases{c,1},'r',r,'cinf',cinf,'maxphi',o.maxphi,'area',sum(eta(:)),'grew',o.maxphi > 0.5); %#ok<AGROW>
        end
    end
end
save(fullfile(cfg.tables, 'nucleus_survival.mat'), 'R');
writetable(struct2table(R), fullfile(cfg.tables, 'nucleus_survival.csv'));
end
