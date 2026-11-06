% RUN_ALL  Run every simulation case used in the paper and the ESI, then rebuild the
% tables and figures. Cases whose .mat file already exists in results/data are skipped,
% so the script can be stopped and restarted. Running all cases takes several days on
% one desktop computer; to run one group only, set GROUPS below, e.g. GROUPS = {'Fig6'}.
cfg = paths_config();
C = case_list();
GROUPS = fieldnames(C);
for g = 1:numel(GROUPS)
    names = C.(GROUPS{g});
    for k = 1:numel(names)
        fprintf('%s: %s\n', GROUPS{g}, names{k});
        run_case(names{k});
    end
end
nucleus_survival();      % ESI Table S2
analysis_1(); analysis_2(); analysis_3(); analysis_4(); analysis_5();
grain_boundary_stats();  % symmetry test (Scope and limitations, ESI)
make_figures(); make_esi_figures(); make_mos2_figure();
