function o = run_case(name, varargin)
% RUN_CASE  Run one simulation case and save it as results/data/<name>.mat.
%   o = run_case('FG_D5_T2_F2')              run and save; loads the saved file if it exists
%   o = run_case(name, 'overwrite', true)    run again even if the file exists
%   o = run_case(name, 'test', true)         short run to t = 0.05 that is not saved
% The parameters of every case are set in case_params.m.
cfg = paths_config();
opt = struct('overwrite', false, 'test', false, 'outdir', cfg.data);
for k = 1:2:numel(varargin), opt.(varargin{k}) = varargin{k+1}; end
f = fullfile(opt.outdir, [name '.mat']);
if ~opt.test && ~opt.overwrite && exist(f, 'file')
    o = getfield(load(f), 'o'); %#ok<GFLD>
    return
end
p = case_params(name);
if opt.test
    p.tend = 0.05; p.snap_times = []; p.snap_every = 0;
    if isfield(p, 'seeds') && ~isempty(p.seeds), p.seeds(:,1) = min(p.seeds(:,1), 0.01); end
else
    p.progfile = fullfile(opt.outdir, [name '.prog']);
end
o = mfpfm_run(p);
if ~opt.test
    save(f, 'o', '-v7.3');
    if exist(p.progfile, 'file'), delete(p.progfile); end
end
end
