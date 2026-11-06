function cfg = paths_config()
% PATHS_CONFIG  Folder layout of this package; also adds the code folders to the MATLAB path.
% The run, analysis and figure scripts call it themselves.
root = fileparts(mfilename('fullpath'));
cfg.root    = root;
cfg.data    = fullfile(root, 'results', 'data');      % one .mat file per simulation case
cfg.figures = fullfile(root, 'results', 'figures');   % figures of the paper and ESI (PNG and PDF)
cfg.tables  = fullfile(root, 'results', 'tables');    % numbers quoted in the paper, ESI and response
addpath(fullfile(root,'solver'), fullfile(root,'cases'), fullfile(root,'analysis'), fullfile(root,'figures_code'));
for f = {cfg.data, cfg.figures, cfg.tables}
    if ~exist(f{1}, 'dir'), mkdir(f{1}); end
end
end
