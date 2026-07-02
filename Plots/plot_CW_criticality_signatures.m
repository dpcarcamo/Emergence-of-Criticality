%% KEEP PLOTS

% Entry point for CW criticality-signature plots using real-data targets.
% The combined implementation lives in signaturespapercurieweiss.m.

scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'signaturespapercurieweiss.m'));
