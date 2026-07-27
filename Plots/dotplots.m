%% KEEP PLOTS

%% Box in (mu,chi) and image in (h,lambda)
% Figure mapping: this script generates the Fig. 3c grid and mapped box
% boundaries in statistics space and model space.
%% Precompute inverse map data
clear; clc; 

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot);

Ns = [20, 100, 1000];
axisLabelFontSize = 30;
tickLabelFontSize = 26;
titleFontSize = 18;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
plotModelPanelsSeparately = false;
insetFontScale = 1.8;
insetAxisLabelFontSize = round(axisLabelFontSize * insetFontScale);
insetTickLabelFontSize = round(tickLabelFontSize * insetFontScale);
forceRecompute = false;
cacheFile = fullfile(plotsDir, 'dotplots_cache.mat');
cacheVersion = 1;

mus  = linspace(-0.9,-0.02,24);
chis = linspace(1,3.5,24); 

edgePts = 50;

muMin = min(mus);
muMax = max(mus);
chiMin = min(chis);
chiMax = max(chis);

% Interior grid
hs = zeros(length(mus), length(chis), length(Ns));
ls = zeros(length(mus), length(chis), length(Ns));

% Box edges
muBottom  = linspace(muMin, muMax, edgePts);
chiBottom = chiMin * ones(size(muBottom));

chiRight = linspace(chiMin, chiMax, edgePts);
muRight  = muMax * ones(size(chiRight));

muTop  = linspace(muMax, muMin, edgePts);
chiTop = chiMax * ones(size(muTop));

chiLeft = linspace(chiMax, chiMin, edgePts);
muLeft  = muMin * ones(size(chiLeft));

muEdges  = {muBottom,  muRight,  muTop,  muLeft};
chiEdges = {chiBottom, chiRight, chiTop, chiLeft};

%% Load or compute inverse map data

useCache = false;

if ~forceRecompute && exist(cacheFile, 'file') == 2
    cached = load(cacheFile);
    useCache = isDotplotCacheValid(cached, Ns, mus, chis, edgePts, ...
        muEdges, chiEdges, cacheVersion);

    if useCache
        disp(['Loading dotplot data from ', cacheFile])
        hs = cached.hs;
        ls = cached.ls;
        hEdges = cached.hEdges;
        lEdges = cached.lEdges;
    else
        disp('Dotplot cache does not match current grid. Recomputing.')
    end
end

if ~useCache
    [hs, ls, hEdges, lEdges] = computeDotplotMap(Ns, mus, chis, muEdges, chiEdges);
    save(cacheFile, 'hs', 'ls', 'hEdges', 'lEdges', ...
        'Ns', 'mus', 'chis', 'edgePts', 'muEdges', 'chiEdges', ...
        'cacheVersion', '-v7.3')
    disp(['Saved dotplot data to ', cacheFile])
end

%%
edgeColors = {
    [0.05 0.25 1 0.55], ...
    [0 0.85 0.15 0.55], ...
    [1 0.05 0.05 0.55], ...
    [0 0.95 0.95 0.55]
};

dotSize = 45;
edgeWidth = 5;

muMin = min(mus);
muMax = max(mus);
chiMin = min(chis);
chiMax = max(chis);

if ~plotModelPanelsSeparately
    figure('Name', 'Mapped dot plots for fixed N', 'Color', 'w')
    modelLayout = tiledlayout(1, length(Ns), ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
end

for a = 1:length(Ns)

    N = Ns(a);

    % Original box in (mu,chi)
    if a == 1
    figure
    hold on
    set(gca, 'TickDir', 'both')
    box on

    for e = 1:4
        plot(muEdges{e}, chiEdges{e}, ...
            'Color', edgeColors{e}, ...
            'LineWidth', edgeWidth);
    end

    % Interior points only
    mus_in  = mus(2:end-1);
    chis_in = chis(2:end-1);

    [MU, CHI] = meshgrid(mus_in, chis_in);
    scatter(MU(:), CHI(:), dotSize, 'k', 'filled');

    xlabel([labelFont 'Average activity \langle\mu\rangle'], ...
        'Interpreter', 'tex', ...
        'FontSize', insetAxisLabelFontSize)
    ylabel([labelFont 'Susceptibility \chi'], ...
        'Interpreter', 'tex', ...
        'FontSize', insetAxisLabelFontSize)
    % title(['(\mu,\chi), N = ', num2str(N)], ...
    %     'FontName', textFontName, ...
    %     'FontSize', titleFontSize)

    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = insetTickLabelFontSize;

    axis square
    xlim([muMin-0.05, muMax+0.05])
    ylim([chiMin-0.1, chiMax+0.1])
    end

    % Transformed box in (h,lambda)
    if plotModelPanelsSeparately
        figure('Name', sprintf('Mapped dot plot N = %d', N), 'Color', 'w')
    else
        nexttile(modelLayout)
    end
    hold on
    

    hs_in = hs(2:end-1, 2:end-1, a);
    ls_in = ls(2:end-1, 2:end-1, a);

    scatter(hs_in(:), ls_in(:), dotSize, 'k', 'filled', MarkerFaceAlpha=0.6);

    for e = 1:4
        plot(hEdges{a,e}, lEdges{a,e}, ...
            'Color', edgeColors{e}, ...
            'LineWidth', edgeWidth);
    end

    plot([0, 0], [1, 3], 'r', 'LineWidth', 2)
    scatter(0, 1, 85, 'r', 'filled')

    xlabel([labelFont 'External field {\ith}'], ...
        'Interpreter', 'tex', ...
        'FontSize', axisLabelFontSize)
    if plotModelPanelsSeparately || a == 1
        ylabel([labelFont 'Interaction strength \lambda'], ...
            'Interpreter', 'tex', ...
            'FontSize', axisLabelFontSize)
    else
        ylabel('')
    end
    % title(['(h,\lambda), N = ', num2str(N)], ...
    %     'FontName', textFontName, ...
    %     'FontSize', titleFontSize)

    axis square

    % xscale log
    % xlim([-1,-0.0005])
    % xticks([-1, -0.1, -0.01, -0.001])
    % xticklabels({'-10^0','-10^{-1}','-10^{-2}','-10^{-3}' })

    ylim([0,3])
    set(gca, 'TickDir', 'both')
    ax = gca;
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    box on
end

function [hs, ls, hEdges, lEdges] = computeDotplotMap(Ns, mus, chis, muEdges, chiEdges)
    hs = zeros(length(mus), length(chis), length(Ns));
    ls = zeros(length(mus), length(chis), length(Ns));
    hEdges = cell(length(Ns),4);
    lEdges = cell(length(Ns),4);

    for a = 1:length(Ns)

        N = Ns(a);

        fprintf('Computing N = %d\n', N)

        % Map grid points
        for j = 1:length(mus)
            for k = 1:length(chis)

                fprintf('N index %d, mu index %d, chi index %d\n', a, j, k)

                mu = mus(j);
                chi = chis(k);

                [h, lambda] = hlambda(mu, chi, N);

                hs(j,k,a) = h;
                ls(j,k,a) = lambda;
            end
        end

        % Map box edges
        for e = 1:4
            hEdge = zeros(size(muEdges{e}));
            lEdge = zeros(size(muEdges{e}));

            for q = 1:length(muEdges{e})
                fprintf('N index %d, edge %d, point %d\n', a, e, q)
                [hTemp, lTemp] = hlambda(muEdges{e}(q), chiEdges{e}(q), N);
                hEdge(q) = hTemp;
                lEdge(q) = lTemp;
            end

            hEdges{a,e} = hEdge;
            lEdges{a,e} = lEdge;
        end
    end
end

function tf = isDotplotCacheValid(cached, Ns, mus, chis, edgePts, muEdges, chiEdges, cacheVersion)
    requiredFields = {'hs', 'ls', 'hEdges', 'lEdges', 'Ns', 'mus', ...
        'chis', 'edgePts', 'muEdges', 'chiEdges', 'cacheVersion'};

    tf = all(isfield(cached, requiredFields));

    if ~tf
        return
    end

    tf = cached.cacheVersion == cacheVersion ...
        && isequal(cached.Ns, Ns) ...
        && isequal(cached.mus, mus) ...
        && isequal(cached.chis, chis) ...
        && isequal(cached.edgePts, edgePts) ...
        && isequal(cached.muEdges, muEdges) ...
        && isequal(cached.chiEdges, chiEdges);
end
