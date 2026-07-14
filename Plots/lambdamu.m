%% KEEP PLOTS
% File to plot the double well theory
% Figure mapping: this script generates the Fig. 5c-d comparison between
% exact finite-size inference and the double-well approximation. It also
% opens auxiliary h/lambda diagnostic figures.

plotsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(plotsDir);
addpath(repoRoot)
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))
%% Plot one population from each dataset: refit h/lambda vs theory

files = ["Allenhldata.mat", "hippomuchidata.mat", "stringerhldata.mat"];

dataMarker = 'o';
colors = ["#2676ad", "#4f4cc4", "#2f682c"];

markerSize = 70;
markerEdgeColor = [0 0 0];
markerLineWidth = 0.5;
markerFaceAlpha = 0.75;
fitLineWidth = 1.5;
theoryLineWidth = 1.5;
axisLabelFontSize = 30;
tickLabelFontSize = 18;
legendFontSize = 14;
textFontName = 'Helvetica';
labelFont = ['\fontname{' textFontName '}'];
xTickValues = 10.^(0:4);
xTickLabels = arrayfun(@(x) sprintf('10^{%d}', x), 0:4, 'UniformOutput', false);

figure(1); hold on
figure(3); hold on

for i = 1:length(files)

    obj = load(files(i));

    Data = obj.Data;
    Data(Data == 0) = nan;

    Nss  = squeeze(mean(squeeze(Data(3,:,:,:)), 1, 'omitnan'));
    mus  = squeeze(mean(squeeze(Data(1,:,:,:)), 1, 'omitnan'));
    chis = squeeze(mean(squeeze(Data(2,:,:,:)), 1, 'omitnan'));

    % Pick one population
    if i == 2
        % Hippocampus only has one population
        Nuse   = Nss(:);
        muse   = mus(:);
        chiuse = chis(:);
    else
        % Allen/Stringer: columns are populations, rows are changing N
        counts = sum(isfinite(Nss), 1);
        [~, popIdx] = max(counts);

        Nuse   = Nss(:, popIdx);
        muse   = mus(:, popIdx);
        chiuse = chis(:, popIdx);
    end

    % Clean selected data
    valid = isfinite(Nuse) & isfinite(muse) & isfinite(chiuse) ...
          & abs(muse) < 1 & chiuse > 0;

    Nuse   = Nuse(valid);
    muse   = muse(valid);
    chiuse = chiuse(valid);

    % Refit h and lambda using hlambda
    hFit = nan(size(Nuse));
    lFit = nan(size(Nuse));

    for j = 1:length(Nuse)
        N   = round(Nuse(j));
        mu  = muse(j);
        chi = chiuse(j);

        try
            [hFit(j), lFit(j)] = hlambda(mu, chi, N);
        catch
            hFit(j) = nan;
            lFit(j) = nan;
        end
    end

    % Theory from exactly the same selected data
    u0 = sqrt(chiuse ./ Nuse + muse.^2);

    hTheory = nan(size(muse));
    lTheory = nan(size(muse));

    goodTheory = isfinite(u0) ...
              & u0 > 0 ...
              & abs(u0) < 1 ...
              & abs(muse ./ u0) < 1;

    lTheory(goodTheory) = atanh(u0(goodTheory)) ./ u0(goodTheory);

    hTheory(goodTheory) = atanh(muse(goodTheory) ./ u0(goodTheory)) ...
                        ./ (Nuse(goodTheory) .* u0(goodTheory));

    % Keep only points where both fit and theory are valid
    good = isfinite(hFit) & isfinite(lFit) ...
         & isfinite(hTheory) & isfinite(lTheory);

    Nuse = Nuse(good);
    muse = muse(good);
    chiuse = chiuse(good);
    hFit = hFit(good);
    lFit = lFit(good);
    hTheory = hTheory(good);
    lTheory = lTheory(good);

    % Sort by N
    [Nuse, idx] = sort(Nuse);
    muse = muse(idx);
    chiuse = chiuse(idx);
    hFit = hFit(idx);
    lFit = lFit(idx);
    hTheory = hTheory(idx);
    lTheory = lTheory(idx);
    % h versus N
    figure(1)

    scatter(Nuse, hFit, markerSize, ...
        'Marker', dataMarker, ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', colors(i), ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)

    plot(Nuse, hTheory, ...
        'Marker', 'none', ...
        'Color', colors(i), ...
        'LineWidth', theoryLineWidth, ...
        'MarkerSize', markerSize*0.7)

    figure(3)

    scatter(Nuse, lFit, markerSize, ...
        'Marker', dataMarker, ...
        'MarkerEdgeColor', markerEdgeColor, ...
        'MarkerFaceColor', colors(i), ...
        'MarkerFaceAlpha', markerFaceAlpha, ...
        'LineWidth', markerLineWidth)

    plot(Nuse, lTheory, ...
        'Marker', 'none', ...
        'Color', colors(i), ...
        'LineWidth', theoryLineWidth, ...
        'MarkerSize', markerSize*0.7)

end

% Format h versus N
figure(1)
set(gca, 'XScale', 'log')
xlabel([labelFont 'Number of neurons {\itN}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'External field {\ith}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
axis square
box on
ax = gca;
ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;
ax.XTick = xTickValues;
ax.XTickLabel = xTickLabels;
ax.TickLabelInterpreter = 'tex';
set(gca, 'TickDir', 'both')

yscale log


figure(3)
xlabel([labelFont 'Number of neurons {\itN}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
ylabel([labelFont 'Interaction strength {\itJ}'], ...
    'Interpreter', 'tex', ...
    'FontSize', axisLabelFontSize)
axis square
box on
%yscale log
xscale log
ax = gca;
ax.FontName = textFontName;
ax.FontSize = tickLabelFontSize;
ax.XTick = xTickValues;
ax.XTickLabel = xTickLabels;
ax.TickLabelInterpreter = 'tex';
set(gca, 'TickDir', 'both')

ylim([0,2.5])
