%% KEEP PLOTS

% Supplementary Curie-Weiss figure panels collected from
% signaturespapercurieweiss.m.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

%% Plot settings

files = ["Allenhldata.mat", "hippomuchidata.mat", "stringerhldata.mat"];
datasetNames = ["Allen", "Hippocampus", "Stringer"];
markers = {'square', 'o', '^'};
colors = ["#2676ad", "#4f4cc4", "#2f682c"];

realTargetSizes = [20 100 1000 -1];
Tzoom = linspace(0.8, 2.5, 250);
alphaVals = linspace(0, 2, 101);
cwTvals = linspace(0.5, 2, 250);
cwDhProbeVals = -logspace(-6, -2, 140);
cwPlotAbsJacobian = true;
cwPlotLogJacobian = true;

lineWidth = 1.5;
guideLineWidth = 1.2;
markerSize = 80;
axisLabelFontSize = 18;
tickLabelFontSize = 14;
titleFontSize = 14;
textFontName = 'Helvetica';

%% Build targets

targets = makeRealDataTargets(files, datasetNames, markers, colors, realTargetSizes);
nTargets = numel(targets);

for a = 1:nTargets
    N = targets(a).N;
    muTarget = targets(a).muTarget;
    chiTarget = targets(a).chiTarget;

    [hN, lambdaN] = hlambda(muTarget, chiTarget, N);

    targets(a).h = hN;
    targets(a).lambda = lambdaN;
end

Calpha = zeros(nTargets, numel(alphaVals));
JacAlpha = zeros(nTargets, numel(alphaVals));

for a = 1:nTargets
    N = targets(a).N;
    muTarget = targets(a).muTarget;
    lambdaN = targets(a).lambda;
    hGuess = targets(a).h;

    for b = 1:numel(alphaVals)
        alpha = alphaVals(b);
        lambdaAlpha = alpha * lambdaN;
        hAlpha = solveFieldForTargetMu(muTarget, lambdaAlpha, N, hGuess);

        [~, ~, ~, ~, C] = muChiExact2Spin(hAlpha, lambdaAlpha, N);
        [Jac, J11, J12, J21, J22] = Jacobian(hAlpha, lambdaAlpha, N);

        Calpha(a, b) = C / N;
        JacAlpha(a, b) = extractJacobianScalar(Jac, J11, J12, J21, J22);

        hGuess = hAlpha;
    end
end

responseResults = computeCWCriticalitySignatures(targets, cwTvals, cwDhProbeVals);

%% Supplementary multi-panel figure

figure('Name', 'Supplementary Curie-Weiss panels', ...
    'Color', 'w', ...
    'Position', [80, 80, 1350, 850]);
tiledlayout(3, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% Specific heat near T = 1

nexttile(1)
hold on

for a = 1:nTargets
    N = targets(a).N;
    hN = targets(a).h;
    lambdaN = targets(a).lambda;
    Cvals = zeros(size(Tzoom));

    for b = 1:numel(Tzoom)
        [~, ~, ~, ~, C] = muChiExact2Spin(hN / Tzoom(b), lambdaN / Tzoom(b), N);
        Cvals(b) = C;
    end

    plot(Tzoom, Cvals / N, ...
        'LineWidth', lineWidth, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xline(1, 'r-', 'T = 1', ...
    'LineWidth', guideLineWidth, ...
    'HandleVisibility', 'off');

xlabel('Dummy temperature T')
ylabel('Specific heat C(T)/N')
formatSupplementAxis(gca, textFontName, tickLabelFontSize, ...
    axisLabelFontSize, titleFontSize)

%% Specific heat at fixed mean activity

nexttile(2)
hold on

for a = 1:nTargets
    plot(alphaVals, Calpha(a, :), ...
        'LineWidth', lineWidth, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xline(1, 'r-', '\alpha = 1', ...
    'LineWidth', guideLineWidth, ...
    'HandleVisibility', 'off');

xlabel('Correlation scale \alpha')
ylabel('Specific heat C(\alpha)/N')
formatSupplementAxis(gca, textFontName, tickLabelFontSize, ...
    axisLabelFontSize, titleFontSize)

%% Jacobian versus dummy T

nexttile(3)
hold on

for a = 1:nTargets
    y = responseResults.jacMat(a, :);

    if cwPlotAbsJacobian
        y = abs(y);
    end

    if cwPlotLogJacobian
        y = max(y, realmin);
        semilogy(cwTvals, y, ...
            'LineWidth', lineWidth, ...
            'Color', targets(a).color, ...
            'DisplayName', targets(a).label);
    else
        plot(cwTvals, y, ...
            'LineWidth', lineWidth, ...
            'Color', targets(a).color, ...
            'DisplayName', targets(a).label);
    end
end

xline(1, 'r-', 'T = 1', ...
    'LineWidth', guideLineWidth, ...
    'HandleVisibility', 'off');

xlabel('Dummy temperature T')

if cwPlotAbsJacobian
    ylabel('|Jac|')
else
    ylabel('Jac')
end

formatSupplementAxis(gca, textFontName, tickLabelFontSize, ...
    axisLabelFontSize, titleFontSize)
yscale log

%% Jacobian at fixed mean activity

nexttile(4)
hold on

for a = 1:nTargets
    y = JacAlpha(a, :);

    if cwPlotAbsJacobian
        y = abs(y);
    end

    if cwPlotLogJacobian
        y = max(y, realmin);
        semilogy(alphaVals, y, ...
            'LineWidth', lineWidth, ...
            'Color', targets(a).color, ...
            'DisplayName', targets(a).label);
    else
        plot(alphaVals, y, ...
            'LineWidth', lineWidth, ...
            'Color', targets(a).color, ...
            'DisplayName', targets(a).label);
    end
end

xline(1, 'r-', '\alpha = 1', ...
    'LineWidth', guideLineWidth, ...
    'HandleVisibility', 'off');

xlabel('Correlation scale \alpha')

if cwPlotAbsJacobian
    ylabel('|Jac|')
else
    ylabel('Jac')
end

formatSupplementAxis(gca, textFontName, tickLabelFontSize, ...
    axisLabelFontSize, titleFontSize)

if cwPlotLogJacobian
    yscale log
end

%% Response around inferred T = 1 model

nexttile(5)
hold on

for a = 1:nTargets
    y = abs(responseResults.deltaMuResponse(a, :));
    y = max(y, realmin);

    loglog(cwDhProbeVals, y, ...
        'LineWidth', lineWidth, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

refTargetIndex = nTargets;
refHIndex = round(0.6 * numel(cwDhProbeVals));
refH = abs(cwDhProbeVals(refHIndex));
refMu = abs(responseResults.deltaMuResponse(refTargetIndex, refHIndex));
guide = refMu * (abs(cwDhProbeVals) / refH);

loglog(cwDhProbeVals, guide, 'k--', ...
    'LineWidth', guideLineWidth, ...
    'DisplayName', '\delta h^{-1} guide');

xlabel('Field perturbation \delta h')
ylabel('Response |\Delta\mu|')

formatSupplementAxis(gca, textFontName, tickLabelFontSize, ...
    axisLabelFontSize, titleFontSize)
xscale log
yscale log

%% Entropy versus energy for varying system sizes

nexttile(6)
hold on

for a = 1:nTargets
    N = targets(a).N;
    hN = targets(a).h;
    lambdaN = targets(a).lambda;

    [ePerSpin, sPerSpin] = cwEntropyEnergyCurve(hN, lambdaN, N);

    plot(ePerSpin, sPerSpin, ...
        'LineWidth', lineWidth, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xlimCurrent = xlim;
xGuide = linspace(xlimCurrent(1), xlimCurrent(2), 100);
plot(xGuide, xGuide, 'k--', ...
    'LineWidth', guideLineWidth, ...
    'DisplayName', 'S = E guide');

xlabel('Energy per spin E/N')
ylabel('Entropy per spin S/N')
formatSupplementAxis(gca, textFontName, tickLabelFontSize, ...
    axisLabelFontSize, titleFontSize)

set(gcf, 'Renderer', 'painters')

%% Local functions

function results = computeCWCriticalitySignatures(targets, Tvals, dhProbeVals)
    nTargets = numel(targets);
    nT = numel(Tvals);
    nH = numel(dhProbeVals);

    results.muMat = zeros(nTargets, nT);
    results.chiMat = zeros(nTargets, nT);
    results.jacMat = zeros(nTargets, nT);
    results.muResponse = zeros(nTargets, nH);
    results.deltaMuResponse = zeros(nTargets, nH);

    for a = 1:nTargets
        N = targets(a).N;
        hN = targets(a).h;
        lambdaN = targets(a).lambda;

        for b = 1:nT
            T = Tvals(b);

            hEff = hN / T;
            lambdaEff = lambdaN / T;

            [mu, chi] = muChiExact2Spin(hEff, lambdaEff, N);
            [Jac, J11, J12, J21, J22] = Jacobian(hEff, lambdaEff, N);

            results.muMat(a, b) = mu;
            results.chiMat(a, b) = chi;
            results.jacMat(a, b) = extractJacobianScalar(Jac, J11, J12, J21, J22);
        end

        [mu0, ~] = muChiExact2Spin(hN, lambdaN, N);

        for b = 1:nH
            hEff = hN + dhProbeVals(b);
            [muPerturbed, ~] = muChiExact2Spin(hEff, lambdaN, N);

            results.muResponse(a, b) = muPerturbed;
            results.deltaMuResponse(a, b) = muPerturbed - mu0;
        end
    end
end

function jacScalar = extractJacobianScalar(Jac, J11, J12, J21, J22)
    if isscalar(Jac)
        jacScalar = Jac;

    elseif isequal(size(Jac), [2, 2])
        jacScalar = det(Jac);

    else
        jacScalar = J11 * J22 - J12 * J21;
    end
end

function targets = makeRealDataTargets(files, datasetNames, markers, colors, requestedSizes)
    targets = struct([]);
    counter = 0;

    for i = 1:numel(files)
        obj = load(files(i));

        if ~isfield(obj, 'Data')
            error('File %s does not contain variable Data.', files(i));
        end

        Data = obj.Data;
        Data(Data == 0) = nan;

        Nss  = squeeze(mean(squeeze(Data(3,:,:,:)), 1, 'omitnan'));
        mus  = squeeze(mean(squeeze(Data(1,:,:,:)), 1, 'omitnan'));
        chis = squeeze(mean(squeeze(Data(2,:,:,:)), 1, 'omitnan'));

        baseColor = hex2rgb(colors(i));

        for v = 1:numel(requestedSizes)
            val = requestedSizes(v);

            if val == -1
                muThis = lastFinite(mus);
                chiThis = lastFinite(chis);
                NThis = lastFinite(Nss);
                labelThis = sprintf('%s full', datasetNames(i));
                colorThis = baseColor;
            else
                [muThis, chiThis, NThis] = getNearestSizePoint(mus, chis, Nss, val);
                labelThis = sprintf('%s N \\approx %d', datasetNames(i), val);

                if val == 20
                    colorThis = lightenColor(baseColor, 0.60);
                elseif val == 100
                    colorThis = lightenColor(baseColor, 0.45);
                elseif val == 1000
                    colorThis = lightenColor(baseColor, 0.19);
                else
                    colorThis = baseColor;
                end
            end

            if ~isfinite(muThis) || ~isfinite(chiThis) || ~isfinite(NThis)
                warning('Skipping invalid target from %s at requested size %g.', files(i), val);
                continue;
            end

            counter = counter + 1;
            targets(counter).label = labelThis;
            targets(counter).dataset = datasetNames(i);
            targets(counter).N = round(NThis);
            targets(counter).Nraw = NThis;
            targets(counter).muTarget = muThis;
            targets(counter).chiTarget = chiThis;
            targets(counter).color = colorThis;
            targets(counter).marker = markers{i};
        end
    end
end

function [muThis, chiThis, NThis] = getNearestSizePoint(mus, chis, Nss, targetN)
    valid = isfinite(mus) & isfinite(chis) & isfinite(Nss);

    if ~any(valid(:))
        muThis = nan;
        chiThis = nan;
        NThis = nan;
        return;
    end

    distance = abs(Nss - targetN);
    distance(~valid) = inf;

    [~, idx] = min(distance(:));

    muThis = mus(idx);
    chiThis = chis(idx);
    NThis = Nss(idx);
end

function val = lastFinite(x)
    x = x(:);
    idx = find(isfinite(x), 1, 'last');

    if isempty(idx)
        val = nan;
    else
        val = x(idx);
    end
end

function [ePerSpin, sPerSpin] = cwEntropyEnergyCurve(h, lambda, N)
    k = (0:N)';
    M = 2 * k - N;

    logOmega = gammaln(N + 1) - gammaln(k + 1) - gammaln(N - k + 1);
    E = -h * M - (lambda / (2 * N)) * M.^2;
    logZ = logsumexpStable(logOmega - E);

    [Euniq, logOmegaUniq] = combineDegenerateEnergies(E, logOmega);

    ePerSpin = (Euniq + logZ) / N;
    sPerSpin = logOmegaUniq / N;

    [ePerSpin, order] = sort(ePerSpin);
    sPerSpin = sPerSpin(order);
end

function [Euniq, logOmegaUniq] = combineDegenerateEnergies(E, logOmega)
    E = E(:);
    logOmega = logOmega(:);

    scale = max(1, max(abs(E)));
    tol = 1e-10 * scale;
    keys = round(E / tol);

    [uniqueKeys, ~, idx] = unique(keys);

    Euniq = zeros(numel(uniqueKeys), 1);
    logOmegaUniq = zeros(numel(uniqueKeys), 1);

    for a = 1:numel(uniqueKeys)
        mask = idx == a;
        Euniq(a) = mean(E(mask));
        logOmegaUniq(a) = logsumexpStable(logOmega(mask));
    end
end

function hOut = solveFieldForTargetMu(muTarget, lambda, N, hGuess)
    f = @(h) muFromExactModel(h, lambda, N) - muTarget;

    if nargin < 4 || isempty(hGuess) || ~isfinite(hGuess)
        hGuess = atanh(max(min(muTarget, 0.999999), -0.999999));
    end

    if abs(f(hGuess)) < 1e-11
        hOut = hGuess;
        return;
    end

    width = 1.0;
    hLow = hGuess - width;
    hHigh = hGuess + width;

    fLow = f(hLow);
    fHigh = f(hHigh);
    expandCount = 0;

    while sign(fLow) == sign(fHigh) && expandCount < 60
        width = 2 * width;
        hLow = hGuess - width;
        hHigh = hGuess + width;

        fLow = f(hLow);
        fHigh = f(hHigh);
        expandCount = expandCount + 1;
    end

    if sign(fLow) ~= sign(fHigh)
        hOut = fzero(f, [hLow, hHigh]);
    else
        hOut = fzero(f, hGuess);
    end
end

function mu = muFromExactModel(h, lambda, N)
    [mu, ~] = muChiExact2Spin(h, lambda, N);
end

function y = logsumexpStable(x)
    x = x(:);
    xmax = max(x);
    y = xmax + log(sum(exp(x - xmax)));
end

function rgb = hex2rgb(hex)
    hex = char(hex);
    hex = strrep(hex, '#', '');

    r = hex2dec(hex(1:2));
    g = hex2dec(hex(3:4));
    b = hex2dec(hex(5:6));

    rgb = [r g b] / 255;
end

function c = lightenColor(baseColor, amount)
    c = baseColor + amount * (1 - baseColor);
    c = min(max(c, 0), 1);
end

function formatSupplementAxis(ax, textFontName, tickLabelFontSize, axisLabelFontSize, titleFontSize)
    box(ax, 'on')
    grid(ax, 'off')
    axis(ax, 'square')
    ax.FontName = textFontName;
    ax.FontSize = tickLabelFontSize;
    ax.LineWidth = 1.2;
    ax.TickDir = 'both';

    ax.XLabel.FontName = textFontName;
    ax.XLabel.FontSize = axisLabelFontSize;
    ax.YLabel.FontName = textFontName;
    ax.YLabel.FontSize = axisLabelFontSize;
    ax.Title.FontName = textFontName;
    ax.Title.FontSize = titleFontSize;
end
