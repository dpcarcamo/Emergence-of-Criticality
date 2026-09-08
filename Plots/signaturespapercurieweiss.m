%% KEEP PLOTS

%% reproduce_TkacikMora_CW_plots.m
%
% Supplementary Curie-Weiss thermodynamic plots.
%
% Curie-Weiss analogues of plots from:
% Tkačik et al., "Thermodynamics and signatures of criticality in a network of neurons"
%
% Requires:
%
%   [h, lambda] = hlambda(muTarget, chiTarget, N)
%
% Optional check uses:
%
%   [mu, chi] = muChiExact2Spin(h, lambda, N)
%
% Curie-Weiss model:
%
%   P(s) proportional to exp[h sum_i s_i + lambda/(2N) (sum_i s_i)^2]
%
% Equivalently:
%
%   P(s) proportional to exp[-E(s)]
%
% with
%
%   E(s) = -h M - lambda/(2N) M^2
%   M = sum_i s_i
%
% Dummy temperature:
%
%   P_T(s) proportional to exp[-E(s)/T]
%
% Equivalent parameter scaling: h -> h/T and lambda -> lambda/T.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, 'Allen'))
addpath(fullfile(repoRoot, 'Hippo'))
addpath(fullfile(repoRoot, 'Stringer'))

%% reproduce_TkacikMora_CW_plots_with_real_targets.m
%
% Curie-Weiss analogues of Tkačik/Mora-style thermodynamic plots.
%
% Requires:
%   [h, lambda] = hlambda(muTarget, chiTarget, N)
%
% Optional:
%   [mu, chi] = muChiExact2Spin(h, lambda, N)
%
% Real data format:
%   obj = load(file)
%   Data = obj.Data
%   Data(1,:,:,:) = mu
%   Data(2,:,:,:) = chi
%   Data(3,:,:,:) = N
%
% Each real-data target becomes:
%   muTarget  = real data mu
%   chiTarget = real data chi
%   N         = real data N
%
% Then:
%   [h_N, lambda_N] = hlambda(muTarget, chiTarget, N)

clear; clc; close all;

%% ---------------- TARGET SOURCE ----------------

% Choose:
%   "manual"
%   "realData"
targetSource = "realData";

%% ---------------- MANUAL TARGET OPTIONS ----------------

manualNs = [20, 40, 60, 80, 100, 120];

manualMuTarget = 0.0;

% Option A: critical-like finite-size target
manualChiTargetFcn = @(N) sqrt(N);

% Option B: fixed chi target
% manualChiTargetFcn = @(N) 10;

%% ---------------- REAL DATA TARGET OPTIONS ----------------

files = ["Allenhldata.mat", "hippomuchidata.mat", "stringerhldata.mat"];

datasetNames = ["Allen", "Hippocampus", "Stringer"];

markers = {'square', 'o', '^'};

colors = ["#1f4fd6", "#c65100", "#167a43"];

% Which real-data points to use.
% Use -1 for the full/final point.
% Examples:
%   realTargetSizes = -1;
%   realTargetSizes = [20 100 1000 -1];
realTargetSizes = [20 100 1000 -1];

%% ---------------- PLOT OPTIONS ----------------

Tvals = linspace(0.5, 4.0, 250);
Tzoom = linspace(0.8, 2.5, 250);

alphaVals = linspace(0, 2, 101);
alphaExamples = [0, 1, 2];

plotLargestNOnlyForSEAlpha = true;

%% ---------------- CW CRITICALITY SIGNATURE OPTIONS ----------------

plotCWCriticalityByDataset = true;

cwTvals = linspace(0.5, 2, 250);
cwDhProbeVals = -logspace(-6, -2, 140);

cwPlotAbsJacobian = true;
cwPlotLogJacobian = true;

%% ---------------- BUILD TARGET LIST ----------------

switch targetSource

    case "manual"

        targets = makeManualTargets(manualNs, manualMuTarget, manualChiTargetFcn);

    case "realData"

        targets = makeRealDataTargets(files, datasetNames, markers, colors, realTargetSizes);

    otherwise

        error('targetSource must be "manual" or "realData".');
end

nTargets = numel(targets);

fprintf('\nNumber of targets = %d\n', nTargets);
%% ---------------- INDEPENDENT MODEL REFERENCES ----------------

addIndependentRefs = targetSource == "realData";

% One independent curve per dataset, using the full-dataset mu.
% N is also taken from the full-dataset point.
if addIndependentRefs
    independentRefs = makeIndependentRefsFromRealData(files, datasetNames, colors);
else
    independentRefs = struct([]);
end


%% ---------------- INFER BASELINE h_N AND lambda_N ----------------

fprintf('\nBaseline inferred models at T = 1\n');
fprintf('%18s %8s %14s %14s %14s %14s %14s\n', ...
    'label', 'N', 'muTarget', 'chiTarget', 'h_N', 'lambda_N', 'check chi');

for a = 1:nTargets

    N = targets(a).N;
    muTarget = targets(a).muTarget;
    chiTarget = targets(a).chiTarget;

    [hN, lambdaN] = hlambda(muTarget, chiTarget, N);

    targets(a).h = hN;
    targets(a).lambda = lambdaN;

    if exist('muChiExact2Spin', 'file') == 2
        [muCheck, chiCheck] = muChiExact2Spin(hN, lambdaN, N);
    else
        stats = cwCanonicalStats(hN, lambdaN, N, 1);
        muCheck = stats.mu;
        chiCheck = stats.chi;
    end

    targets(a).muCheck = muCheck;
    targets(a).chiCheck = chiCheck;

    fprintf('%18s %8d %14.6g %14.6g %14.6g %14.6g %14.6g\n', ...
        targets(a).label, N, muTarget, chiTarget, hN, lambdaN, chiCheck);
end

%% ---------------- FIGURE 1: ENTROPY-ENERGY AND HEAT CAPACITY ----------------

figure('Color', 'w', 'Position', [80, 100, 1550, 430]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Panel 1: entropy versus energy

nexttile;
hold on;

for a = 1:nTargets

    N = targets(a).N;
    hN = targets(a).h;
    lambdaN = targets(a).lambda;

    [ePerSpin, sPerSpin] = cwEntropyEnergyCurve(hN, lambdaN, N);

    plot(ePerSpin, sPerSpin, ...
        'LineWidth', 1.5, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xlimCurrent = xlim;
xGuide = linspace(xlimCurrent(1), xlimCurrent(2), 100);
plot(xGuide, xGuide, 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'S = E guide');

xlabel('energy per spin, E/N');
ylabel('entropy per spin, S/N');
title('Curie-Weiss entropy versus energy');
%legend('Location', 'best');
grid on;
box on;
axis square
%% Panel 2: heat capacity C(T)

nexttile;
hold on;

for a = 1:nTargets

    N = targets(a).N;
    hN = targets(a).h;
    lambdaN = targets(a).lambda;

    Cvals = zeros(size(Tvals));

    for b = 1:numel(Tvals)
        stats = cwCanonicalStats(hN, lambdaN, N, Tvals(b));
        Cvals(b) = stats.C;
    end

    plot(Tvals, Cvals, ...
        'LineWidth', 1.5, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xline(1, 'r-', 'T = 1', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

xlabel('dummy T');
ylabel('C(T)');
title('Heat capacity');
%legend('Location', 'best');
grid on;
box on;
axis square
% Independent-model reference curves
if addIndependentRefs
    for r = 1:numel(independentRefs)

        Nind = independentRefs(r).N;
        hind = independentRefs(r).h;
        lambdaInd = 0;

        Cind = zeros(size(Tvals));

        for b = 1:numel(Tvals)
            stats = cwCanonicalStats(hind, lambdaInd, Nind, Tvals(b));
            Cind(b) = stats.C;
        end

        plot(Tvals, Cind, '--', ...
            'LineWidth', 1.8, ...
            'Color', independentRefs(r).color, ...
            'DisplayName', sprintf('%s independent', independentRefs(r).dataset));
    end
end

%% Panel 3: specific heat C(T)/N

nexttile;
hold on;

for a = 1:nTargets

    N = targets(a).N;
    hN = targets(a).h;
    lambdaN = targets(a).lambda;

    Cvals = zeros(size(Tzoom));

    for b = 1:numel(Tzoom)
        stats = cwCanonicalStats(hN, lambdaN, N, Tzoom(b));
        Cvals(b) = stats.C;
    end

    plot(Tzoom, Cvals / N, ...
        'LineWidth', 1.5, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xline(1, 'r-', 'T = 1', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

xlabel('dummy T');
ylabel('C(T)/N');
title('Specific heat near T = 1');
%legend('Location', 'best');
grid on;
box on;
axis square
sgtitle(sprintf('Curie-Weiss thermodynamic plots, targetSource = %s', targetSource));


% Independent-model reference curves
if addIndependentRefs
    for r = 1:numel(independentRefs)

        Nind = independentRefs(r).N;
        hind = independentRefs(r).h;
        lambdaInd = 0;

        Cind = zeros(size(Tzoom));

        for b = 1:numel(Tzoom)
            stats = cwCanonicalStats(hind, lambdaInd, Nind, Tzoom(b));
            Cind(b) = stats.C;
        end

        plot(Tzoom, Cind / Nind, '--', ...
            'LineWidth', 1.8, ...
            'Color', independentRefs(r).color, ...
            'DisplayName', sprintf('%s independent', independentRefs(r).dataset));
    end
end

%% ---------------- FIGURE 2: ALPHA ENSEMBLE AT FIXED MEAN ----------------

% In the paper, alpha changes correlation strength while fields are adjusted
% to keep spike rates fixed.
%
% For Curie-Weiss:
%
%   lambda_alpha = alpha * lambda_N
%
% with h_alpha solved so that:
%
%   mu(h_alpha, lambda_alpha, N) = muTarget.

Salpha = zeros(nTargets, numel(alphaVals));
Calpha = zeros(nTargets, numel(alphaVals));
halphaMat = zeros(nTargets, numel(alphaVals));

fprintf('\nComputing alpha ensembles at fixed mu...\n');

for a = 1:nTargets

    N = targets(a).N;
    muTarget = targets(a).muTarget;
    lambdaN = targets(a).lambda;
    hGuess = targets(a).h;

    fprintf('  %s, N = %d\n', targets(a).label, N);

    for b = 1:numel(alphaVals)

        alpha = alphaVals(b);

        lambdaAlpha = alpha * lambdaN;
        hAlpha = solveFieldForTargetMu(muTarget, lambdaAlpha, N, hGuess);

        stats = cwCanonicalStats(hAlpha, lambdaAlpha, N, 1);

        Salpha(a, b) = stats.S / N;
        Calpha(a, b) = stats.C / N;
        halphaMat(a, b) = hAlpha;

        hGuess = hAlpha;
    end
end

figure('Color', 'w', 'Position', [120, 140, 1550, 430]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Panel 1: entropy versus energy for alpha examples

nexttile;
hold on;

if plotLargestNOnlyForSEAlpha
    plotTargetList = nTargets;
else
    plotTargetList = 1:nTargets;
end

for a = plotTargetList

    N = targets(a).N;
    muTarget = targets(a).muTarget;
    lambdaN = targets(a).lambda;

    for q = 1:numel(alphaExamples)

        alpha = alphaExamples(q);
        lambdaAlpha = alpha * lambdaN;

        hAlpha = solveFieldForTargetMu(muTarget, lambdaAlpha, N, targets(a).h);

        [ePerSpin, sPerSpin] = cwEntropyEnergyCurve(hAlpha, lambdaAlpha, N);

        plot(ePerSpin, sPerSpin, ...
            'LineWidth', 1.7, ...
            'DisplayName', sprintf('%s, \\alpha = %.1f', targets(a).label, alpha));
    end
end

xlimCurrent = xlim;
xGuide = linspace(xlimCurrent(1), xlimCurrent(2), 100);
plot(xGuide, xGuide, 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'S = E guide');

xlabel('energy per spin, E/N');
ylabel('entropy per spin, S/N');
title('Entropy versus energy for \alpha ensemble');
legend('Location', 'best');
grid on;
box on;
axis square

%% Panel 2: entropy per spin versus alpha

nexttile;
hold on;

for a = 1:nTargets
    plot(alphaVals, Salpha(a, :), ...
        'LineWidth', 1.5, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xline(1, 'r-', '\alpha = 1', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

xlabel('\alpha');
ylabel('S(\alpha)/N');
title('Entropy at fixed activity');
%legend('Location', 'best');
grid on;
box on;
axis square
%% Panel 3: heat capacity per spin versus alpha

nexttile;
hold on;

for a = 1:nTargets
    plot(alphaVals, Calpha(a, :), ...
        'LineWidth', 1.5, ...
        'Color', targets(a).color, ...
        'DisplayName', targets(a).label);
end

xline(1, 'r-', '\alpha = 1', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

xlabel('\alpha');
ylabel('C(\alpha)/N');
title('Heat capacity at fixed activity');
%legend('Location', 'best');
grid on;
box on;
axis square
sgtitle('Curie-Weiss analogue of changing correlations at fixed activity');

%% ---------------- FIGURE 3: TARGETS IN MU-CHI SPACE ----------------

figure('Color', 'w', 'Position', [180, 180, 520, 460]);
hold on;

for a = 1:nTargets
    scatter(targets(a).muTarget, targets(a).chiTarget, 120, ...
        'Marker', targets(a).marker, ...
        'MarkerEdgeColor', targets(a).color, ...
        'LineWidth', 1.5, ...
        'DisplayName', targets(a).label);
end

x = linspace(-1, 1, 1000);
plot(x, 1 - x.^2, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'independent bound');
plot(x, x.*(1-x.^2)./(x-atanh(x).*(1-x.^2)), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'independent bound');
xlim([-1,-0.5])
xlabel('m');
ylabel('\chi');
title('Targets used for Curie-Weiss inversion');
set(gca, 'YScale', 'log');
grid on;
box on;
axis square;
%legend('Location', 'best');

%% ---------------- FIGURE 4+: CW CRITICALITY SIGNATURES BY DATASET ----------------

if plotCWCriticalityByDataset
    datasets = unique([targets.dataset], 'stable');

    for d = 1:numel(datasets)

        datasetName = datasets(d);
        datasetTargets = targets([targets.dataset] == datasetName);

        if isempty(datasetTargets)
            continue;
        end

        fprintf('\nComputing CW criticality signatures for %s\n', datasetName);

        results = computeCWCriticalitySignatures(datasetTargets, cwTvals, cwDhProbeVals);
        nDatasetTargets = numel(datasetTargets);

        figure('Color', 'w', 'Position', [100, 100, 1550, 430]);
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        hold on;

        for a = 1:nDatasetTargets
            plot(cwTvals, results.chiMat(a, :), ...
                'LineWidth', 1.7, ...
                'Color', datasetTargets(a).color, ...
                'DisplayName', targetDisplayName(datasetTargets(a)));
        end

        xline(1, 'k--', 'T = 1', 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');

        xlabel('dummy T');
        ylabel('\chi');
        title('\chi versus dummy T');
        legend('Location', 'best');
        grid on;
        box on;

        nexttile;
        hold on;

        for a = 1:nDatasetTargets

            y = results.jacMat(a, :);

            if cwPlotAbsJacobian
                y = abs(y);
            end

            if cwPlotLogJacobian
                y = max(y, realmin);
                semilogy(cwTvals, y, ...
                    'LineWidth', 1.7, ...
                    'Color', datasetTargets(a).color, ...
                    'DisplayName', targetDisplayName(datasetTargets(a)));
            else
                plot(cwTvals, y, ...
                    'LineWidth', 1.7, ...
                    'Color', datasetTargets(a).color, ...
                    'DisplayName', targetDisplayName(datasetTargets(a)));
            end
        end

        xline(1, 'k--', 'T = 1', 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');

        xlabel('dummy T');

        if cwPlotAbsJacobian
            ylabel('|Jac|');
        else
            ylabel('Jac');
        end

        title('Jacobian versus dummy T');
        legend('Location', 'best');
        grid on;
        box on;

        nexttile;
        hold on;

        for a = 1:nDatasetTargets

            y = abs(results.deltaMuResponse(a, :));
            y = max(y, realmin);

            loglog(cwDhProbeVals, y, ...
                'LineWidth', 1.7, ...
                'Color', datasetTargets(a).color, ...
                'DisplayName', targetDisplayName(datasetTargets(a)));
        end

        refTargetIndex = nDatasetTargets;
        refHIndex = round(0.6 * numel(cwDhProbeVals));

        refH = abs(cwDhProbeVals(refHIndex));
        refMu = abs(results.deltaMuResponse(refTargetIndex, refHIndex));
        
        guide = refMu * (abs(cwDhProbeVals) / refH).^(1);

        loglog(cwDhProbeVals, guide, 'k--', 'LineWidth', 1.4, ...
            'DisplayName', '\sim |\delta h^{-1}|');

        xlabel('field perturbation \delta h');
        ylabel('|\Delta m|');
        title('Response around inferred T = 1 model');
        legend('Location', 'best');
        grid on;
        box on;
        xscale log;
        yscale log;

        sgtitle(sprintf('Curie-Weiss criticality signatures: %s data targets', datasetName));
    end
end

%% ---------------- LOCAL FUNCTIONS ----------------

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

        fprintf('  %s, N = %d\n', targets(a).label, N);

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

            dh = dhProbeVals(b);

            hEff = hN + dh;
            lambdaEff = lambdaN;

            [muPerturbed, ~] = muChiExact2Spin(hEff, lambdaEff, N);

            results.muResponse(a, b) = muPerturbed;
            results.deltaMuResponse(a, b) = muPerturbed - mu0;
        end
    end
end

function label = targetDisplayName(target)

    label = sprintf('N = %d', target.N);
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

function targets = makeManualTargets(Ns, muTarget, chiTargetFcn)

    targets = struct([]);

    for a = 1:numel(Ns)

        N = Ns(a);
        chiTarget = chiTargetFcn(N);

        targets(a).label = sprintf('N = %d', N);
        targets(a).dataset = "manual";
        targets(a).N = round(N);
        targets(a).Nraw = N;
        targets(a).muTarget = muTarget;
        targets(a).chiTarget = chiTarget;
        targets(a).color = linesColor(a);
        targets(a).marker = 'o';
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

function stats = cwCanonicalStats(h, lambda, N, T)
    % Exact canonical statistics for Curie-Weiss at dummy temperature T.
    %
    % Energy:
    %
    %   E(k) = -h M - lambda/(2N) M^2
    %
    % Temperature deformation:
    %
    %   P_T(k) proportional to binom(N,k) exp[-E(k)/T]
    %
    % Heat capacity:
    %
    %   C(T) = Var_T(E) / T^2

    k = (0:N)';
    M = 2 * k - N;
    m = M / N;

    logOmega = gammaln(N + 1) - gammaln(k + 1) - gammaln(N - k + 1);

    E = -h * M - (lambda / (2 * N)) * M.^2;

    logWeight = logOmega - E / T;
    logZ = logsumexpStable(logWeight);

    p = exp(logWeight - logZ);

    meanM = sum(p .* M);
    meanm = meanM / N;

    meanm2 = sum(p .* m.^2);
    chi = N * (meanm2 - meanm^2);

    meanE = sum(p .* E);
    meanE2 = sum(p .* E.^2);
    varE = meanE2 - meanE^2;

    C = varE / T^2;

    S = logZ + meanE / T;

    stats.mu = meanm;
    stats.chi = chi;
    stats.E = meanE;
    stats.varE = varE;
    stats.C = C;
    stats.S = S;
    stats.logZ = logZ;
    stats.pK = p;
    stats.k = k;
    stats.m = m;
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

    f = @(h) cwCanonicalStats(h, lambda, N, 1).mu - muTarget;

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
    % amount = 0 gives baseColor.
    % amount = 1 gives white.

    c = baseColor + amount * (1 - baseColor);
    c = min(max(c, 0), 1);
end

function c = linesColor(i)

    L = lines(max(i, 7));
    c = L(i, :);
end

function independentRefs = makeIndependentRefsFromRealData(files, datasetNames, colors)
    % Build one independent-model reference per dataset.
    %
    % Uses the full/final dataset mean activity:
    %
    %   h_ind = atanh(mu_full)
    %   lambda_ind = 0
    %
    % For N, this uses the full/final dataset N.

    independentRefs = struct([]);

    for i = 1:numel(files)

        obj = load(files(i));

        if ~isfield(obj, 'Data')
            error('File %s does not contain variable Data.', files(i));
        end

        Data = obj.Data;
        Data(Data == 0) = nan;

        Nss  = squeeze(mean(squeeze(Data(3,:,:,:)), 1, 'omitnan'));
        mus  = squeeze(mean(squeeze(Data(1,:,:,:)), 1, 'omitnan'));

        % Dataset subsampling axis:
        % most datasets use dim = 1, hippocampus uses dim = 2.
        subSampleDim = 1;
        if i == 2
            subSampleDim = 2;
        end

        I = finalSubsampleIndices(mus, subSampleDim);

        muFullVals = mus(I);
        NFullVals = Nss(I);

        muFull = mean(muFullVals, 'omitnan');
        NFull = round(mean(NFullVals, 'omitnan'));

        % Avoid infinite atanh if mu is numerically too close to +/- 1.
        muClipped = min(max(muFull, -0.999999999), 0.999999999);

        hInd = atanh(muClipped);

        independentRefs(i).dataset = datasetNames(i);
        independentRefs(i).mu = muFull;
        independentRefs(i).N = NFull;
        independentRefs(i).h = hInd;
        independentRefs(i).lambda = 0;
        independentRefs(i).color = hex2rgb(colors(i));
    end
end

function I = finalSubsampleIndices(A, dim)
    % Return linear indices for the final subsample along dimension dim.

    sz = size(A);

    if isvector(A)
        valid = find(isfinite(A(:)));
        if isempty(valid)
            I = [];
        else
            I = valid(end);
        end
        return;
    end

    if dim == 1

        lastRow = sz(1);
        cols = 1:sz(2);
        I = sub2ind(sz, lastRow * ones(size(cols)), cols);

    elseif dim == 2

        rows = 1:sz(1);
        lastCol = sz(2);
        I = sub2ind(sz, rows, lastCol * ones(size(rows)));

    else

        error('Only dim = 1 or dim = 2 is supported here.');
    end

    I = I(:);
end
