%% fit_RB_Ising_Allen_NN_gradient_jacobian_specificheat.m
%
% Standalone script for fitting random-bond pairwise Ising models from
% Allen-style nearest-neighbor spatial clusters.
%
% Spins are in plus/minus notation:
%
%   s_i = +/- 1
%
% Model:
%
%   P(s) proportional to exp[ sum_i h_i s_i + sum_{i<j} J_ij s_i s_j ]
%
% Target sufficient statistics:
%
%   <s_i>
%   <s_i s_j>
%
% Fitting:
%
%   Uses the standard exponential-family gradient.
%
%   Objective:
%       L(theta) = log Z(theta) - theta' * T_data
%
%   Gradient:
%       grad L = <T>_model - T_data
%
%   Update:
%       theta <- theta - stepSize * grad
%
% There is no ridge regression or inverse-Jacobian update in the inference.
%
% Jacobian diagnostic:
%
%   After fitting, the script estimates
%
%       J_ab = d<T_a>/d theta_b = Cov(T_a,T_b)
%
%   along the fictitious-temperature path
%
%       h -> h/T
%       J -> J/T
%
%   and plots:
%
%       C(T) = Var_T(E)/T^2
%       log det Jacobian
%       (1/p) log det Jacobian
%
% The small diagonal jitter in estimateJacobianLogdet is only for numerical
% stabilization of the determinant estimate from finite MCMC samples. It is
% not used in inference.
%
% Energy convention for the specific heat:
%
%   E(s) = -sum_i h_i s_i - sum_{i<j} J_ij s_i s_j
%
% The fictitious temperature ensemble is:
%
%   P_T(s) proportional to exp[-E(s)/T]
%
% so
%
%   C(T) = Var_T(E)/T^2.

clear; clc; close all;

%% ---------------- PATH SETUP ----------------

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
dataFolder = fullfile(scriptDir, 'Allen_data_share');

addpath(repoRoot);
addpath(dataFolder);

listing = dir(dataFolder);
isMat = endsWith({listing.name}, ".mat");
listing = listing(isMat);

%% ---------------- RUN OPTIONS ----------------

datasetName = "Allen";

% Model sizes to fit.
fitNs = [10, 20, 40];

% To keep runtime reasonable, start with one recording and one cluster.
% Increase these after checking the script runs.
maxRecordingsToUse = 5;
nClustersPerNPerRecording = 5;

recordingIndices = 1:min(maxRecordingsToUse, length(listing));

% Center selection:
%   "random"
%   "highestVariance"
%   "first"
centerSelection = "random";



%% ---------------- FITTING OPTIONS ----------------

fitOpts.maxIter = 200;
fitOpts.tol = 2e-3;

% Gradient descent step size for negative log likelihood.
% Update is:
%
%   theta <- theta - stepSize * (<T>_model - T_data)
fitOpts.stepSize = 0.05;

fitOpts.displayEvery = 1;

% Initialization:
%   "independent"
%   "zero"
%   "naiveMeanField"
fitOpts.initMode = "independent";
fitOpts.maxAbsInitialJ = 0.25;

%% ---------------- MCMC OPTIONS FOR FITTING ----------------
% These samples are used only to estimate <T>_model during gradient descent.

mcFit.nChains = 12;
mcFit.nSweeps = 5000;
mcFit.burnSweeps = 1000;
mcFit.thinSweeps = 5;
mcFit.initMode = "random";

%% ---------------- TEMPERATURE JACOBIAN AND CHI OPTIONS ----------------

Tvals = linspace(0.5, 2.5, 50);

mcJac.nChains = 16;
mcJac.nSweeps = 3500;
mcJac.burnSweeps = 1000;
mcJac.thinSweeps = 5;
mcJac.initMode = "random";

% Small diagonal jitter for finite-sample log determinant estimation only.
% This does NOT affect fitting.
logdetJitterFrac = 1e-6;

% Choose:
%
%   "effectiveParams"
%       Jacobian is d stats / d(h/T,J/T) = Cov_T(T_a,T_b).
%
%   "baselineParams"
%       Jacobian is d stats / d(h,J) after applying h/T,J/T.
%       This equals (1/T) Cov_T, so logdet shifts by -p log T.
jacobianParameterization = "effectiveParams";

%% ---------------- BUILD CLUSTERS AND FIT MODELS ----------------

models = struct([]);
modelCounter = 0;

for rr = 1:numel(recordingIndices)

    namenum = recordingIndices(rr);

    fprintf("\n==============================\n");
    fprintf("Recording %d of %d\n", namenum, length(listing));
    fprintf("==============================\n");

    filename = listing(namenum).name;
    fullFileName = fullfile(dataFolder, filename);

    rec = loadAllenRecordingMoments(fullFileName);

    datamean = rec.datamean;
    datacorr = rec.datacorr;
    Locations = rec.Locations;
    num_neurons = rec.num_neurons;

    distanceSMetric = squareform(pdist(Locations));

    for nIdx = 1:numel(fitNs)

        N = fitNs(nIdx);

        if num_neurons < N
            warning("Skipping recording %d, N=%d because only %d neurons.", ...
                namenum, N, num_neurons);
            continue;
        end

        centers = chooseClusterCenters(datamean, nClustersPerNPerRecording, centerSelection);

        for cc = 1:numel(centers)

            centerSpin = centers(cc);

            region = nearestNeighborRegion(distanceSMetric, centerSpin, N);

            reducedMean = datamean(region);
            reducedCorr = datacorr(region, region);

            [pairI, pairJ] = pairList(N);
            qTarget = reducedCorr(sub2ind([N, N], pairI, pairJ));

            targetStats = [reducedMean; qTarget];

            modelCounter = modelCounter + 1;

            fprintf("\n------------------------------\n");
            fprintf("Fitting model %d\n", modelCounter);
            fprintf("Recording %d, cluster %d, N = %d\n", namenum, cc, N);
            fprintf("------------------------------\n");

            [h, J, fitInfo] = fitIsingByGradientMatching( ...
                targetStats, N, pairI, pairJ, fitOpts, mcFit);

            models(modelCounter).dataset = datasetName;
            models(modelCounter).recordingIndex = namenum;
            models(modelCounter).filename = filename;
            models(modelCounter).clusterIndex = cc;
            models(modelCounter).centerSpin = centerSpin;
            models(modelCounter).region = region;
            models(modelCounter).N = N;
            models(modelCounter).muTarget = reducedMean;
            models(modelCounter).qTarget = qTarget;
            models(modelCounter).qTargetMat = reducedCorr;
            models(modelCounter).covTargetMat = reducedCorr - reducedMean * reducedMean';
            models(modelCounter).h = h;
            models(modelCounter).J = J;
            models(modelCounter).fitInfo = fitInfo;
        end
    end
end

fprintf("\nFinished fitting %d models.\n", numel(models));

%% ---------------- ESTIMATE JACOBIAN LOG DET AND SPECIFIC HEAT VERSUS T ----------------

for a = 1:numel(models)

    N = models(a).N;
    h = models(a).h;
    J = models(a).J;

    [pairI, pairJ] = pairList(N);
    p = N + numel(pairI);

    logdetVals = nan(size(Tvals));
    logdetPerParamVals = nan(size(Tvals));
    specificHeatVals = nan(size(Tvals));
    specificHeatPerSpinVals = nan(size(Tvals));

    fprintf("\nEstimating log det Jacobian and specific heat for model %d, N = %d, p = %d\n", ...
        a, N, p);

    for t = 1:numel(Tvals)

        T = Tvals(t);

        fprintf("  T = %.3f\n", T);

        samples = sampleIsingGlauber(h, J, T, mcJac);

        [specificHeatVals(t), specificHeatPerSpinVals(t)] = estimateSpecificHeat(samples, h, J, T);

        [logdetJ, logdetPerParam] = estimateJacobianLogdet(samples, pairI, pairJ, logdetJitterFrac);

        if jacobianParameterization == "baselineParams"
            logdetJ = logdetJ - p * log(T);
            logdetPerParam = logdetJ / p;
        end

        logdetVals(t) = logdetJ;
        logdetPerParamVals(t) = logdetPerParam;
    end

    models(a).Tvals = Tvals;
    models(a).specificHeatT = specificHeatVals;
    models(a).specificHeatPerSpinT = specificHeatPerSpinVals;
    models(a).logdetJ = logdetVals;
    models(a).logdetJPerParam = logdetPerParamVals;
    models(a).jacobianParameterization = jacobianParameterization;
end

%% ---------------- PLOTS ----------------

figure("Color", "w", "Position", [100, 100, 1550, 450]);
tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");

%% Panel 1: specific heat versus fictitious temperature

nexttile;
hold on;

for a = 1:numel(models)
    plot(models(a).Tvals, models(a).specificHeatPerSpinT, ...
        "LineWidth", 1.7, ...
        "DisplayName", modelLabel(models(a)));
end

xline(1, "k--", "T = 1", "HandleVisibility", "off");

xlabel("fictitious temperature T");
ylabel("C(T)/N");
title("Specific heat per spin");
legend("Location", "best");
grid on;
box on;

%% Panel 2: log determinant of Jacobian

nexttile;
hold on;

for a = 1:numel(models)
    plot(models(a).Tvals, models(a).logdetJ, ...
        "LineWidth", 1.7, ...
        "DisplayName", modelLabel(models(a)));
end

xline(1, "k--", "T = 1", "HandleVisibility", "off");

xlabel("fictitious temperature T");
ylabel("log det Jacobian");
title("Jacobian determinant");
legend("Location", "best");
grid on;
box on;

%% Panel 3: normalized log determinant

nexttile;
hold on;

for a = 1:numel(models)
    plot(models(a).Tvals, models(a).logdetJPerParam, ...
        "LineWidth", 1.7, ...
        "DisplayName", modelLabel(models(a)));
end

xline(1, "k--", "T = 1", "HandleVisibility", "off");

xlabel("fictitious temperature T");
ylabel("(1/p) log det Jacobian");
title("Size-normalized log determinant");
legend("Location", "best");
grid on;
box on;

sgtitle(sprintf("%s random-bond Ising, Jacobian convention: %s", ...
    datasetName, jacobianParameterization));

%% ---------------- FIT DIAGNOSTICS ----------------

figure("Color", "w", "Position", [150, 150, 1200, 430]);
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");

nexttile;
hold on;

for a = 1:numel(models)
    plot(models(a).fitInfo.maxAbsResidualHistory, ...
        "LineWidth", 1.5, ...
        "DisplayName", modelLabel(models(a)));
end

xlabel("iteration");
ylabel("max abs moment residual");
title("Moment-matching convergence");
legend("Location", "best");
grid on;
box on;

nexttile;
hold on;

for a = 1:numel(models)
    plot(models(a).fitInfo.rmsResidualHistory, ...
        "LineWidth", 1.5, ...
        "DisplayName", modelLabel(models(a)));
end

xlabel("iteration");
ylabel("RMS moment residual");
title("RMS residual");
legend("Location", "best");
grid on;
box on;

%% ---------------- SAVE RESULTS ----------------

saveFile = sprintf("RB_Ising_Allen_NN_gradient_jacobian_specificheat_%s.mat", ...
    matlab.lang.makeValidName(datasetName));

save(saveFile, "models", "Tvals", "fitOpts", "mcFit", "mcJac", ...
    "logdetJitterFrac", "jacobianParameterization", "-v7.3");

fprintf("\nSaved results to %s\n", saveFile);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function rec = loadAllenRecordingMoments(fullFileName)

    load(fullFileName);

    xs = cell_data.anterior_posterior_ccf_coordinate;
    ys = cell_data.dorsal_ventral_ccf_coordinate;
    zs = cell_data.left_right_ccf_coordinate;

    spiking_patterns = 2 .* X - 1;

    num_bins = size(spiking_patterns, 2);
    num_neurons = size(spiking_patterns, 1);

    if exist("num_cells", "var") && num_neurons ~= num_cells
        error("num_neurons ~= num_cells in file %s.", fullFileName);
    end

    % This is the pairwise second moment <s_i s_j>, not covariance.
    datacorr = spiking_patterns * spiking_patterns.' / num_bins;

    % Mean spin <s_i>.
    datamean = mean(spiking_patterns, 2);

    Locations = [xs, ys, zs];

    valid = all(isfinite(Locations), 2) & isfinite(datamean);

    datamean = datamean(valid);
    datacorr = datacorr(valid, valid);
    Locations = Locations(valid, :);

    rec.datamean = datamean(:);
    rec.datacorr = datacorr;
    rec.Locations = Locations;
    rec.num_neurons = numel(datamean);
end

function centers = chooseClusterCenters(datamean, nCenters, mode)

    n = numel(datamean);

    switch mode

        case "random"

            centers = randperm(n, min(nCenters, n));

        case "highestVariance"

            % For +/-1 variables, variance is 1 - mean^2.
            variances = 1 - datamean.^2;
            [~, order] = sort(variances, "descend");
            centers = order(1:min(nCenters, n));

        case "first"

            centers = 1:min(nCenters, n);

        otherwise

            error('centerSelection must be "random", "highestVariance", or "first".');
    end
end

function region = nearestNeighborRegion(distanceSMetric, centerSpin, N)

    distanceMetric = distanceSMetric(centerSpin, :);
    [~, I] = sort(distanceMetric, "ascend");

    % Exactly N nearest cells, including the center cell.
    region = I(1:N);
    region = region(:);
end

function [h, J, info] = fitIsingByGradientMatching(targetStats, N, pairI, pairJ, fitOpts, mcFit)
    % Fits pairwise Ising by gradient descent on the negative log likelihood.
    %
    % Objective:
    %   L(theta) = log Z(theta) - theta' * targetStats
    %
    % Gradient:
    %   grad L = <T>_model - targetStats
    %
    % Update:
    %   theta <- theta - stepSize * grad

    nPairs = numel(pairI);
    p = N + nPairs;

    muTarget = targetStats(1:N);

    switch fitOpts.initMode

        case "independent"

            h = atanh(max(min(muTarget, 0.999999), -0.999999));
            J = zeros(N, N);

        case "zero"

            h = zeros(N, 1);
            J = zeros(N, N);

        case "naiveMeanField"

            qTarget = zeros(N, N);
            qTarget(sub2ind([N, N], pairI, pairJ)) = targetStats(N+1:end);
            qTarget = qTarget + qTarget';
            qTarget(1:N+1:end) = 1;

            covTarget = qTarget - muTarget * muTarget';
            covTarget = 0.5 * (covTarget + covTarget');

            invC = pinv(covTarget);

            J = -invC;
            J(1:N+1:end) = 0;

            J = max(min(J, fitOpts.maxAbsInitialJ), -fitOpts.maxAbsInitialJ);
            J = 0.5 * (J + J');

            h = atanh(max(min(muTarget, 0.999999), -0.999999)) - J * muTarget;

        otherwise

            error('fitOpts.initMode must be "independent", "zero", or "naiveMeanField".');
    end

    theta = HJToTheta(h, J, pairI, pairJ);

    maxAbsResidualHistory = nan(fitOpts.maxIter, 1);
    rmsResidualHistory = nan(fitOpts.maxIter, 1);
    gradNormHistory = nan(fitOpts.maxIter, 1);
    thetaNormHistory = nan(fitOpts.maxIter, 1);

    for iter = 1:fitOpts.maxIter

        [h, J] = thetaToHJ(theta, N, pairI, pairJ);

        samples = sampleIsingGlauber(h, J, 1, mcFit);

        modelStats = estimateStatsOnly(samples, pairI, pairJ);

        % Gradient of negative log likelihood.
        grad = modelStats - targetStats;

        maxAbsRes = max(abs(grad));
        rmsRes = sqrt(mean(grad.^2));

        maxAbsResidualHistory(iter) = maxAbsRes;
        rmsResidualHistory(iter) = rmsRes;
        gradNormHistory(iter) = norm(grad);
        thetaNormHistory(iter) = norm(theta);

        if mod(iter, fitOpts.displayEvery) == 0
            fprintf("iter %4d, max abs residual %.4g, RMS residual %.4g, grad norm %.4g\n", ...
                iter, maxAbsRes, rmsRes, norm(grad));
        end

        if maxAbsRes < fitOpts.tol
            fprintf("Converged at iteration %d.\n", iter);
            break;
        end

        theta = theta - fitOpts.stepSize * grad;
    end

    [h, J] = thetaToHJ(theta, N, pairI, pairJ);

    info.theta = theta;
    info.maxAbsResidualHistory = maxAbsResidualHistory;
    info.rmsResidualHistory = rmsResidualHistory;
    info.gradNormHistory = gradNormHistory;
    info.thetaNormHistory = thetaNormHistory;
    info.nParams = p;
    info.N = N;
    info.stepSize = fitOpts.stepSize;
    info.initMode = fitOpts.initMode;
end

function statsVec = estimateStatsOnly(samples, pairI, pairJ)

    F = sufficientStatsPM(samples, pairI, pairJ);
    statsVec = mean(F, 1)';
end

function samples = sampleIsingGlauber(h, J, T, mc)

    N = numel(h);

    hT = h / T;
    JT = J / T;

    samplesPerChain = floor((mc.nSweeps - mc.burnSweeps) / mc.thinSweeps);
    totalSamples = mc.nChains * samplesPerChain;

    samples = zeros(totalSamples, N);
    sampleCounter = 0;

    for c = 1:mc.nChains

        switch mc.initMode

            case "random"
                s = 2 * randi([0, 1], 1, N) - 1;

            case "plus"
                s = ones(1, N);

            case "minus"
                s = -ones(1, N);

            otherwise
                error('mc.initMode must be "random", "plus", or "minus".');
        end

        for sweep = 1:mc.nSweeps

            order = randperm(N);

            for kk = 1:N

                i = order(kk);

                % Since J has zero diagonal, s * J(:,i) excludes self-term.
                localField = hT(i) + s * JT(:, i);

                pPlus = 1 / (1 + exp(-2 * localField));

                if rand < pPlus
                    s(i) = 1;
                else
                    s(i) = -1;
                end
            end

            if sweep > mc.burnSweeps && mod(sweep - mc.burnSweeps, mc.thinSweeps) == 0
                sampleCounter = sampleCounter + 1;
                samples(sampleCounter, :) = s;
            end
        end
    end

    samples = samples(1:sampleCounter, :);
end

function [C, CperSpin] = estimateSpecificHeat(samples, h, J, T)
    % Specific heat for the fitted +/-1 pairwise Ising model.
    %
    % Energy convention:
    %
    %   E(s) = -sum_i h_i s_i - sum_{i<j} J_ij s_i s_j
    %
    % The fictitious-temperature ensemble is:
    %
    %   P_T(s) proportional to exp[-E(s)/T]
    %
    % Therefore:
    %
    %   C(T) = Var_T(E) / T^2
    %
    % and the specific heat per spin is C(T)/N.

    N = size(samples, 2);

    fieldEnergy = -samples * h(:);

    % Since J is symmetric with zero diagonal, samples * J gives all
    % pair contributions twice when multiplied elementwise by samples.
    pairEnergy = -0.5 * sum((samples * J) .* samples, 2);

    E = fieldEnergy + pairEnergy;

    C = var(E, 1) / (T^2);
    CperSpin = C / N;
end

function [logdetJ, logdetPerParam] = estimateJacobianLogdet(samples, pairI, pairJ, jitterFrac)

    F = sufficientStatsPM(samples, pairI, pairJ);

    C = cov(F, 1);
    C = 0.5 * (C + C');

    p = size(C, 1);

    scale = mean(diag(C));
    jitter = max(jitterFrac * scale, 1e-12);

    Creg = C + jitter * eye(p);

    [R, flag] = chol(Creg);

    if flag == 0
        logdetJ = 2 * sum(log(diag(R)));
    else
        eigvals = eig(C);
        eigvals = max(real(eigvals), 0);
        logdetJ = sum(log(eigvals + jitter));
    end

    logdetPerParam = logdetJ / p;
end

function F = sufficientStatsPM(samples, pairI, pairJ)

    nSamples = size(samples, 1);
    N = size(samples, 2);
    nPairs = numel(pairI);

    F = zeros(nSamples, N + nPairs);

    F(:, 1:N) = samples;

    for e = 1:nPairs
        F(:, N + e) = samples(:, pairI(e)) .* samples(:, pairJ(e));
    end
end

function [pairI, pairJ] = pairList(N)

    [I, J] = find(triu(ones(N), 1));
    pairI = I(:);
    pairJ = J(:);
end

function theta = HJToTheta(h, J, pairI, pairJ)

    N = numel(h);
    nPairs = numel(pairI);

    theta = zeros(N + nPairs, 1);
    theta(1:N) = h(:);

    for e = 1:nPairs
        theta(N + e) = J(pairI(e), pairJ(e));
    end
end

function [h, J] = thetaToHJ(theta, N, pairI, pairJ)

    h = theta(1:N);

    Ju = theta(N + 1:end);

    J = zeros(N, N);

    idx1 = sub2ind([N, N], pairI, pairJ);
    idx2 = sub2ind([N, N], pairJ, pairI);

    J(idx1) = Ju;
    J(idx2) = Ju;
end

function label = modelLabel(model)

    label = sprintf("rec %d, c %d, N=%d", ...
        model.recordingIndex, model.clusterIndex, model.N);
end
