% Parse salamander retinal spike patterns and fit an approximate
% random-bond Ising model with Monte Carlo moment matching.

clear
clc

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(repoRoot);

Data =  matfile(fullfile(scriptDir, "salamander_processed_20ms.mat"));

%% Load spike patterns and compute empirical moments

spiking_patterns = Data.spike_patterns;
spiking_patterns = spiking_patterns(:,:)*2-1;
trial_labels = Data.trial_labels;
num_bins = Data.num_bins;
num_nuerons = size(spiking_patterns,1);

datacorr = full(spiking_patterns*spiking_patterns.')/num_bins;
datamean = full(mean(spiking_patterns.'));
mu = mean(datamean)
chi = (sum(datacorr,'all') - sum(datamean,"all").^2)/num_nuerons
(1-mu^2)*num_nuerons    

Ts = 1./linspace(0.2,2 ,2000);

[h1, l] = hlambda(mu,chi,num_nuerons);

%% Curie-Weiss thermodynamic curves matched to empirical mu/chi
CC = [];
X = [];
M = [];

for i = 1:length(Ts)
    T = Ts(i);
    [mt,xt,~,~,ct] = muChiExact2Spin(h1/T,l/T ,num_nuerons);
    M = [M,mt];
    X = [X,xt/T];
    CC = [CC, ct];

end
%

figure(1)
plot(Ts, CC/num_nuerons, 'r')
hold on
plot(Ts, X/num_nuerons, 'b')


figure(2)
hold on
plot(M, X, 'o')
%plot(x, x.*(1-x.^2)./(x-atanh(x).*(1-x.^2)), 'LineWidth', 2, 'Color', 'red')
ylim([0,50])

figure(2)
hold on
[mt,xt,~,~,ct, wk] = muChiExact2Spin(h1,l ,num_nuerons);
plot(mt, xt, 'square')


figure(3)
hold on
plot(h1./Ts, l./Ts, 'ob', 'MarkerSize', 5)
axis square
xlabel('h', 'FontSize', 18)
ylabel('\lambda ', 'FontSize', 18)
plot(h1,l, 'or')



%% Fit random-bond Ising model to empirical moments

[h, J, info] = fitRandomBondIsingMC(datamean, datacorr, ...
    'LearningRate', 0.01, ...
    'Momentum', 0.8, ...
    'MaxIter', 5000, ...
    'Tol', 1e-5, ...
    'NumChains', 50, ...
    'BurnInSweeps', 500, ...
    'SweepsPerSample', 10, ...
    'NumSamples', 50, ...
    'Verbose', true);



plot(datacorr(:), info.modelPairMeans(:), '.')
hold on
plot(datacorr(:), datacorr(:), '.')


%% Pairwise Ising thermodynamics
Tlist = linspace(0.5, 3, 80);

thermo = isingThermoVsT_MC(h, J, Tlist, ...
    'NumChains', 8, ...
    'BurnInSweeps', 500, ...
    'SweepsPerSample', 10, ...
    'NumSamples', 500, ...
    'Verbose', true);

figure;
plot(thermo.T, thermo.chiPerSpin, 'LineWidth', 2);
xlabel('T');
ylabel('\chi / N');
title('Susceptibility per spin');

figure;
plot(thermo.T, thermo.heatPerSpin, 'LineWidth', 2);
xlabel('T');
ylabel('C / N');
title('Specific heat per spin');

function [h, J, info] = fitRandomBondIsingMC(means, pairMeans, varargin)
% fitRandomBondIsingMC
%
% Fit a random-bond Ising model in +-1 notation for large systems using
% Monte Carlo moment matching with Glauber dynamics.
%
% Model:
%   P(s) propto exp( sum_i h_i s_i + 1/2 sum_{i,j} J_ij s_i s_j )
%
% Inputs
%   means      : Nx1 vector, means(i) = <s_i>
%   pairMeans  : NxN matrix, pairMeans(i,j) = <s_i s_j>
%
% Optional name-value pairs
%   'LearningRate'   : gradient step size for h and J (default 0.02)
%   'Momentum'       : inertia term in [0,1) (default 0.9)
%   'MaxIter'        : maximum fitting iterations (default 500)
%   'Tol'            : target max abs moment mismatch (default 1e-5)
%   'Verbose'        : true/false (default true)
%   'NumChains'      : number of MC chains for model estimates (default 8)
%   'BurnInSweeps'   : burn-in sweeps per fitting iteration (default 200)
%   'SweepsPerSample': MC sweeps between saved samples (default 10)
%   'NumSamples'     : saved samples per chain per fitting iteration (default 50)
%   'JInitScale'     : scale for initial coupling guess (default 0.05)
%   'ClipField'      : clip |h_i| to this value (default 4)
%   'ClipCoupling'   : clip |J_ij| to this value (default 1)
%   'UsePersistentChains' : reuse chains between iterations (default true)
%
% Outputs
%   h     : Nx1 fitted fields
%   J     : NxN symmetric coupling matrix, diagonal zero
%   info  : diagnostics struct
%
% Notes
%   - This is approximate, not exact.
%   - For large N, pairMeans is still NxN, so memory can become large.
%   - Noisy pairMeans may require a looser tolerance than 1e-5.

    p = inputParser;
    addParameter(p, 'LearningRate', 0.02);
    addParameter(p, 'Momentum', 0.9);
    addParameter(p, 'MaxIter', 500);
    addParameter(p, 'Tol', 1e-5);
    addParameter(p, 'Verbose', true);

    addParameter(p, 'NumChains', 8);
    addParameter(p, 'BurnInSweeps', 200);
    addParameter(p, 'SweepsPerSample', 10);
    addParameter(p, 'NumSamples', 50);

    addParameter(p, 'JInitScale', 0.05);
    addParameter(p, 'ClipField', 4);
    addParameter(p, 'ClipCoupling', 10);
    addParameter(p, 'UsePersistentChains', true);

    parse(p, varargin{:});

    lr = p.Results.LearningRate;
    mom = p.Results.Momentum;
    maxIter = p.Results.MaxIter;
    tol = p.Results.Tol;
    verbose = p.Results.Verbose;

    numChains = p.Results.NumChains;
    burnInSweeps = p.Results.BurnInSweeps;
    sweepsPerSample = p.Results.SweepsPerSample;
    numSamples = p.Results.NumSamples;

    JInitScale = p.Results.JInitScale;
    clipField = p.Results.ClipField;
    clipCoupling = p.Results.ClipCoupling;
    usePersistentChains = p.Results.UsePersistentChains;

    means = means(:);
    N = length(means);

    if ~isequal(size(pairMeans), [N N])
        error('pairMeans must be NxN with N = length(means).');
    end

    pairMeans = 0.5 * (pairMeans + pairMeans.');
    pairMeans(1:N+1:end) = 1;

    if any(abs(means) >= 1)
        warning('Some means are very close to ±1. Initial h values will be clipped.');
    end

    % Initial field guess from independent model
    mClip = min(max(means, -0.999), 0.999);
    h = atanh(mClip);

    % Initial J guess from connected correlations
    C = pairMeans - means * means.';
    J = zeros(N, N);
    for i = 1:N
        for j = i+1:N
            denom = sqrt(max(1 - means(i)^2, 1e-12) * max(1 - means(j)^2, 1e-12));
            rho = C(i,j) / denom;
            rho = min(max(rho, -0.95), 0.95);
            J(i,j) = JInitScale * atanh(rho);
            J(j,i) = J(i,j);
        end
    end
    J(1:N+1:end) = 0;

    h = max(min(h, clipField), -clipField);
    J = max(min(J, clipCoupling), -clipCoupling);
    J = 0.5 * (J + J.');
    J(1:N+1:end) = 0;


    % Momentum buffers
    vh = zeros(N,1);
    vJ = zeros(N,N);

    % Persistent chains
    if usePersistentChains
        states = sign(randn(N, numChains));
        states(states == 0) = 1;
    else
        states = [];
    end

    errHist = nan(maxIter,1);
    meanErrHist = nan(maxIter,1);
    pairErrHist = nan(maxIter,1);

    for iter = 1:maxIter
        if ~usePersistentChains || isempty(states)
            states = sign(randn(N, numChains));
            states(states == 0) = 1;
        end

        [modelMeans, modelPairMeans, states] = estimateModelMomentsMC( ...
            h, J, states, burnInSweeps, sweepsPerSample, numSamples);

        grad_h = means - modelMeans;
        grad_J = pairMeans - modelPairMeans;

        grad_J = 0.5 * (grad_J + grad_J.');
        grad_J(1:N+1:end) = 0;

        vh = mom * vh + lr * grad_h;
        vJ = mom * vJ + lr * grad_J;

        h = h + vh;
        J = J + vJ;

        % Symmetry and clipping
        h = max(min(h, clipField), -clipField);
        J = 0.5 * (J + J.');
        J(1:N+1:end) = 0;
        J = max(min(J, clipCoupling), -clipCoupling);
        J(1:N+1:end) = 0;

        meanErr = max(abs(grad_h));
        pairErr = max(abs(grad_J(:)));
        err = max(meanErr, pairErr);

        errHist(iter) = err;
        meanErrHist(iter) = meanErr;
        pairErrHist(iter) = pairErr;

        if verbose && (iter == 1 || mod(iter,10) == 0)
            fprintf('iter %4d | max mean err = %.3e | max pair err = %.3e | max err = %.3e\n | %.3e\n', ...
                iter, meanErr, pairErr, err, mean(abs(grad_J(:))) + mean(abs(grad_h(:))));
        end

        if err < tol
            if verbose
                fprintf('Converged at iter %d with max error %.3e\n', iter, err);
            end
            errHist = errHist(1:iter);
            meanErrHist = meanErrHist(1:iter);
            pairErrHist = pairErrHist(1:iter);
            break;
        end
    end

    if iter == maxIter
        errHist = errHist(1:iter);
        meanErrHist = meanErrHist(1:iter);
        pairErrHist = pairErrHist(1:iter);
        if verbose
            fprintf('Stopped at MaxIter. Final max error = %.3e\n', errHist(end));
        end
    end

    info = struct();
    info.iterations = iter;
    info.finalError = errHist(end);
    info.errorHistory = errHist;
    info.meanErrorHistory = meanErrHist;
    info.pairErrorHistory = pairErrHist;
    info.modelMeans = modelMeans;
    info.modelPairMeans = modelPairMeans;
end


function [modelMeans, modelPairMeans, states] = estimateModelMomentsMC( ...
    h, J, states, burnInSweeps, sweepsPerSample, numSamples)
% Estimate <s_i> and <s_i s_j> from Glauber MC.
%
% states is N x numChains, each column one chain.

    [N, numChains] = size(states);

    % Burn-in
    states = glauberSweepMulti(states, h, J, burnInSweeps);

    sumS = zeros(N,1);
    sumSS = zeros(N,N);

    totalSaved = 0;

    for t = 1:numSamples
        states = glauberSweepMulti(states, h, J, sweepsPerSample);

        sumS = sumS + sum(states, 2);
        for c = 1:numChains
            s = states(:,c);
            sumSS = sumSS + (s * s.');
        end
        totalSaved = totalSaved + numChains;
    end

    modelMeans = sumS / totalSaved;
    modelPairMeans = sumSS / totalSaved;
    modelPairMeans = 0.5 * (modelPairMeans + modelPairMeans.');
    modelPairMeans(1:N+1:end) = 1;
end


function states = glauberSweepMulti(states, h, J, nSweeps)
% Run nSweeps Glauber sweeps on multiple chains.
%
% states: N x numChains, entries ±1

    [N, numChains] = size(states);

    for sweep = 1:nSweeps
        for step = 1:N
            i = randi(N);

            localField = h(i) + J(i,:) * states;    % 1 x numChains
            pPlus = 1 ./ (1 + exp(-2 * localField));

            states(i,:) = 2 * (rand(1, numChains) < pPlus) - 1;
        end
    end
end


function thermo = isingThermoVsT_MC(h, J, Tlist, varargin)
% isingThermoVsT_MC
%
% Estimate susceptibility and specific heat vs temperature using Glauber MC
% for the Ising model in +-1 notation.
%
% Model:
%   P_T(s) propto exp( [sum_i h_i s_i + 1/2 sum_{i,j} J_ij s_i s_j] / T )
%
% Hamiltonian:
%   H(s) = -sum_i h_i s_i - 1/2 sum_{i,j} J_ij s_i s_j
%
% Inputs
%   h     : Nx1 fields
%   J     : NxN symmetric couplings, diagonal zero
%   Tlist : vector of temperatures
%
% Optional name-value pairs
%   'NumChains'       : number of chains per temperature (default 8)
%   'BurnInSweeps'    : burn-in sweeps (default 500)
%   'SweepsPerSample' : sweeps between saved samples (default 10)
%   'NumSamples'      : number of saved samples per chain (default 200)
%   'Verbose'         : true/false (default true)
%   'UseAnnealing'    : use previous T state as init for next T (default true)
%
% Output
%   thermo struct with fields:
%       .T
%       .meanM
%       .meanm
%       .energy
%       .energyPerSpin
%       .chiTotal
%       .chiPerSpin
%       .heatCapacity
%       .heatPerSpin

    p = inputParser;
    addParameter(p, 'NumChains', 8);
    addParameter(p, 'BurnInSweeps', 500);
    addParameter(p, 'SweepsPerSample', 10);
    addParameter(p, 'NumSamples', 200);
    addParameter(p, 'Verbose', true);
    addParameter(p, 'UseAnnealing', true);
    parse(p, varargin{:});

    numChains = p.Results.NumChains;
    burnInSweeps = p.Results.BurnInSweeps;
    sweepsPerSample = p.Results.SweepsPerSample;
    numSamples = p.Results.NumSamples;
    verbose = p.Results.Verbose;
    useAnnealing = p.Results.UseAnnealing;

    h = h(:);
    N = length(h);

    if ~isequal(size(J), [N N])
        error('J must be NxN with N = length(h).');
    end

    J = 0.5 * (J + J.');
    J(1:N+1:end) = 0;

    Tlist = Tlist(:);
    nT = length(Tlist);

    meanM = zeros(nT,1);
    meanm = zeros(nT,1);
    meanE = zeros(nT,1);
    meanEPerSpin = zeros(nT,1);
    chiTotal = zeros(nT,1);
    chiPerSpin = zeros(nT,1);
    heatCapacity = zeros(nT,1);
    heatPerSpin = zeros(nT,1);

    states = sign(randn(N, numChains));
    states(states == 0) = 1;

    for k = 1:nT
        T = Tlist(k);
        if T <= 0
            error('Temperatures must be positive.');
        end

        if verbose
            fprintf('Temperature %d / %d, T = %.5g\n', k, nT, T);
        end

        hT = h / T;
        JT = J / T;

        if ~useAnnealing || k == 1
            states = sign(randn(N, numChains));
            states(states == 0) = 1;
        end

        states = glauberSweepMulti(states, hT, JT, burnInSweeps);

        Mvals = zeros(numChains * numSamples, 1);
        Evals = zeros(numChains * numSamples, 1);

        idx = 1;
        for sIdx = 1:numSamples
            states = glauberSweepMulti(states, hT, JT, sweepsPerSample);

            for c = 1:numChains
                s = states(:,c);
                Mvals(idx) = sum(s);
                Evals(idx) = -h.' * s - 0.5 * s.' * J * s;
                idx = idx + 1;
            end
        end

        EM = mean(Mvals);
        EM2 = mean(Mvals.^2);
        EE = mean(Evals);
        EE2 = mean(Evals.^2);

        varM = EM2 - EM^2;
        varE = EE2 - EE^2;

        meanM(k) = EM;
        meanm(k) = EM / N;
        meanE(k) = EE;
        meanEPerSpin(k) = EE / N;

        chiTotal(k) = varM / T;
        chiPerSpin(k) = varM / (N * T);

        heatCapacity(k) = varE / (T^2);
        heatPerSpin(k) = varE / (N * T^2);
    end

    thermo = struct();
    thermo.T = Tlist;
    thermo.meanM = meanM;
    thermo.meanm = meanm;
    thermo.energy = meanE;
    thermo.energyPerSpin = meanEPerSpin;
    thermo.chiTotal = chiTotal;
    thermo.chiPerSpin = chiPerSpin;
    thermo.heatCapacity = heatCapacity;
    thermo.heatPerSpin = heatPerSpin;
end


