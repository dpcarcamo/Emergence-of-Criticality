%% KEEP PLOTS

%% RG trajectories from fixed initial (mu0, chi0)
% Trajectory plotting script for the supplement.
% The dynamical system is defined in DhDNtot.m and DlamDNtot.m.

clearvars

%% Plot settings

N0 = 10;
Nf = 1000;
Nspan = [N0 Nf];

% Use every pair in the mu0 x chi0 grid below.
mu0Values = linspace(-0.8, -0.5, 4);
chi0Values = linspace(0.8, 3.5, 4);
initialConditions = makeInitialConditionGrid(mu0Values, chi0Values);

% Alternatively, specify exact (mu0, chi0) rows by hand:
% initialConditions = [
%     -0.30, 1.20
%     -0.50, 0.80
%     -0.60, 2.00
% ];

muColorLevels = unique(initialConditions(:, 1), 'stable');
muColors = lines(max(numel(muColorLevels), 1));

solver = @ode15s;
odeOptions = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-9);
axisLabelFontSize = 26;
tickLabelFontSize = 18;

% Invariant checks call muExact/chiExact repeatedly, so keep them off while
% exploring trajectory shapes.
checkInvariants = false;
invariantSampleCount = 25;

%% Plot initial conditions in mu0-chi0 plane

muSpan = max(initialConditions(:, 1)) - min(initialConditions(:, 1));
muPadding = max(0.05, 0.10*muSpan);
muMinPlot = max(-0.999, min(initialConditions(:, 1)) - muPadding);
muMaxPlot = min(0.999, max(initialConditions(:, 1)) + muPadding);

if muMinPlot >= muMaxPlot
    muCenter = min(max(mean(initialConditions(:, 1)), -0.9), 0.9);
    muMinPlot = max(-0.999, muCenter - 0.1);
    muMaxPlot = min(0.999, muCenter + 0.1);
end

muCurve = linspace(muMinPlot, muMaxPlot, 1000);
chiIndependent = 1 - muCurve.^2;
chiDenominator = muCurve - atanh(muCurve).*(1 - muCurve.^2);
chiCritical = muCurve.*(1 - muCurve.^2)./chiDenominator;

chiMinPlot = min(0, 0.95*min(initialConditions(:, 2)));
chiScale = max([initialConditions(:, 2); chiIndependent(:)]);
criticalForLimit = chiCritical( ...
    isfinite(chiCritical) & ...
    chiCritical > chiMinPlot & ...
    chiCritical < 5*chiScale);


chiMaxPlot = 1.15*chiScale;


if chiMaxPlot <= chiMinPlot
    chiMaxPlot = chiMinPlot + 1;
end

figure('Name', 'Initial conditions in mu0-chi0 plane');
hold on
grid on
box on

shadeBottom = max(chiCritical, chiMinPlot);
shadeMask = isfinite(shadeBottom) & shadeBottom < chiMaxPlot;
shadeIdx = find(shadeMask);

if ~isempty(shadeIdx)
    segmentBreaks = [0, find(diff(shadeIdx) > 1), numel(shadeIdx)];
    for segment = 1:numel(segmentBreaks)-1
        idx = shadeIdx(segmentBreaks(segment)+1:segmentBreaks(segment+1));
        xPatch = [muCurve(idx), fliplr(muCurve(idx))];
        yPatch = [shadeBottom(idx), chiMaxPlot*ones(size(idx))];
        fill(xPatch, yPatch, 'r', ...
            'FaceAlpha', 0.12, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

plot(muCurve, chiIndependent, ...
    'b-', ...
    'LineWidth', 1.25, ...
    'DisplayName', '1 - \mu^2');
plot(muCurve, chiCritical, ...
    'r-', ...
    'LineWidth', 1.25, ...
    'DisplayName', '\mu(1-\mu^2)/(\mu - atanh(\mu)(1-\mu^2))');
muLegendShown = false(size(muColorLevels));
for idx = 1:size(initialConditions, 1)
    mu0 = initialConditions(idx, 1);
    chi0 = initialConditions(idx, 2);
    muColorIdx = find(abs(muColorLevels - mu0) < ...
        10*eps(max(1, abs(mu0))), 1);

    if muLegendShown(muColorIdx)
        label = '';
        visibility = 'off';
    else
        label = sprintf('\\mu_0 = %.3g', mu0);
        visibility = 'on';
        muLegendShown(muColorIdx) = true;
    end

    plot(mu0, chi0, ...
        'o', ...
        'Color', muColors(muColorIdx, :), ...
        'MarkerFaceColor', muColors(muColorIdx, :), ...
        'DisplayName', label, ...
        'HandleVisibility', visibility);
end

xlim([muMinPlot, muMaxPlot])
ylim([chiMinPlot, chiMaxPlot])
set(gca, 'FontSize', tickLabelFontSize)
xlabel('\mu_0', 'FontSize', axisLabelFontSize)
ylabel('\chi_0', 'FontSize', axisLabelFontSize)
title('Initial conditions')

%legend('Location', 'best')
hold off
drawnow

axis square

%% Integrate trajectories

trajectories = integrateTrajectoryFamily( ...
    initialConditions, Nspan, solver, odeOptions, ...
    checkInvariants, invariantSampleCount);

%% Plot trajectories in h-lambda plane

figure('Name', 'RG trajectories in h-lambda plane');
hold on
grid on
box on


okIdx = find(strcmp({trajectories.status}, 'ok'));
legendShown = false(size(muColorLevels));

for k = 1:numel(okIdx)
    idx = okIdx(k);
    trajectory = trajectories(idx);
    muColorIdx = find(abs(muColorLevels - trajectory.mu0) < ...
        10*eps(max(1, abs(trajectory.mu0))), 1);
    color = muColors(muColorIdx, :);

    if legendShown(muColorIdx)
        label = '';
        visibility = 'off';
    else
        label = sprintf('\\mu_0 = %.3g', trajectory.mu0);
        visibility = 'on';
        legendShown(muColorIdx) = true;
    end

    plot(trajectory.h, trajectory.lambda, ...
        'LineWidth', 1.5, ...
        'Color', color, ...
        'DisplayName', label, ...
        'HandleVisibility', visibility);
    plot(trajectory.h(1), trajectory.lambda(1), ...
        'o', ...
        'Color', color, ...
        'MarkerFaceColor', color, ...
        'HandleVisibility', 'off');
    plot(trajectory.h(end), trajectory.lambda(end), ...
        'x', ...
        'Color', color, ...
        'LineWidth', 1.25, ...
        'HandleVisibility', 'off');
end

xlabel('h', 'FontSize', axisLabelFontSize)
ylabel('\lambda', 'FontSize', axisLabelFontSize)
set(gca, 'FontSize', tickLabelFontSize)
title(sprintf('RG trajectories, N = %.4g to %.4g', Nspan(1), Nspan(2)))

if ~isempty(okIdx) && numel(muColorLevels) <= 12
    legend('Location', 'best')
end
axis square
xscale log

hold off

if checkInvariants
    fprintf('\nInvariant drift summary:\n');
    for idx = 1:numel(trajectories)
        trajectory = trajectories(idx);
        if ~strcmp(trajectory.status, 'ok')
            fprintf('  mu0 = %.6g, chi0 = %.6g: failed (%s)\n', ...
                trajectory.mu0, trajectory.chi0, trajectory.message);
            continue
        end

        fprintf( ...
            '  mu0 = %.6g, chi0 = %.6g: max |dmu| = %.3e, max |dchi| = %.3e\n', ...
            trajectory.mu0, trajectory.chi0, ...
            trajectory.maxAbsMuDrift, trajectory.maxAbsChiDrift);
    end

    okInvariantIdx = find(strcmp({trajectories.status}, 'ok') & ...
        arrayfun(@(trajectory) ~isempty(trajectory.invariantN), trajectories));

    if ~isempty(okInvariantIdx)
        invariantColors = lines(numel(okInvariantIdx));
        figure('Name', 'Invariant drift');

        subplot(2, 1, 1)
        hold on
        grid on
        box on
        for k = 1:numel(okInvariantIdx)
            trajectory = trajectories(okInvariantIdx(k));
            plot(trajectory.invariantN, trajectory.muDrift, ...
                'Color', invariantColors(k, :), ...
                'LineWidth', 1.25);
        end
        xlabel('N', 'FontSize', axisLabelFontSize)
        ylabel('\mu(N) - \mu_0', 'FontSize', axisLabelFontSize)
        set(gca, 'FontSize', tickLabelFontSize)
        title('\mu drift')
        hold off

        subplot(2, 1, 2)
        hold on
        grid on
        box on
        for k = 1:numel(okInvariantIdx)
            trajectory = trajectories(okInvariantIdx(k));
            plot(trajectory.invariantN, trajectory.chiDrift, ...
                'Color', invariantColors(k, :), ...
                'LineWidth', 1.25);
        end
        xlabel('N', 'FontSize', axisLabelFontSize)
        ylabel('\chi(N) - \chi_0', 'FontSize', axisLabelFontSize)
        set(gca, 'FontSize', tickLabelFontSize)
        title('\chi drift')
        hold off
    end
end
axis square

%% Local functions

function initialConditions = makeInitialConditionGrid(mu0Values, chi0Values)
    [muGrid, chiGrid] = meshgrid(mu0Values(:), chi0Values(:));
    initialConditions = [muGrid(:), chiGrid(:)];
end

function trajectories = integrateTrajectoryFamily( ...
        initialConditions, Nspan, solver, odeOptions, ...
        checkInvariants, invariantSampleCount)

    nTrajectories = size(initialConditions, 1);
    trajectories = repmat(emptyTrajectory(), nTrajectories, 1);

    for idx = 1:nTrajectories
        mu0 = initialConditions(idx, 1);
        chi0 = initialConditions(idx, 2);

        fprintf( ...
            'Integrating %d/%d: mu0 = %.6g, chi0 = %.6g\n', ...
            idx, nTrajectories, mu0, chi0);

        trajectories(idx).mu0 = mu0;
        trajectories(idx).chi0 = chi0;

        try
            trajectories(idx) = integrateSingleTrajectory( ...
                mu0, chi0, Nspan, solver, odeOptions, ...
                checkInvariants, invariantSampleCount);
        catch err
            trajectories(idx).status = 'failed';
            trajectories(idx).message = err.message;
            warning( ...
                'RGderivatives:TrajectoryFailed', ...
                'Trajectory failed for mu0 = %.6g, chi0 = %.6g: %s', ...
                mu0, chi0, err.message);
        end
    end
end

function trajectory = integrateSingleTrajectory( ...
        mu0, chi0, Nspan, solver, odeOptions, ...
        checkInvariants, invariantSampleCount)

    trajectory = emptyTrajectory();
    trajectory.mu0 = mu0;
    trajectory.chi0 = chi0;

    [h0, lambda0] = hlambda(mu0, chi0, Nspan(1));
    initialState = real([h0; lambda0]);

    rhs = @(N, state) rgRhs(N, state);
    [N, state] = solver(rhs, Nspan, initialState, odeOptions);

    trajectory.N = N;
    trajectory.h = real(state(:, 1));
    trajectory.lambda = real(state(:, 2));
    trajectory.h0 = initialState(1);
    trajectory.lambda0 = initialState(2);
    trajectory.status = 'ok';

    if checkInvariants
        trajectory = addInvariantDiagnostics(trajectory, invariantSampleCount);
    end
end

function dState = rgRhs(N, state)
    h = state(1);
    lambda = state(2);

    dState = [
        DhDNtot(h, lambda, N)
        DlamDNtot(h, lambda, N)
    ];

    dState = real(dState);

    if any(~isfinite(dState))
        error('RGderivatives:NonFiniteRhs', ...
            'The RG right-hand side returned a non-finite value.');
    end
end

function trajectory = addInvariantDiagnostics(trajectory, invariantSampleCount)
    sampleCount = min(invariantSampleCount, numel(trajectory.N));
    sampleIdx = unique(round(linspace(1, numel(trajectory.N), sampleCount)));

    Nsample = trajectory.N(sampleIdx);
    hSample = trajectory.h(sampleIdx);
    lambdaSample = trajectory.lambda(sampleIdx);

    mu = arrayfun(@muExact, hSample, lambdaSample, Nsample);
    chi = arrayfun(@chiExact, hSample, lambdaSample, Nsample);

    trajectory.invariantN = Nsample;
    trajectory.mu = real(mu);
    trajectory.chi = real(chi);
    trajectory.muDrift = trajectory.mu - trajectory.mu0;
    trajectory.chiDrift = trajectory.chi - trajectory.chi0;
    trajectory.maxAbsMuDrift = max(abs(trajectory.muDrift));
    trajectory.maxAbsChiDrift = max(abs(trajectory.chiDrift));
end

function trajectory = emptyTrajectory()
    trajectory = struct( ...
        'mu0', NaN, ...
        'chi0', NaN, ...
        'N', [], ...
        'h', [], ...
        'lambda', [], ...
        'h0', NaN, ...
        'lambda0', NaN, ...
        'status', '', ...
        'message', '', ...
        'invariantN', [], ...
        'mu', [], ...
        'chi', [], ...
        'muDrift', [], ...
        'chiDrift', [], ...
        'maxAbsMuDrift', NaN, ...
        'maxAbsChiDrift', NaN);
end
