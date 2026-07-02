function [h, lambda] = hlambda(muTarget, chiTarget, N)
%HLAMBDA Invert (mu, chi) -> (h, lambda) using safeguarded continuation.
%
% Assumes MUCHI is available and returns:
%   [muVal, chiVal, dmu_dlambda] = muchi(h, lambda, N)
%
% Inputs:
%   muTarget   target mean magnetization
%   chiTarget  target susceptibility
%   N          system size
%
% Outputs:
%   h, lambda  inferred parameters
%
% Method:
%   1) Start near lambda = 0, where h ~ atanh(mu)
%   2) March along the constant-mu manifold using
%         dh/dlambda = -(dmu/dlambda)/(dmu/dh)
%   3) After each predictor step, do damped Newton corrections in h
%      to keep mu(h,lambda) close to muTarget
%   4) Adapt step size based on whether chi moves toward target
%
% Safeguards prioritize stable continuation over aggressive step growth.

    %% ---------------- user-tunable parameters ----------------
    maxIters      = 50000;   % max continuation steps
    tolChi        = 1e-8;    % tolerance on chi
    tolMu         = min(max(1-muTarget.^2, 1e-12),10e-7) ;% tolerance on mu
    dlambda0      = 5e-4;    % initial lambda step size
    dlambdaMin    = 1e-8;    % smallest allowed step size
    dlambdaMax    = 2e-4;    % largest allowed step size
    growthFactor  = 1.0;    % step growth after successful steps
    shrinkFactor  = 0.5;     % step shrink on bad steps
    dhMaxPred     = 0.01;    % max predictor correction in h
    dhMaxCorr     = 0.01;    % max Newton corrector step in h
    maxCorrIters  = 20;       % number of h-corrector iterations per step
    tiny          = 1e-14;   % small-number safeguard
    epsMu         = 1e-12;   % keeps atanh away from infinities
    %% ---------------------------------------------------------

    % Keep mu away from exactly +-1
    muTarget = min(max(muTarget, -1 + epsMu), 1 - epsMu);

    % Initial point near lambda = 0
    h = atanh(muTarget);
    lambda = 1e-4;

    % Choose initial side of lambda=0 using the lambda=0 susceptibility
    % At lambda = 0, chi = 1 - mu^2 under this normalization.
    if chiTarget < 1 - muTarget^2
        lambda = -1e-4;
    end

    % Initial lambda step
    dlambda = sign(lambda) * dlambda0;
    if dlambda == 0
        dlambda = dlambda0;
    end

    % Evaluate initial point
    [muNow, chiNow, dmu_dlambda] = muchi(h, lambda, N);
    
    if any(~isfinite([muNow, chiNow, dmu_dlambda]))
        error('Initial MUCHI evaluation returned non-finite values.')
    end

    errNow = chiNow - chiTarget;

    % Best point seen so far
    bestErr = abs(errNow);
    bestMuErr = abs(muNow - muTarget);
    bestH = h;
    bestLambda = lambda;

    % Immediate convergence check
    if abs(errNow) < tolChi && abs(muNow - muTarget) < tolMu
        return
    end

    for iter = 1:maxIters
        
        % dmu/dh equals chi in this parameterization.
        dmu_dh = chiNow;

        % If tangent is ill-defined, shrink step and retry
        if ~isfinite(dmu_dh) || abs(dmu_dh) < tiny
            dlambda = -sign(dlambda) * max(abs(dlambda) * shrinkFactor, dlambdaMin);
            if abs(dlambda) <= dlambdaMin
                break
            end
            continue
        end

        % ---------- predictor step along constant-mu manifold ----------
        dhPred = -(dmu_dlambda / dmu_dh) * dlambda;

        if ~isfinite(dhPred)
            dlambda = -sign(dlambda) * max(abs(dlambda) * shrinkFactor, dlambdaMin);
            if abs(dlambda) <= dlambdaMin
                break
            end
            continue
        end

        % Cap predictor step in h
        dhPred = sign(dhPred) * min(abs(dhPred), dhMaxPred);

        hTrial = h + dhPred;
        lambdaTrial = lambda + dlambda;

        if ~isfinite(hTrial) || ~isfinite(lambdaTrial)
            dlambda = -sign(dlambda) * max(abs(dlambda) * shrinkFactor, dlambdaMin);
            if abs(dlambda) <= dlambdaMin
                break
            end
            continue
        end

        % ---------- corrector: project back toward muTarget ----------
        goodTrial = true;
        for k = 1:maxCorrIters
            [muTrial, chiTrial, dmu_dlambdaTrial] = muchi(hTrial, lambdaTrial, N);

            if any(~isfinite([muTrial, chiTrial, dmu_dlambdaTrial]))
                goodTrial = false;
                break
            end

            muErr = muTrial - muTarget;

            if abs(muErr) < tolMu
                break
            end

            % Newton corrector in h at fixed lambda
            if abs(chiTrial) < tiny
                goodTrial = false;
                break
            end

            dhCorr = -muErr / chiTrial;

            if ~isfinite(dhCorr)
                goodTrial = false;
                break
            end

            % Cap corrector step
            dhCorr = sign(dhCorr) * min(abs(dhCorr), dhMaxCorr);
            hTrial = hTrial + dhCorr;

            if ~isfinite(hTrial)
                goodTrial = false;
                break
            end
        end

        if ~goodTrial
            dlambda = -sign(dlambda) * max(abs(dlambda) * shrinkFactor, dlambdaMin);
            if abs(dlambda) <= dlambdaMin
                break
            end
            continue
        end

        % Final evaluation of corrected trial point
        [muTrial, chiTrial, dmu_dlambdaTrial] = muchi(hTrial, lambdaTrial, N);

        if any(~isfinite([muTrial, chiTrial, dmu_dlambdaTrial]))
            dlambda = -sign(dlambda) * max(abs(dlambda) * shrinkFactor, dlambdaMin);
            if abs(dlambda) <= dlambdaMin
                break
            end
            continue
        end

        errTrial = chiTrial - chiTarget;
        muErrTrial = muTrial - muTarget;

        % Update best point seen
        if abs(errTrial) < bestErr || ...
           (abs(errTrial) == bestErr && abs(muErrTrial) < bestMuErr)
            bestErr = abs(errTrial);
            bestMuErr = abs(muErrTrial);
            bestH = hTrial;
            bestLambda = lambdaTrial;
        end

        % Convergence
        if abs(errTrial) < tolChi && abs(muErrTrial) < tolMu
            h = hTrial;
            lambda = lambdaTrial;
            return
        end

        % Decide whether step was helpful
        improvedChi = abs(errTrial) < abs(errNow);

        if improvedChi
            % Accept step
            h = hTrial;
            lambda = lambdaTrial;
            muNow = muTrial;
            chiNow = chiTrial;
            dmu_dlambda = dmu_dlambdaTrial;
            errNow = errTrial;

            % Grow step cautiously
            dlambda = sign(dlambda) * min(abs(dlambda) * growthFactor, dlambdaMax);
        else
            % If chi got worse, reverse and shrink
            dlambda = -sign(dlambda) * max(abs(dlambda) * shrinkFactor, dlambdaMin);

            if abs(dlambda) <= dlambdaMin
                break
            end
        end
    end

    % Fall back to best point found
    h = bestH;
    lambda = bestLambda;
end

% function [h, lambda] = hlambda(mu, chi, N)
% 
%     % Initial point near lambda = 0
%     h0 = atanh(mu);
%     lambda0 = 1e-4;
% 
%     % Decide which side of lambda=0 to start on
%     if chi < 1 - mu^2
%         lambda0 = -1e-4;
%     end
% 
%     h = h0;
%     lambda = lambda0;
% 
%     % Initial step size in lambda
%     dlambda = 5e-4;
%     dlambda_min = 1e-5;
%     dlambda_max = 1e-2;
% 
%     maxIters = 10000;
%     tolChi = 1e-2;
% 
%     % Track best point found
%     bestErr = inf;
%     bestH = h;
%     bestLambda = lambda;
% 
%     % Evaluate initial point
%     [~, xt, dudl] = muchi(h, lambda, N);
%     errPrev = xt - chi;
% 
%     if abs(errPrev) < bestErr
%         bestErr = abs(errPrev);
%         bestH = h;
%         bestLambda = lambda;
%     end
% 
%     for i = 1:maxIters
% 
%         % Tangent step along constant-mu manifold:
%         % dh/dlambda = -(dmu/dlambda)/(dmu/dh)
%         dudh = xt;
% 
%         % Safety check
%         if abs(dudh) < 1e-14
%             warning('dmu/dh became too small; stopping.');
%             break
%         end
% 
%         % Predictor step
%         hTrial = h - (dudl / dudh) * dlambda;
%         lambdaTrial = lambda + dlambda;
% 
%         % Evaluate trial point
%         [~, xtTrial, dudlTrial] = muchi(hTrial, lambdaTrial, N);
%         errTrial = xtTrial - chi;
% 
%         % Save best point seen so far
%         if abs(errTrial) < bestErr
%             bestErr = abs(errTrial);
%             bestH = hTrial;
%             bestLambda = lambdaTrial;
%         end
% 
%         % Converged
%         if abs(errTrial) < tolChi
%             disp(i)
%             h = real(hTrial);
%             lambda = real(lambdaTrial);
%             return
%         end
% 
%         % If the step crosses the target chi, reverse direction and shrink.
%         if sign(errTrial) ~= sign(errPrev)
%             dlambda = -dlambda / 2;
% 
%             % Stop if step is tiny
%             if abs(dlambda) < dlambda_min
%                 break
%             end
% 
%         else
%             % Accept the step
%             h = hTrial;
%             lambda = lambdaTrial;
%             xt = xtTrial;
%             dudl = dudlTrial;
%             errPrev = errTrial;
% 
%             % Increase the step size after accepted steps.
%             %dlambda = sign(dlambda) * min(1.2 * abs(dlambda), dlambda_max);
%         end
%     end
% 
%     % Return the best point found
%     h = real(bestH);
%     lambda = real(bestLambda);
% 
% end

% function [h,lambda] = hlambda(mu,chi, N)
% 
%     h0 = atanh(mu);
%     lambda0 = 0.0001;
% 
%     if chi < 1-mu^2
%         lambda0 = -0.0001;
%     end
% 
%     h = h0;
%     lambda = lambda0;
%     xs = [];
%     hs = [];
%     ls = [];
% 
%     iters = 5000;
% 
%     deltal = 0.0005; % Could come back to make this a dynamic search
%     for i = 1:iters
% 
%         [mut, xt, dudl] = muchi(h, lambda, N);
% 
%         dudh = xt;
% 
%         h = h - dudl/dudh * deltal;
% 
%         lambda = lambda + sign(lambda0)*deltal;
% 
%         xs = [xs,xt];
%         hs = [hs,h];
%         ls = [ls, lambda];
% 
%         if sign(lambda0) > 0
%             if 0 > 1.01*chi - xt
%                 break
%             end
%         else
%             if xt < 0.99*chi
%                 break
%             end
%         end
% 
%     end
% 
%     [M,I] = min(abs(chi-xs)) ;
% 
%     lambda = real(ls(I));
%     h = real(hs(I));
% 
% end
