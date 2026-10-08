function [rangeRMSE, velRMSE, crbRangeStd, crbVelStd] = simSensing(waveform, channelType, snrdB, N, bps, P, maxDelay, maxDoppler, guardWidth, samplingFreq, carrierFreq, numTrials, progressFcn, dopplerRes)
% [rangeRMSE, velRMSE, crbRangeStd, crbVelStd] = simSensing(waveform, channelType, snrdB, N, bps, P, maxDelay, maxDoppler, guardWidth, samplingFreq, carrierFreq, numTrials, progressFcn, dopplerRes)
% Monte Carlo range/velocity estimation error of the dominant (target, path 1) scatterer,
% extracted from the effective channel matrix via a matched-filter delay-Doppler search,
% together with the theoretical CRB-derived standard deviation for comparison.
%   Inputs: see simBERSER.m for shared parameters, plus
%      - carrierFreq : carrier frequency [Hz], used to convert Doppler to velocity
%      - dopplerRes  : Doppler search grid resolution passed to estimateDelayDoppler
%                       (default 0.02 - see estimateDelayDoppler.m). The estimator's
%                       physical Doppler/velocity resolution is dopplerRes*(samplingFreq/N),
%                       so for a sensing-only configuration a much finer value is normally
%                       needed to avoid the RMSE floor being dominated by grid quantization
%                       rather than noise.
%   Outputs:
%      - rangeRMSE, velRMSE       : 1 x numSNR empirical RMSE curves [m], [m/s]
%      - crbRangeStd, crbVelStd   : 1 x numSNR CRB-derived standard deviation curves [m], [m/s]

if nargin < 13 || isempty(progressFcn)
    progressFcn = @(x) [];
end
if nargin < 14 || isempty(dopplerRes)
    dopplerRes = 0.02;
end

lightspeed = physconst("lightspeed");
carrierWavelength = lightspeed/carrierFreq;

M = 2^bps;
numSNR = length(snrdB);
rangeRMSE = zeros(1, numSNR);
velRMSE = zeros(1, numSNR);
crbRangeStd = zeros(1, numSNR);
crbVelStd = zeros(1, numSNR);

totalSteps = numSNR*numTrials;
step = 0;
repS = [];

for k = 1:numSNR
    sqErrRange = 0;
    sqErrVel = 0;
    for t = 1:numTrials
        X_i = randi([0 M-1], N, 1);
        X = qammod(X_i, M, 'UnitAveragePower', true, 'InputType', 'integer');
        [S, cfg] = waveformMod(waveform, X, N, maxDoppler, guardWidth);

        if isempty(repS)
            repS = S;
        end

        if cfg.includeCCP
            [H, pathDelays, pathDopplers] = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, true, false, cfg.c1);
        else
            [H, pathDelays, pathDopplers] = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, false, false);
        end

        R = awgn(H*S, snrdB(k), 'measured');

        if cfg.includeCCP
            [delayIdxHat, dopplerHat] = estimateDelayDoppler(R, S, N, maxDelay, maxDoppler, true, cfg.c1, dopplerRes);
        else
            [delayIdxHat, dopplerHat] = estimateDelayDoppler(R, S, N, maxDelay, maxDoppler, false, 0, dopplerRes);
        end

        delayHat = delayIdxHat/samplingFreq;
        dopplerHzHat = dopplerHat/N*samplingFreq;

        rangeHat = lightspeed*delayHat/2;
        velHat = carrierWavelength*dopplerHzHat/2;

        rangeTrue = lightspeed*pathDelays(1)/2;
        velTrue = carrierWavelength*pathDopplers(1)/2;

        sqErrRange = sqErrRange + (rangeHat - rangeTrue)^2;
        sqErrVel = sqErrVel + (velHat - velTrue)^2;

        step = step + 1;
        progressFcn(step/totalSteps);
    end
    rangeRMSE(k) = sqrt(sqErrRange/numTrials);
    velRMSE(k) = sqrt(sqErrVel/numTrials);

    snrLinear = 10^(snrdB(k)/10);
    [crbDelay, crbDoppler] = crbDelayDoppler(repS, samplingFreq, snrLinear);
    crbRangeStd(k) = lightspeed/2*sqrt(crbDelay);
    crbVelStd(k) = carrierWavelength/2*sqrt(crbDoppler);
end

end
