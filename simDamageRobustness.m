function results = simDamageRobustness(waveform, snrdB, N, bps, P, maxDelay, maxDoppler, guardWidth, samplingFreq, numTrials, damageType, damageParams, equalizerMode, numTaps, corruptThresholdFrac, progressFcn)
% results = simDamageRobustness(waveform, snrdB, N, bps, P, maxDelay, maxDoppler, guardWidth, ...
%           samplingFreq, numTrials, damageType, damageParams, equalizerMode, numTaps, ...
%           corruptThresholdFrac, progressFcn)
%
% Paired ("common random numbers") Monte Carlo comparison of BER/SER and
% per-symbol corruption spread with vs. without a deterministic
% channel-damage impairment (applyChannelDamage.m) layered on the usual
% P-path equal-strength Rayleigh doubly-dispersive channel (genDDChannel.m).
% Each trial reuses the SAME channel realization, transmit data, and noise
% for the damaged and control (undamaged) runs, isolating the damage's
% effect from Monte Carlo variance.
%
% Equalization ALWAYS uses genie CSI of the TRUE Rayleigh channel H only -
% the receiver does NOT know about the deterministic damage. This isolates
% each waveform's intrinsic structural resilience to the damage (the
% property under test), rather than conflating it with channel-estimation
% quality (already covered by simBERSER.m's genie/estimated CSI axis).
%
%   Inputs:
%      - waveform            : 'AFDM', 'OTFS' or 'OFDM'
%      - snrdB               : 1 x numSNR vector of SNR points [dB]
%      - N, bps, P, maxDelay, maxDoppler, guardWidth, samplingFreq, numTrials
%                             : as in simBERSER.m (channel is always 'Rayleigh')
%      - damageType          : 'freqNotch' | 'timeBurst' (see applyChannelDamage.m)
%      - damageParams        : struct passed to applyChannelDamage.m
%      - equalizerMode       : 'ideal' (default) or 'banded' (see equalizeChannel.m)
%      - numTaps             : banded mode only (default N)
%      - corruptThresholdFrac: fraction of mean |X|^2 defining "meaningfully
%                              corrupted" for the per-symbol metric (default 0.5)
%      - progressFcn         : optional progress callback
%
%   Output: results struct, each field 1 x numSNR unless noted:
%      - BER_control, BER_damaged, SER_control, SER_damaged
%      - newErrorRate               : fraction of symbol-slots correct in the
%                                      control run but wrong in the damaged run
%                                      (denominator N*numTrials) - symbols
%                                      newly broken by the damage
%      - corruptedSymbolRate        : fraction of symbol-slots whose deltaSqErr
%                                      exceeds corruptThresholdFrac*mean(|X|^2)
%      - meanDeltaSqErrPerCorrupted : average deltaSqErr among corrupted slots only
%      - totalDeltaSqErrEnergy      : mean(deltaSqErr) over all symbol-slots
%      - avgDeltaSqErrPerIndex      : numSNR x N matrix, average deltaSqErr at
%                                      each of the N symbol-index positions -
%                                      the direct "spread vs concentrated"
%                                      visualization data (compare against
%                                      check_damage_basis.m's structural PR
%                                      prediction)

if nargin < 13 || isempty(equalizerMode)
    equalizerMode = 'ideal';
end
if nargin < 14 || isempty(numTaps)
    numTaps = N;
end
if nargin < 15 || isempty(corruptThresholdFrac)
    corruptThresholdFrac = 0.5;
end
if nargin < 16 || isempty(progressFcn)
    progressFcn = @(x) [];
end

M = 2^bps;
numSNR = length(snrdB);

fields = {'BER_control','BER_damaged','SER_control','SER_damaged','newErrorRate', ...
    'corruptedSymbolRate','meanDeltaSqErrPerCorrupted','totalDeltaSqErrEnergy'};
for f = 1:length(fields)
    results.(fields{f}) = zeros(1, numSNR);
end
results.avgDeltaSqErrPerIndex = zeros(numSNR, N);

totalSteps = numSNR*numTrials;
step = 0;

for k = 1:numSNR
    snrLinear = 10^(snrdB(k)/10);
    bitErrsC = 0; bitErrsD = 0; symErrsC = 0; symErrsD = 0; newErrCount = 0;
    corruptCount = 0; sumDeltaSqErr = 0; sumDeltaSqErrCorrupt = 0;
    sumDeltaSqErrPerIndex = zeros(N,1);

    for t = 1:numTrials
        X_i = randi([0 M-1], N, 1);
        X = qammod(X_i, M, 'UnitAveragePower', true, 'InputType', 'integer');

        [S, cfg] = waveformMod(waveform, X, N, maxDoppler, guardWidth);

        if cfg.includeCCP
            H = genDDChannel('Rayleigh', N, samplingFreq, maxDelay, maxDoppler, P, true, true, cfg.c1);
        else
            H = genDDChannel('Rayleigh', N, samplingFreq, maxDelay, maxDoppler, P, false, true);
        end

        Sclean = H*S;
        noiseVar = mean(abs(Sclean).^2)/snrLinear;
        noise = sqrt(noiseVar/2)*(randn(N,1) + 1j*randn(N,1));   % SAME noise for both runs

        R_control = Sclean + noise;
        R_damaged = applyChannelDamage(R_control, samplingFreq, damageType, damageParams);

        Y_control = waveformDemod(waveform, R_control, N, cfg);
        Y_damaged = waveformDemod(waveform, R_damaged, N, cfg);

        G = waveformChannelMatrix(waveform, H, N, maxDoppler, guardWidth);   % genie H, oblivious to damage

        Xhat_c = equalizeChannel(G, Y_control, equalizerMode, numTaps, noiseVar);
        Xhat_d = equalizeChannel(G, Y_damaged, equalizerMode, numTaps, noiseVar);

        Xi_c = qamdemod(Xhat_c, M, 'UnitAveragePower', true, 'OutputType', 'integer');
        Xi_d = qamdemod(Xhat_d, M, 'UnitAveragePower', true, 'OutputType', 'integer');

        bitErrsC = bitErrsC + biterr(X_i, Xi_c, bps);
        bitErrsD = bitErrsD + biterr(X_i, Xi_d, bps);
        symErrsC = symErrsC + symerr(X_i, Xi_c);
        symErrsD = symErrsD + symerr(X_i, Xi_d);

        newErrorMask = (Xi_d ~= X_i) & (Xi_c == X_i);
        newErrCount = newErrCount + sum(newErrorMask);

        sqErrC = abs(Xhat_c - X).^2;
        sqErrD = abs(Xhat_d - X).^2;
        deltaSqErr = sqErrD - sqErrC;
        corruptedMask = deltaSqErr > corruptThresholdFrac*mean(abs(X).^2);

        corruptCount = corruptCount + sum(corruptedMask);
        sumDeltaSqErr = sumDeltaSqErr + sum(deltaSqErr);
        sumDeltaSqErrCorrupt = sumDeltaSqErrCorrupt + sum(deltaSqErr(corruptedMask));
        sumDeltaSqErrPerIndex = sumDeltaSqErrPerIndex + deltaSqErr;

        step = step + 1;
        progressFcn(step/totalSteps);
    end

    results.BER_control(k) = bitErrsC/(N*bps*numTrials);
    results.BER_damaged(k) = bitErrsD/(N*bps*numTrials);
    results.SER_control(k) = symErrsC/(N*numTrials);
    results.SER_damaged(k) = symErrsD/(N*numTrials);
    results.newErrorRate(k) = newErrCount/(N*numTrials);
    results.corruptedSymbolRate(k) = corruptCount/(N*numTrials);
    results.meanDeltaSqErrPerCorrupted(k) = sumDeltaSqErrCorrupt/max(corruptCount,1);
    results.totalDeltaSqErrEnergy(k) = sumDeltaSqErr/(N*numTrials);
    results.avgDeltaSqErrPerIndex(k,:) = sumDeltaSqErrPerIndex'/numTrials;
end

end
