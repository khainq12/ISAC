function [BER, SER] = simBERSER_DopplerMode(waveform, channelType, snrdB, N, bps, P, maxDelay, maxDoppler, guardWidth, samplingFreq, numTrials, equalizerMode, numTaps, csiMode, dopplerMode, progressFcn)
% [BER, SER] = simBERSER_DopplerMode(waveform, channelType, snrdB, N, bps, P, maxDelay, maxDoppler, ...
%           guardWidth, samplingFreq, numTrials, equalizerMode, numTaps, csiMode, dopplerMode, progressFcn)
%
% Same Monte Carlo BER/SER simulation as simBERSER.m (see that file for the
% established per-trial pattern this mirrors exactly), with one addition:
% dopplerMode controls whether each path's Doppler shift is left as the
% generically continuous (almost surely non-integer) value genDDChannel.m
% produces ('fractional', matches simBERSER.m's normal behavior and every
% other experiment in this project), or is rounded to the nearest INTEGER
% normalized Doppler bin before the channel matrix is built ('integer').
%
% This isolates the theoretical distinction from the theory doc
% ("LyThuyet_AFDM_OTFS_OFDM.docx", Muc 2.2, Eq 2.6 vs Eq 2.7): integer
% Doppler gives a pure grid shift with zero leakage (no ICI/IDI at all),
% fractional Doppler gives localized leakage - a property specific to how
% each waveform's own domain represents the channel, distinct from the
% multipath diversity-order argument (Muc 3-5).
%
%   Inputs: same as simBERSER.m, plus
%      - dopplerMode : 'fractional' (default, same as simBERSER.m) or
%                      'integer' (rounds every path's normalized Doppler to
%                      the nearest integer bin before building H, reusing
%                      cconvLTVChannel.m directly - no new channel-generation
%                      logic, only a post-processing rounding step on
%                      genDDChannel.m's own returned path parameters)
%   Outputs:
%      - BER, SER : 1 x numSNR bit/symbol error rate curves

if nargin < 14 || isempty(csiMode)
    csiMode = 'genie';
end
if nargin < 15 || isempty(dopplerMode)
    dopplerMode = 'fractional';
end
if nargin < 16 || isempty(progressFcn)
    progressFcn = @(x) [];
end

M = 2^bps;
numSNR = length(snrdB);
BER = zeros(1, numSNR);
SER = zeros(1, numSNR);

totalSteps = numSNR*numTrials;
step = 0;

for k = 1:numSNR
    bitErrs = 0;
    symErrs = 0;
    snrLinear = 10^(snrdB(k)/10);
    for t = 1:numTrials
        X_i = randi([0 M-1], N, 1);
        X = qammod(X_i, M, 'UnitAveragePower', true, 'InputType', 'integer');

        [S, cfg] = waveformMod(waveform, X, N, maxDoppler, guardWidth);

        if cfg.includeCCP
            [H, pathDelays, pathDopplers, fadingCoefs] = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, true, true, cfg.c1);
        else
            [H, pathDelays, pathDopplers, fadingCoefs] = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, false, true);
        end

        if strcmpi(dopplerMode, 'integer')
            normDopplersInt = round(N*pathDopplers/samplingFreq);
            pathDopplersInt = normDopplersInt*samplingFreq/N;
            if cfg.includeCCP
                H = cconvLTVChannel(N, pathDelays, pathDopplersInt, fadingCoefs, samplingFreq, true, cfg.c1);
            else
                H = cconvLTVChannel(N, pathDelays, pathDopplersInt, fadingCoefs, samplingFreq, false);
            end
        end

        Sclean = H*S;
        noiseVar = mean(abs(Sclean).^2)/snrLinear;
        R = Sclean + sqrt(noiseVar/2)*(randn(N,1) + 1j*randn(N,1));

        Y = waveformDemod(waveform, R, N, cfg);

        if strcmpi(csiMode, 'estimated')
            Xpilot = ones(N, 1);
            [Spilot, ~] = waveformMod(waveform, Xpilot, N, maxDoppler, guardWidth);
            Rpilot = H*Spilot + sqrt(noiseVar/2)*(randn(N,1) + 1j*randn(N,1));
            if cfg.includeCCP
                Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, true, cfg.c1);
            else
                Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, false, 0);
            end
            G = waveformChannelMatrix(waveform, Hhat, N, maxDoppler, guardWidth);
        else
            G = waveformChannelMatrix(waveform, H, N, maxDoppler, guardWidth);
        end
        X_hat = equalizeChannel(G, Y, equalizerMode, numTaps, noiseVar);

        X_hat_i = qamdemod(X_hat, M, 'UnitAveragePower', true, 'OutputType', 'integer');

        bitErrs = bitErrs + biterr(X_i, X_hat_i, bps);
        symErrs = symErrs + symerr(X_i, X_hat_i);

        step = step + 1;
        progressFcn(step/totalSteps);
    end
    BER(k) = bitErrs/(N*bps*numTrials);
    SER(k) = symErrs/(N*numTrials);
end

end
