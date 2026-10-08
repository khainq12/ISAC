function [BER] = simDiversityML_FixedPaths(waveform, snrdB, N, maxDoppler, guardWidth, ell1, ell2, fd1, fd2, numTrials, progressFcn)
% BER = simDiversityML_FixedPaths(waveform, snrdB, N, maxDoppler, guardWidth, ell1, ell2, fd1, fd2, numTrials, progressFcn)
%
% Same exhaustive-ML BPSK diversity demo as simDiversityML.m (mirrors its
% per-trial pattern exactly - fading h1,h2 random each trial, genie CSI
% only), EXCEPT the two paths' delay/Doppler (ell1, ell2, fd1, fd2) are
% FIXED inputs instead of drawn randomly each trial.
%
% Purpose: directly test the exact mechanism from AFDM Theorem 1's proof
% (theory doc "LyThuyet_AFDM_OTFS_OFDM.docx", Muc 4, Buoc 1): each path p
% gets assigned a position loc_p = fd_p + 2*N*c1*ell_p in the chirp domain;
% the theorem's guarantee of full diversity (rho=P) requires loc_1 != loc_2
% for AFDM, which in turn requires BOTH paths' Doppler to stay within the
% design bound maxDoppler that c1 was computed for. By choosing ell1,fd1,
% ell2,fd2 explicitly (design bound maxDoppler/guardWidth stays fixed and
% feeds afdmParams.m's c1 exactly as normal), this can construct BOTH a
% "safe" scenario (loc1 != loc2, matches theory doc's own worked example)
% and a "collision" scenario (loc1 == loc2 exactly, by deliberately using a
% Doppler value for one path that exceeds the design bound maxDoppler used
% to compute c1) to empirically confirm the mechanism, rather than relying
% on a velocity sweep's random clutter path to stumble into a collision by
% chance.
%
%   Inputs:
%      - waveform            : 'AFDM', 'OTFS' or 'OFDM'
%      - snrdB                : 1 x numSNR vector of SNR points [dB]
%      - N                    : number of subcarriers
%      - maxDoppler, guardWidth : DESIGN bound fed to afdmParams.m for c1
%                                (may be smaller than |fd1| or |fd2| on
%                                purpose, to construct a guard-violating
%                                scenario)
%      - ell1, ell2           : fixed integer normalized delays of path 1, 2
%      - fd1, fd2             : fixed normalized Doppler shifts of path 1, 2
%                                (same units as maxDoppler)
%      - numTrials            : number of Monte Carlo blocks per SNR point
%                                (only fading h1,h2 vary trial to trial)
%      - progressFcn          : optional progress callback
%   Output:
%      - BER : 1 x numSNR bit error rate curve (BPSK, exhaustive ML, genie CSI)

if nargin < 11 || isempty(progressFcn)
    progressFcn = @(x) [];
end

samplingFreq = 2e6;
numSNR = length(snrdB);
BER = zeros(1, numSNR);

allBits = dec2bin(0:(2^N - 1), N) - '0';
AllX = 1 - 2*allBits.';      % N x 2^N, BPSK candidates in {+1,-1}

totalSteps = numSNR*numTrials;
step = 0;

for k = 1:numSNR
    snrLinear = 10^(snrdB(k)/10);
    bitErrs = 0;
    for t = 1:numTrials
        X = 1 - 2*randi([0 1], N, 1);

        [S, cfg] = waveformMod(waveform, X, N, maxDoppler, guardWidth);

        h1 = (randn+1j*randn)/sqrt(2); h2 = (randn+1j*randn)/sqrt(2);   % equal-strength paths
        if cfg.includeCCP
            H = cconvLTVChannel(N, [ell1 ell2]/samplingFreq, [fd1 fd2]/N*samplingFreq, [h1 h2], samplingFreq, true, cfg.c1);
        else
            H = cconvLTVChannel(N, [ell1 ell2]/samplingFreq, [fd1 fd2]/N*samplingFreq, [h1 h2], samplingFreq, false);
        end

        Sclean = H*S;
        noiseVar = mean(abs(Sclean).^2)/snrLinear;
        R = Sclean + sqrt(noiseVar/2)*(randn(N,1) + 1j*randn(N,1));

        Y = waveformDemod(waveform, R, N, cfg);
        G = waveformChannelMatrix(waveform, H, N, maxDoppler, guardWidth);

        GX = G*AllX;
        dists = sum(abs(Y - GX).^2, 1);
        [~, bestCol] = min(dists);
        xHat = AllX(:, bestCol);

        bitErrs = bitErrs + sum(sign(X) ~= sign(xHat));

        step = step + 1;
        progressFcn(step/totalSteps);
    end
    BER(k) = bitErrs/(N*numTrials);
end

end
