function [BER, SER] = simBERSER(waveform, channelType, snrdB, N, bps, P, maxDelay, maxDoppler, guardWidth, samplingFreq, numTrials, equalizerMode, numTaps, csiMode, progressFcn)
% [BER, SER] = simBERSER(waveform, channelType, snrdB, N, bps, P, maxDelay, maxDoppler, guardWidth, samplingFreq, numTrials, equalizerMode, numTaps, csiMode, progressFcn)
% Monte Carlo BER/SER of AFDM/OTFS/OFDM over an AWGN or Rayleigh doubly-dispersive channel.
% Equalization happens in the symbol domain against the effective channel matrix
% G = Demod*H*Modulate (see waveformChannelMatrix.m) - this is the matrix whose
% structure (concentrated for AFDM, Doppler-leaked for OTFS, ICI-smeared for OFDM)
% actually differs by waveform, unlike the raw sample-domain channel H.
%   Inputs:
%      - waveform     : 'AFDM', 'OTFS' or 'OFDM'
%      - channelType  : 'AWGN' or 'Rayleigh'
%      - snrdB        : 1 x numSNR vector of SNR points [dB]
%      - N            : number of subcarriers (perfect square if waveform is 'OTFS')
%      - bps          : bits per QAM symbol
%      - P            : number of channel paths (P == 1 for 'AWGN')
%      - maxDelay     : maximum normalised integer delay index of any path
%      - maxDoppler   : maximum normalised digital Doppler shift of any path
%      - guardWidth   : AFDM guard width
%      - samplingFreq : sampling frequency [Hz]
%      - numTrials    : number of Monte Carlo blocks per SNR point
%      - equalizerMode: 'ideal' (full zero-forcing on G, perfect CSI - waveform-invariant
%                       BER by construction) or 'banded' (fixed tap-budget receiver that
%                       only keeps the numTaps strongest entries of G per row -
%                       reveals each waveform's energy concentration in its own domain)
%      - numTaps      : number of strongest per-row entries of G kept (banded mode only)
%      - csiMode      : 'genie' (default, equalize against the true channel matrix G) or
%                       'estimated' (realistic: a known pilot block is sent through the
%                       same channel realization, the dominant path is estimated via
%                       estimateChannelSinglePath.m - a single-path, low-complexity
%                       receiver assumption - and the resulting incomplete/imperfect
%                       Ghat is used for equalization instead of the true G)
%      - progressFcn  : optional function handle called with fraction complete in [0,1]
%   Outputs:
%      - BER, SER : 1 x numSNR bit/symbol error rate curves

if nargin < 14 || isempty(csiMode)
    csiMode = 'genie';
end
if nargin < 15 || isempty(progressFcn)
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
            H = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, true, true, cfg.c1);
        else
            H = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, false, true);
        end

        Sclean = H*S;
        noiseVar = mean(abs(Sclean).^2)/snrLinear;
        R = Sclean + sqrt(noiseVar/2)*(randn(N,1) + 1j*randn(N,1));

        Y = waveformDemod(waveform, R, N, cfg);            % linear demod only, no equalization yet

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
