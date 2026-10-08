function [BER] = simDiversityML(waveform, snrdB, N, maxDelay, maxDoppler, guardWidth, numTrials, csiMode, progressFcn)
% BER = simDiversityML(waveform, snrdB, N, maxDelay, maxDoppler, guardWidth, numTrials, csiMode, progressFcn)
% Small-scale exhaustive-ML BER demo isolating each waveform's diversity order under a
% genuine 2-path equal-strength Rayleigh doubly-dispersive channel. Unlike zero-forcing
% (which is provably waveform-invariant, see simBERSER.m), Maximum Likelihood detection
% achieves the true diversity order supported by each waveform's structure, matching the
% "AFDM achieves guaranteed full diversity" claim in the literature (e.g. Bemani et al.,
% cited by arXiv:2401.07700 Table I). BPSK is used and N kept small (<=16) so exhaustive
% ML search (2^N candidates) stays computationally tractable.
%   Inputs:
%      - waveform     : 'AFDM', 'OTFS' or 'OFDM'
%      - snrdB        : 1 x numSNR vector of SNR points [dB]
%      - N            : number of subcarriers, <= ~20 for tractable exhaustive ML
%                       (perfect square if waveform is 'OTFS')
%      - maxDelay     : maximum normalised integer delay index of either path
%      - maxDoppler   : maximum normalised digital Doppler shift of either path
%      - guardWidth   : AFDM guard width
%      - numTrials    : number of Monte Carlo blocks per SNR point
%      - csiMode      : 'genie' (default, ML distances computed against the true channel
%                       matrix G) or 'estimated' (realistic: dominant path estimated from
%                       a known pilot via estimateChannelSinglePath.m, ML distances use the
%                       resulting incomplete/imperfect Ghat instead)
%      - progressFcn  : optional function handle called with fraction complete in [0,1]
%   Output:
%      - BER : 1 x numSNR bit error rate curve (BPSK, exhaustive ML)

if nargin < 8 || isempty(csiMode)
    csiMode = 'genie';
end
if nargin < 9 || isempty(progressFcn)
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

        ell1 = randi([0 maxDelay]); ell2 = randi([0 maxDelay]);
        fd1 = (2*rand-1)*maxDoppler; fd2 = (2*rand-1)*maxDoppler;
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

        if strcmpi(csiMode, 'estimated')
            Xpilot = ones(N, 1);
            [Spilot, ~] = waveformMod(waveform, Xpilot, N, maxDoppler, guardWidth);
            noiseVarPilot = mean(abs(H*Spilot).^2)/snrLinear;
            Rpilot = H*Spilot + sqrt(noiseVarPilot/2)*(randn(N,1) + 1j*randn(N,1));
            if cfg.includeCCP
                Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, true, cfg.c1);
            else
                Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, false, 0);
            end
            G = waveformChannelMatrix(waveform, Hhat, N, maxDoppler, guardWidth);
        else
            G = waveformChannelMatrix(waveform, H, N, maxDoppler, guardWidth);
        end

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
