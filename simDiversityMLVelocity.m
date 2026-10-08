function [BER] = simDiversityMLVelocity(waveform, velocitiesKmh, snrdB, N, maxDelay, designMaxDoppler, guardWidth, carrierFreq, samplingFreq, numTrials, csiMode, progressFcn, dopplerMode)
% BER = simDiversityMLVelocity(waveform, velocitiesKmh, snrdB, N, maxDelay, designMaxDoppler, guardWidth, carrierFreq, samplingFreq, numTrials, csiMode, progressFcn, dopplerMode)
% Exhaustive-ML BER of a single waveform versus TARGET VELOCITY (converted to a fixed
% Doppler shift via velocityToDoppler.m, monostatic convention f_D = 2v/lambda), at a
% set of SNR points. Path 1 (the target) has its Doppler fixed to the value implied by
% each tested velocity; path 2 is an independent random-Doppler Rayleigh clutter path
% (equal average strength), same two-path structure as simDiversityML.m. designMaxDoppler
% sizes the waveform's own AFDM guard band / search range and must cover the largest
% tested velocity's normalised Doppler.
%   Inputs:
%      - waveform          : 'AFDM', 'OTFS' or 'OFDM'
%      - velocitiesKmh     : 1 x numVel vector of target velocities [km/h]
%      - snrdB             : 1 x numSNR vector of SNR points [dB]
%      - N                 : number of subcarriers, <= ~20 for tractable exhaustive ML
%      - maxDelay          : maximum normalised integer delay index of either path
%      - designMaxDoppler  : normalised Doppler design bound (AFDM guard band sizing);
%                            must be >= max(|velocityToDoppler(velocitiesKmh)|)
%      - guardWidth        : AFDM guard width
%      - carrierFreq       : carrier frequency [Hz], used for velocity -> Doppler
%      - samplingFreq      : sampling frequency [Hz]
%      - numTrials         : number of Monte Carlo blocks per (velocity, SNR) point
%      - csiMode           : 'genie' (default, true channel used for ML distances) or
%                            'estimated' (dominant path estimated from a known pilot via
%                            estimateChannelSinglePath.m, ML distances use the resulting
%                            incomplete/imperfect Ghat instead)
%      - progressFcn       : optional function handle called with fraction complete in [0,1]
%      - dopplerMode       : 'fractional' (default: target Doppler = true continuous value
%                            implied by the tested velocity, i.e. the physically realistic
%                            case) or 'integer' (target Doppler is rounded to the nearest
%                            integer before building the channel - an idealised comparison
%                            case where the target always lands exactly on the Doppler grid).
%                            Clutter path 2's Doppler is always continuous/random in both modes.
%   Output:
%      - BER : numVel x numSNR bit error rate matrix (BPSK, exhaustive ML)

if nargin < 11 || isempty(csiMode)
    csiMode = 'genie';
end
if nargin < 12 || isempty(progressFcn)
    progressFcn = @(x) [];
end
if nargin < 13 || isempty(dopplerMode)
    dopplerMode = 'fractional';
end

numVel = length(velocitiesKmh);
numSNR = length(snrdB);
BER = zeros(numVel, numSNR);

allBits = dec2bin(0:(2^N - 1), N) - '0';
AllX = 1 - 2*allBits.';      % N x 2^N, BPSK candidates in {+1,-1}

totalSteps = numVel*numSNR*numTrials;
step = 0;

for v = 1:numVel
    [~, fdTarget] = velocityToDoppler(velocitiesKmh(v), carrierFreq, N, samplingFreq);
    if strcmpi(dopplerMode, 'integer')
        fdTarget = round(fdTarget);
    end
    if abs(fdTarget) > designMaxDoppler
        warning("velocityToDoppler: %.0f km/h exceeds designMaxDoppler (%.4f > %.4f); results may alias.", ...
            velocitiesKmh(v), abs(fdTarget), designMaxDoppler);
    end

    for k = 1:numSNR
        snrLinear = 10^(snrdB(k)/10);
        bitErrs = 0;
        for t = 1:numTrials
            X = 1 - 2*randi([0 1], N, 1);

            [S, cfg] = waveformMod(waveform, X, N, designMaxDoppler, guardWidth);

            ell1 = randi([0 maxDelay]); ell2 = randi([0 maxDelay]);
            fd1 = fdTarget;                                  % target: fixed to this velocity
            fd2 = (2*rand-1)*designMaxDoppler;                % clutter: random Doppler
            h1 = (randn+1j*randn)/sqrt(2); h2 = (randn+1j*randn)/sqrt(2);
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
                [Spilot, ~] = waveformMod(waveform, Xpilot, N, designMaxDoppler, guardWidth);
                noiseVarPilot = mean(abs(H*Spilot).^2)/snrLinear;
                Rpilot = H*Spilot + sqrt(noiseVarPilot/2)*(randn(N,1) + 1j*randn(N,1));
                if cfg.includeCCP
                    Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, designMaxDoppler, samplingFreq, true, cfg.c1);
                else
                    Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, designMaxDoppler, samplingFreq, false, 0);
                end
                G = waveformChannelMatrix(waveform, Hhat, N, designMaxDoppler, guardWidth);
            else
                G = waveformChannelMatrix(waveform, H, N, designMaxDoppler, guardWidth);
            end

            GX = G*AllX;
            dists = sum(abs(Y - GX).^2, 1);
            [~, bestCol] = min(dists);
            xHat = AllX(:, bestCol);

            bitErrs = bitErrs + sum(sign(X) ~= sign(xHat));

            step = step + 1;
            progressFcn(step/totalSteps);
        end
        BER(v,k) = bitErrs/(N*numTrials);
    end
end

end
