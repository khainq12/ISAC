function Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, includeCCP, c1, dopplerRes)
% Hhat = estimateChannelSinglePath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, includeCCP, c1, dopplerRes)
% Realistic (non-genie) channel estimate from a single known pilot block, reusing the
% matched-filter delay/Doppler estimator (estimateDelayDoppler.m, already validated against
% the CRB in simSensing.m). Only the DOMINANT path is resolved - a low-complexity receiver
% assumption - so on a true P>1-path channel this estimate is deliberately incomplete
% (missing weaker paths entirely), unlike genie CSI which uses the exact H.
%   Inputs:
%      - Rpilot     : N x 1 received pilot samples (through the true channel + noise)
%      - Spilot     : N x 1 known transmitted pilot samples
%      - N          : number of samples
%      - maxDelay   : maximum normalised integer delay index to search
%      - maxDoppler : maximum normalised digital Doppler shift to search
%      - samplingFreq : sampling frequency [Hz]
%      - includeCCP : whether Spilot was transmitted with AFDM's chirp-cyclic prefix
%      - c1         : AFDM chirp frequency (only used if includeCCP)
%      - dopplerRes : Doppler search grid resolution (default 0.02, see estimateDelayDoppler.m)
%   Output:
%      - Hhat : N x N single-path estimated channel matrix (unit path count, may omit
%               weaker true paths entirely - this is the point of "imperfect CSI")

if nargin < 9 || isempty(dopplerRes)
    dopplerRes = 0.02;
end

[ellHat, dopplerHatNorm] = estimateDelayDoppler(Rpilot, Spilot, N, maxDelay, maxDoppler, includeCCP, c1, dopplerRes);

delayHatSec = ellHat/samplingFreq;
dopplerHatHz = dopplerHatNorm*samplingFreq/N;

if includeCCP
    Hmodel = cconvLTVChannel(N, delayHatSec, dopplerHatHz, 1, samplingFreq, true, c1);
else
    Hmodel = cconvLTVChannel(N, delayHatSec, dopplerHatHz, 1, samplingFreq, false);
end

Smodel = Hmodel*Spilot;
gHat = (Smodel'*Rpilot)/(Smodel'*Smodel);   % LS single-tap complex gain estimate

Hhat = gHat*Hmodel;

end
