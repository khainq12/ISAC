function [crbDelay, crbDoppler] = crbDelayDoppler(S, samplingFreq, snrLinear)
% [crbDelay, crbDoppler] = crbDelayDoppler(S, samplingFreq, snrLinear)
% Cramer-Rao bound for delay [s] and Doppler [Hz] estimation of a single unknown
% complex-amplitude path in additive white Gaussian noise, given the known transmit
% samples S (classical matched-filter/ambiguity-function CRB, e.g. Kay, "Fundamentals
% of Statistical Signal Processing: Estimation Theory", using the block SNR = Es/N0
% definition consistent with awgn(...,'measured').
%   Inputs:
%      - S            : N x 1 known transmit samples
%      - samplingFreq : sampling frequency [Hz]
%      - snrLinear    : linear SNR (Es/N0), scalar or vector
%   Outputs:
%      - crbDelay   : CRB on delay variance [s^2], same size as snrLinear
%      - crbDoppler : CRB on Doppler variance [Hz^2], same size as snrLinear

N = length(S);
t = (0:N-1)'/samplingFreq;
Es = sum(abs(S).^2);

tBar = sum(t .* abs(S).^2)/Es;
T2 = sum(((t - tBar).^2) .* abs(S).^2)/Es;         % mean-square duration [s^2]

Sf = fftshift(fft(S));
f = (-N/2:N/2-1)'/N*samplingFreq;
Ef = sum(abs(Sf).^2);
fBar = sum(f .* abs(Sf).^2)/Ef;
beta2 = sum(((f - fBar).^2) .* abs(Sf).^2)/Ef;     % mean-square bandwidth [Hz^2]

crbDelay = 1./(8*pi^2*beta2*N*snrLinear);
crbDoppler = 1./(8*pi^2*T2*N*snrLinear);

end
