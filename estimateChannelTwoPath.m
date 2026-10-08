function Hhat = estimateChannelTwoPath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, includeCCP, c1, dopplerRes)
% Hhat = estimateChannelTwoPath(Rpilot, Spilot, N, maxDelay, maxDoppler, samplingFreq, includeCCP, c1, dopplerRes)
% Two-path channel estimate from a single known pilot block, via successive
% interference cancellation (SIC) on top of the existing single-path matched-filter
% estimator (estimateDelayDoppler.m, estimateChannelSinglePath.m): estimate the
% dominant path first, reconstruct and subtract its pilot contribution, then estimate
% the residual (second) path from what remains. Unlike estimateChannelSinglePath.m,
% this does NOT deliberately omit the weaker path, so at high SNR the estimate
% converges toward the true P=2 channel (and toward genie-CSI performance).
%   Inputs/outputs: same convention as estimateChannelSinglePath.m.

if nargin < 9 || isempty(dopplerRes)
    dopplerRes = 0.02;
end

% Path 1 (dominant)
[ell1, dopplerNorm1] = estimateDelayDoppler(Rpilot, Spilot, N, maxDelay, maxDoppler, includeCCP, c1, dopplerRes);
delaySec1 = ell1/samplingFreq;
dopplerHz1 = dopplerNorm1*samplingFreq/N;
if includeCCP
    Hmodel1 = cconvLTVChannel(N, delaySec1, dopplerHz1, 1, samplingFreq, true, c1);
else
    Hmodel1 = cconvLTVChannel(N, delaySec1, dopplerHz1, 1, samplingFreq, false);
end
Smodel1 = Hmodel1*Spilot;
g1 = (Smodel1'*Rpilot)/(Smodel1'*Smodel1);

% Cancel path 1's reconstructed contribution, then estimate the residual path
Rresidual = Rpilot - g1*Smodel1;
[ell2, dopplerNorm2] = estimateDelayDoppler(Rresidual, Spilot, N, maxDelay, maxDoppler, includeCCP, c1, dopplerRes);
delaySec2 = ell2/samplingFreq;
dopplerHz2 = dopplerNorm2*samplingFreq/N;
if includeCCP
    Hmodel2 = cconvLTVChannel(N, delaySec2, dopplerHz2, 1, samplingFreq, true, c1);
else
    Hmodel2 = cconvLTVChannel(N, delaySec2, dopplerHz2, 1, samplingFreq, false);
end
Smodel2 = Hmodel2*Spilot;
g2 = (Smodel2'*Rresidual)/(Smodel2'*Smodel2);

Hhat = g1*Hmodel1 + g2*Hmodel2;

end
