function [delayIdxHat, dopplerHat] = estimateDelayDopplerRefined(R, S, N, maxDelay, maxDoppler, includeCCP, ccpf, dopplerRes)
% [delayIdxHat, dopplerHat] = estimateDelayDopplerRefined(R, S, N, maxDelay, maxDoppler, includeCCP, ccpf, dopplerRes)
% Same coarse matched-filter grid search as estimateDelayDoppler.m (delayIdxHat is
% still an integer sample index - delay is exactly grid-aligned by construction in this
% codebase's channel model, so no refinement is needed there), PLUS a standard quadratic
% (parabolic) interpolation refinement of the Doppler estimate using the correlation
% magnitude at the coarse peak and its two neighbouring grid points. This is a textbook
% sub-bin peak-refinement technique (e.g. used for FFT-based frequency estimation) - it
% does not change the coarse search at all, only adds a continuous-valued correction on
% top, so it is a strict improvement in resolution, not a different algorithm.
%   Inputs/Outputs: identical to estimateDelayDoppler.m.

if nargin < 8 || isempty(dopplerRes)
    dopplerRes = 0.02;
end

dopplerCandidates = -maxDoppler:dopplerRes:maxDoppler;
w = exp(-1j*2*pi*(0:N-1)'/N);
E = w.^dopplerCandidates;

bestScore = -1;
delayIdxHat = 0;
dopplerHat = 0;
bestCorr = [];
bestIdx = 1;

for ell = 0:maxDelay
    shifted = circshift(S, ell);
    if includeCCP && ell > 0
        ccpVec = [exp(-1j*2*pi*ccpf*((N^2) - 2*N*(ell - (0:(ell-1))))), ones(1, N-ell)];
        shifted = ccpVec(:) .* shifted;
    end
    z = conj(R) .* shifted;
    corr = abs(z.' * E);
    [mx, idx] = max(corr);
    if mx > bestScore
        bestScore = mx;
        delayIdxHat = ell;
        dopplerHat = dopplerCandidates(idx);
        bestCorr = corr;
        bestIdx = idx;
    end
end

% Parabolic sub-bin refinement, only when the peak is strictly interior (not at the
% edge of the search range, where the 3-point formula would need an out-of-range sample)
if bestIdx > 1 && bestIdx < length(bestCorr)
    yL = bestCorr(bestIdx-1); yC = bestCorr(bestIdx); yR = bestCorr(bestIdx+1);
    denom = (yL - 2*yC + yR);
    if abs(denom) > 1e-12
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5); % keep correction within one grid cell
        dopplerHat = dopplerHat + offset*dopplerRes;
    end
end

end
