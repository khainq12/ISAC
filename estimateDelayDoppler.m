function [delayIdxHat, dopplerHat] = estimateDelayDoppler(R, S, N, maxDelay, maxDoppler, includeCCP, ccpf, dopplerRes)
% [delayIdxHat, dopplerHat] = estimateDelayDoppler(R, S, N, maxDelay, maxDoppler, includeCCP, ccpf, dopplerRes)
% Matched-filter (cross-ambiguity) grid search for the dominant path's normalised
% integer delay index and digital Doppler shift, given the known transmit samples S.
%   Inputs:
%      - R          : N x 1 received samples
%      - S          : N x 1 known transmit samples
%      - N          : number of samples
%      - maxDelay   : maximum normalised integer delay index to search
%      - maxDoppler : maximum normalised digital Doppler shift to search
%      - includeCCP : whether S was transmitted with AFDM's chirp-cyclic prefix
%      - ccpf       : AFDM chirp frequency c1 (only used if includeCCP)
%      - dopplerRes : Doppler search grid resolution (default 0.02)
%   Outputs:
%      - delayIdxHat : estimated normalised integer delay index
%      - dopplerHat  : estimated normalised digital Doppler shift (can be fractional)

if nargin < 8 || isempty(dopplerRes)
    dopplerRes = 0.02;
end

dopplerCandidates = -maxDoppler:dopplerRes:maxDoppler;
w = exp(-1j*2*pi*(0:N-1)'/N);
E = w.^dopplerCandidates;  % N x numDoppler, matches cconvLTVChannel's W^nu (principal branch for fractional nu)

bestScore = -1;
delayIdxHat = 0;
dopplerHat = 0;

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
    end
end

end
