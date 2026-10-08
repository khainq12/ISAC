function [H, H_p] = cconvLTVChannelFast(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, includeCCP, varargin)
% [H, H_p] = cconvLTVChannelFast(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, includeCCP, [c1])
% Numerically IDENTICAL to cconvLTVChannel.m (validated against it in
% validate_cconvLTVChannelFast.m across many random cases, integer and fractional
% Doppler, with/without CCP), but O(N) per path instead of O(N^3): the original builds
% W^nu and Pi^ell as full N x N matrices and calls generic mpower/matrix-multiply on
% them; since W is exactly diagonal and Pi is exactly a permutation, this version builds
% them as sparse matrices directly (W^nu via elementwise scalar power on the diagonal -
% mathematically identical to mpower's eigendecomposition route for a diagonal matrix,
% same principal-branch convention for fractional nu; Pi^ell via direct index
% construction, exact for integer ell) and multiplies them as sparse, then returns a
% full N x N H for drop-in compatibility with every other function in this codebase.
%   Inputs/Outputs: identical to cconvLTVChannel.m.

if length(pathDelays) == length(pathDopplers) && length(pathDopplers) == length(fadingCoefs)
    P = length(pathDelays);
else
    error("Invalid number of target parameters");
end

normDelays = round(pathDelays*samplingFreq);
normDopplers = N*pathDopplers/samplingFreq;

n = (0:N-1)';
wBase = exp(-1j*2*pi*n/N);   % diagonal entries of W

if includeCCP
    ccpf = varargin{1};
else
    ccpf = 0;
end

H = sparse(N, N);
H_p = zeros(N, N, P);
for p = 1:P
    ell = normDelays(p);
    cppDiag = [exp(-1j*2*pi*ccpf*((N^2)-2*N*(ell - (0:(ell-1))))), ones(1, N - ell)];
    CPP_p = spdiags(cppDiag(:), 0, N, N);

    wnuDiag = wBase.^normDopplers(p);   % scalar elementwise power, same principal branch as mpower on a diagonal matrix
    Wnu = spdiags(wnuDiag, 0, N, N);

    % Pi^ell as an exact sparse permutation matrix (ell is always a non-negative
    % integer here, so this is an exact identity, not an approximation)
    rowIdx = (1:N)';
    colIdx = mod(rowIdx - 1 - ell, N) + 1;
    Piell = sparse(rowIdx, colIdx, 1, N, N);

    Hp_sparse = fadingCoefs(p) * (CPP_p * Wnu * Piell);
    H = H + Hp_sparse;
    if nargout > 1
        H_p(:,:,p) = full(Hp_sparse);
    end
end
H = full(H);

end
