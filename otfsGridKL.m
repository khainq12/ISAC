function [K, L] = otfsGridKL(N)
% [K, L] = otfsGridKL(N) Square K x L delay-Doppler grid (K = L = sqrt(N)) for OTFS.
%   Input:
%      - N : number of subcarriers, must be a perfect square
%   Output:
%      - K, L : OTFS grid dimensions, K*L = N

L = sqrt(N);
if L ~= floor(L)
    error("N must be a perfect square for OTFS (e.g. 16, 36, 64, 100, 144, 256).");
end
K = L;

end
