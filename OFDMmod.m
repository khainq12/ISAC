function [S] = OFDMmod(X)
% S = OFDMmod(X) OFDM modulation (IDFT)
%   Inputs:
%      - X  : N x T transmit symbol matrix
%   Output:
%      - S  : N x T matrix of OFDM modulated samples

sizeX = size(X);
N = sizeX(1);

F = dftmtx(N)/sqrt(N);

S = F'*X;

end
