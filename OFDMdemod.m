function [X] = OFDMdemod(S)
% X = OFDMdemod(S) OFDM demodulation (DFT)
%   Inputs:
%      - S  : N x T sampled sequence of an OFDM modulated signal
%   Output:
%      - X  : N x T matrix of OFDM demodulated samples

sizeS = size(S);
N = sizeS(1);

F = dftmtx(N)/sqrt(N);

X = F*S;

end
