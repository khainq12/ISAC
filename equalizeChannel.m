function [S_hat] = equalizeChannel(H, R, mode, numTaps, noiseVar)
% S_hat = equalizeChannel(H, R, mode, numTaps, noiseVar)
% Channel equalization for the received vector R = H*S + n.
%   Inputs:
%      - H        : N x N channel matrix (known, perfect CSI)
%      - R        : N x 1 received samples
%      - mode     : 'ideal' (full zero-forcing) or 'banded' (low-complexity,
%                   receiver only keeps the numTaps largest-magnitude entries per
%                   row, i.e. a fixed tap budget shared across waveforms)
%      - numTaps  : number of strongest entries kept per row (banded mode only)
%      - noiseVar : per-sample noise variance, used for MMSE regularisation (banded mode only)
%   Output:
%      - S_hat : N x 1 equalized sample estimate
%
% 'ideal' full ZF is unitarily invariant to waveform choice (same BER for
% AFDM/OTFS/OFDM on any channel), since it exactly inverts H regardless of how
% concentrated/spread its energy is. 'banded' truncates each row of H to only
% its numTaps strongest entries before inverting (matched-filter/message-passing
% style low-complexity receivers, as commonly used in practice) - since a path's
% energy sits on a delay-offset diagonal (not necessarily the main one), truncating
% by magnitude rather than a fixed position around i==j is what actually captures
% "how concentrated is this waveform's channel representation" - this is where a
% waveform's ability to keep a path's energy concentrated in few taps (AFDM by
% design) versus spreading it out (OTFS Doppler leakage, OFDM ICI) shows up in BER/SER.

switch mode
    case 'ideal'
        S_hat = H\R;
    case 'banded'
        N = size(H, 1);
        Hb = zeros(N);
        for i = 1:N
            [~, idx] = maxk(abs(H(i,:)), numTaps);
            Hb(i, idx) = H(i, idx);
        end
        S_hat = (Hb'*Hb + noiseVar*eye(N)) \ (Hb'*R);
    otherwise
        error("Unknown equalizer mode: %s", mode);
end

end
