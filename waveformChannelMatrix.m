function [G] = waveformChannelMatrix(waveform, H, N, maxDoppler, guardWidth)
% G = waveformChannelMatrix(waveform, H, N, maxDoppler, guardWidth)
% Effective symbol-domain channel matrix G = Demod * H * Modulate, i.e. the same
% "effective channel" matrix visualised by app_DDChannel.mlapp. This is where
% AFDM/OTFS/OFDM differ structurally (AFDM stays concentrated near the diagonal,
% OTFS Doppler-leaks across parallel diagonals, OFDM ICI-smears near the diagonal)
% -- the raw sample-domain channel H itself is waveform-independent and always
% maximally sparse (single nonzero per row), so structural differences only show
% up after transforming into each waveform's own symbol domain.
%   Inputs:
%      - waveform   : 'AFDM', 'OTFS' or 'OFDM'
%      - H          : N x N sample-domain channel matrix
%      - N          : number of subcarriers
%      - maxDoppler : maximum normalised digital Doppler shift (AFDM only)
%      - guardWidth : AFDM guard width (AFDM only)
%   Output:
%      - G : N x N effective symbol-domain channel matrix

F = dftmtx(N)/sqrt(N);

switch upper(waveform)
    case 'AFDM'
        [c1, c2] = afdmParams(N, maxDoppler, guardWidth);
        L1 = diag(exp(-1j*2*pi*c1*((0:N-1).^2)));
        L2 = diag(exp(-1j*2*pi*c2*((0:N-1).^2)));
        G = L2*F*L1 * H * L1'*F'*L2';
    case 'OTFS'
        [K, L] = otfsGridKL(N);
        F_L = dftmtx(L)/sqrt(L);
        G = kron(F_L, eye(K)) * H * kron(F_L', eye(K));
    case 'OFDM'
        G = F*H*F';
    otherwise
        error("Unknown waveform: %s", waveform);
end

end
