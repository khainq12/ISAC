% Reproduction of Fig. 8a in Bemani, Ksairi, Kountouris,
% "Affine Frequency Division Multiplexing for Next Generation Wireless
% Communications," IEEE TWC 2023 (arXiv:2204.12798).
% AFDM, BPSK, N=16, ML detection, P=2/3/4 paths, integer Doppler, alpha_max=1.
clear; clc; close all;
rng(1);

N = 16;
alpha_max = 1;
c1 = (2*alpha_max+1)/(2*N);
c2 = 1/(N^2);

SNR_dB = 0:2:20;
Pvals = [2 3 4];
numMonte = getenv_num('NUMMONTE', 10000);

S_all = [-1 1];
% precompute all 2^N BPSK hypotheses as columns of S_hyp (N x 2^N)
nHyp = 2^N;
S_hyp = zeros(N, nHyp, 'like', 1+1i);
for idx = 0:nHyp-1
    bits = bitget(idx, N:-1:1);
    S_hyp(:,idx+1) = 1 - 2*bits; % 0->+1, 1->-1
end
% also compute this via DAFT domain: transmit AFDM-modulated bit hypotheses
A = daft_matrix(N, c1, c2);
X_hyp = A' * S_hyp;  % time-domain AFDM signal for each hypothesis (N x nHyp)

BER = zeros(numel(Pvals), numel(SNR_dB));

for ip = 1:numel(Pvals)
    P = Pvals(ip);
    lmax = P-1;
    for isnr = 1:numel(SNR_dB)
        snr_lin = 10^(SNR_dB(isnr)/10);
        N0 = 1/snr_lin;  % unit tx power per symbol (BPSK +-1)
        errBits = 0; totBits = 0;
        for imc = 1:numMonte
            gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
            delays = randperm(lmax+1, min(P, lmax+1)) - 1;
            if numel(delays) < P
                delays = [delays, randi([0 lmax], 1, P-numel(delays))];
            end
            dopplers_int = randi([-alpha_max alpha_max], 1, P);

            H = Gen_channel_mtx(N, P, gains, delays, dopplers_int/N, c1);

            xBits = randi([0 1], N, 1);
            xSym = 1 - 2*xBits;
            s = A' * xSym;

            noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));
            r = H*s + noise;

            % ML: minimize || r - H*X_hyp(:,k) ||^2 over all hypotheses
            HX = H * X_hyp; % N x nHyp
            dist2 = sum(abs(r - HX).^2, 1);
            [~, kBest] = min(dist2);
            xHatBits = bitget(kBest-1, N:-1:1).';

            errBits = errBits + sum(xHatBits ~= xBits);
            totBits = totBits + N;
        end
        BER(ip, isnr) = errBits / totBits;
        fprintf('P=%d SNR=%2d dB: BER=%.4e (numMonte=%d)\n', P, SNR_dB(isnr), BER(ip,isnr), numMonte);
    end
end

save('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8a_results.mat', 'BER', 'SNR_dB', 'Pvals', 'numMonte');

figure('Color','w');
semilogy(SNR_dB, BER(1,:), 'o-', 'DisplayName','P = 2'); hold on;
semilogy(SNR_dB, BER(2,:), 's-', 'DisplayName','P = 3');
semilogy(SNR_dB, BER(3,:), '^-', 'DisplayName','P = 4');
grid on; xlabel('SNR in dB'); ylabel('BER');
title(sprintf('Fig 8a reproduction: AFDM diversity (N=%d, numMonte=%d)', N, numMonte));
legend('show'); ylim([1e-6 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8a_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function A = daft_matrix(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end

function v = getenv_num(name, default)
s = getenv(name);
if isempty(s)
    v = default;
else
    v = str2double(s);
end
end
