% Reproduction of Fig. 8b in Bemani, Ksairi, Kountouris (arXiv:2204.12798).
% OFDM/OCDM/OTFS/AFDM, BPSK, N=16, ML detection, 3-path channel (fixed
% delays 0,1,2 i.e. lmax=2), alpha_max=1, Jakes Doppler spectrum (fractional).
clear; clc; close all;
rng(4);

N = 16;
P = 3; lmax = 2;
alpha_max = 1;
xi_nu = 1;
c1_afdm = (2*(alpha_max+xi_nu)+1)/(2*N);
c2_afdm = 1/(N^2);
c1_ocdm = 1/(2*N);
c2_ocdm = 1/(2*N);
NOTFS = 4; MOTFS = 4;

delays_fixed = 0:lmax; % fixed delay taps as stated in the paper text for Fig 8b

SNR_dB = 0:2:20;
numMonte = 3000;

% U_list stores the TX (modulation) operator directly: s = U*sym.
F = dftmtx(N)/sqrt(N);
U_ofdm = F';                            % IFFT, already TX-form
U_ocdm = daft(N, c1_ocdm, c2_ocdm)';    % A' = IDAFT
U_afdm = daft(N, c1_afdm, c2_afdm)';    % A' = IDAFT
FN = dftmtx(NOTFS)/sqrt(NOTFS);
U_otfs = kron(FN', eye(MOTFS));         % already TX-form (Zak/block operator)
waveNames = {'OFDM','OCDM','OTFS','AFDM'};
U_list = {U_ofdm, U_ocdm, U_otfs, U_afdm};
c1_list = [0, c1_ocdm, 0, c1_afdm];

nHyp = 2^N;
S_hyp = zeros(N, nHyp);
for idx = 0:nHyp-1
    bits = bitget(idx, N:-1:1);
    S_hyp(:,idx+1) = 1 - 2*bits;
end
X_hyp = cell(4,1);
for w=1:4
    X_hyp{w} = U_list{w} * S_hyp;
end

BER = zeros(4, numel(SNR_dB));

for isnr = 1:numel(SNR_dB)
    snr_lin = 10^(SNR_dB(isnr)/10);
    N0 = 1/snr_lin;
    errBits = zeros(4,1); totBits = zeros(4,1);
    for imc = 1:numMonte
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        theta = -pi + 2*pi*rand(1,P);
        doppler_norm = alpha_max * cos(theta); % Jakes spectrum, fractional

        xBits = randi([0 1], N, 1);
        xSym = 1 - 2*xBits;
        noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));

        for w = 1:4
            Uw = U_list{w};
            Hw = Gen_channel_mtx(N, P, gains, delays_fixed, doppler_norm/N, c1_list(w));
            s = Uw * xSym;
            r = Hw*s + noise;
            HX = Hw * X_hyp{w};
            dist2 = sum(abs(r - HX).^2, 1);
            [~, kBest] = min(dist2);
            xHatBits = bitget(kBest-1, N:-1:1).';
            errBits(w) = errBits(w) + sum(xHatBits ~= xBits);
            totBits(w) = totBits(w) + N;
        end
    end
    BER(:,isnr) = errBits ./ totBits;
    fprintf('SNR=%2d dB: OFDM=%.4e OCDM=%.4e OTFS=%.4e AFDM=%.4e\n', SNR_dB(isnr), BER(1,isnr), BER(2,isnr), BER(3,isnr), BER(4,isnr));
end

save('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8b_results.mat', 'BER', 'SNR_dB', 'waveNames', 'numMonte');

figure('Color','w');
mk = {'-o','-s','-d','-^'};
for w=1:4
    semilogy(SNR_dB, BER(w,:), mk{w}, 'DisplayName', waveNames{w}); hold on;
end
grid on; xlabel('SNR in dB'); ylabel('BER');
title(sprintf('Fig 8b reproduction: ML, N=%d, numMonte=%d', N, numMonte));
legend('show'); ylim([1e-5 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8b_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
