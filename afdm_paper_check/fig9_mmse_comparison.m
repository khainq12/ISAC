% Reproduction of Fig. 9a/9b in Bemani, Ksairi, Kountouris, "Affine
% Frequency Division Multiplexing for Next Generation Wireless
% Communications," IEEE TWC 2023 (arXiv:2204.12798).
% OFDM/OCDM/OTFS/AFDM, QPSK, N=256, MMSE detection, 3-path channel,
% lmax=2, alpha_max=2, Jakes Doppler spectrum.
clear; clc; close all;
rng(2);

N = 256;
NOTFS = 16; MOTFS = 16; % NOTFS*MOTFS = N
P = 3; lmax = 2;
alpha_max = 2;
xi_nu = 1; % guard for fractional Doppler (paper eq. 43), assumed minimal value

c1_afdm = (2*(alpha_max+xi_nu)+1)/(2*N);
c2_afdm = 1/(2*N^2);
c1_ocdm = 1/(2*N);
c2_ocdm = 1/(2*N);

SNR_dB = 0:2.5:30;
numMonte = 3000;

F = dftmtx(N)/sqrt(N);
n = (0:N-1).';

% U_list stores the TX (modulation) operator directly: s = U*sym.
% For OCDM/AFDM the paper's analysis matrix is A = L2*F*L1 and modulation
% is x = A'*X (IDAFT), so the TX operator is A' (Hermitian transpose).
U_ofdm = F'; % IFFT
U_ocdm = daft(N, c1_ocdm, c2_ocdm)';
U_afdm = daft(N, c1_afdm, c2_afdm)';
FN = dftmtx(NOTFS)/sqrt(NOTFS);
U_otfs = kron(FN', eye(MOTFS)); % already in TX form (Zak/block operator)

waveNames = {'OFDM','OCDM','OTFS','AFDM'};
U_list = {U_ofdm, U_ocdm, U_otfs, U_afdm};

BER = zeros(4, numel(SNR_dB));

for isnr = 1:numel(SNR_dB)
    snr_lin = 10^(SNR_dB(isnr)/10);
    N0 = 1/snr_lin; % unit avg symbol power (QPSK normalized)
    errBits = zeros(4,1); totBits = zeros(4,1);
    for imc = 1:numMonte
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        delays = randperm(lmax+1, P) - 1;
        theta = -pi + 2*pi*rand(1,P);
        doppler_norm = alpha_max * cos(theta); % fractional Doppler, Jakes spectrum

        % NOTE: the CPP correction (eq. 23) depends on each waveform's own
        % c1, so H is regenerated per waveform inside the loop below.
        bits = randi([0 1], N*2, 1);
        sym = qammod(bits, 4, 'gray', 'InputType','bit', 'UnitAveragePower', true);

        noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));

        for w = 1:4
            switch w
                case 1, c1w = 0;         % OFDM: no CPP needed (2N*c1 integer trivially with c1=0)
                case 2, c1w = c1_ocdm;
                case 3, c1w = 0;         % OTFS block operator: no CPP phase (uses CP only)
                case 4, c1w = c1_afdm;
            end
            Hw = Gen_channel_mtx(N, P, gains, delays, doppler_norm/N, c1w);
            Uw = U_list{w}; % TX operator; RX operator is Uw' (unitary)
            s = Uw * sym;
            r = Hw*s + noise;
            Heff = Uw' * Hw * Uw;
            symHat = (Heff'*Heff + N0*eye(N)) \ (Heff' * (Uw'*r));
            bitsHat = qamdemod(symHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
            errBits(w) = errBits(w) + sum(bitsHat ~= bits);
            totBits(w) = totBits(w) + numel(bits);
        end
    end
    BER(:,isnr) = errBits ./ totBits;
    fprintf('SNR=%5.1f dB : OFDM=%.4e OCDM=%.4e OTFS=%.4e AFDM=%.4e\n', SNR_dB(isnr), BER(1,isnr), BER(2,isnr), BER(3,isnr), BER(4,isnr));
end

save('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig9a_results.mat', 'BER', 'SNR_dB', 'waveNames', 'numMonte');

figure('Color','w');
mk = {'-o','-s','-d','-^'};
for w=1:4
    semilogy(SNR_dB, BER(w,:), mk{w}, 'DisplayName', waveNames{w}); hold on;
end
grid on; xlabel('SNR in dB'); ylabel('BER');
title(sprintf('Fig 9a reproduction: MMSE, N=%d, numMonte=%d', N, numMonte));
legend('show'); ylim([1e-6 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig9a_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
