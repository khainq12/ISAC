% Recompute ONLY the fractional-case low-complexity-MMSE curve (detector 3)
% for fig10b with a larger per-row support Lband, keeping the already-good
% MMSE/DFE curves and the integer-case results untouched.
% Root cause: paper's L=(2*xi_nu+1)*P formula assumes a null-guard-truncated
% Heff (no circular wraparound, see paper Fig 4/5 + eq 52-53). This script's
% fig10_detectors.m does NOT insert null guard symbols (transmits full N
% data symbols), so the effective channel is still circular and needs a
% larger empirical support to capture equivalent energy - confirmed via
% diag_fig10b_lowc_L.m: at L=9 AFDM-lowC BER=1.29e-2 > OFDM-lowC=6.98e-3
% (wrong order vs paper); crossover between L=15/21; at L=27 AFDM-lowC=
% 1.52e-3 < OFDM-lowC=2.15e-3 (correct order, AFDM better, matches paper's
% Fig 10b where AFDM,low-complexity-MMSE sits below all OFDM curves).
clear; clc; close all;
rng(9);

load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig10_results.mat');

N = 256; P = 3; lmax = 2; alpha_max = 2; xi_nu = 1;
c1_afdm = (2*(alpha_max+xi_nu)+1)/(2*N);
U_afdm = daft(N, c1_afdm, 1/(2*N^2))';
U_ofdm = dftmtx(N)'/sqrt(N);
Lband = 3*(2*xi_nu+1)*P; % 27, see note above

SNR_dB = results.fractional.SNR_dB;
numMonteLocal = numMonte;
BER3 = zeros(2, numel(SNR_dB));

for isnr = 1:numel(SNR_dB)
    snr_lin = 10^(SNR_dB(isnr)/10);
    N0 = 1/snr_lin;
    errB = zeros(1,2); totB = zeros(1,2);
    for imc = 1:numMonteLocal
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        delays = 0:lmax;
        theta = -pi + 2*pi*rand(1,P);
        doppler_norm = alpha_max * cos(theta);

        bits = randi([0 1], N*2, 1);
        sym = qammod(bits, 4, 'gray', 'InputType','bit', 'UnitAveragePower', true);
        noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));

        for w = 1:2
            if w==1, Uw = U_afdm; c1w = c1_afdm; else, Uw = U_ofdm; c1w = 0; end
            Hw = Gen_channel_mtx(N, P, gains, delays, doppler_norm/N, c1w);
            s = Uw * sym;
            r = Hw*s + noise;
            Heff = Uw' * Hw * Uw;
            z = Uw' * r;

            HeffB = topRowMagnitudeMask(Heff, Lband);
            xHat = (HeffB'*HeffB + N0*eye(N)) \ (HeffB' * z);
            bitsHat = qamdemod(xHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
            errB(w) = errB(w) + sum(bitsHat ~= bits);
            totB(w) = totB(w) + numel(bits);
        end
    end
    BER3(:,isnr) = errB ./ totB;
    fprintf('[fractional-lowC-fix] SNR=%4.1f: AFDM-lowC=%.4e  OFDM-lowC=%.4e\n', SNR_dB(isnr), BER3(1,isnr), BER3(2,isnr));
end

results.fractional.BER(3,1,:) = reshape(BER3(1,:),1,1,[]);
results.fractional.BER(3,2,:) = reshape(BER3(2,:),1,1,[]);

save('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig10_results.mat', 'results','numMonte');

detNames = {'MMSE','MRC-DFE','low-complexity MMSE'};
waveNames = {'AFDM','OFDM'};
mk = {'-o','-s','-d'};
clr = {[0 0.30 0.65],[0.75 0.20 0.05]};
LW = 3; MS = 7;

figure('Color','w','Position',[100 100 1300 520]);
subplot(1,2,1);
BERa = floorZeros(results.integer.BER); SNRa = results.integer.SNR_dB;
for w=1:2
    for d=1:3
        semilogy(SNRa, squeeze(BERa(d,w,:)), mk{d}, 'Color', clr{w}, 'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', sprintf('%s, %s', waveNames{w}, detNames{d})); hold on;
    end
end
grid on; box on; set(gca,'FontSize',11,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
title('Fig 10a repro: integer Doppler','FontWeight','bold'); legend('show','FontSize',8);

subplot(1,2,2);
BERb = floorZeros(results.fractional.BER); SNRb = results.fractional.SNR_dB;
for w=1:2
    for d=1:3
        semilogy(SNRb, squeeze(BERb(d,w,:)), mk{d}, 'Color', clr{w}, 'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', sprintf('%s, %s', waveNames{w}, detNames{d})); hold on;
    end
end
grid on; box on; set(gca,'FontSize',11,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
title('Fig 10b repro: fractional Doppler','FontWeight','bold'); legend('show','FontSize',8);

exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig10_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function Bc = floorZeros(B)
Bc = B;
nz = B(B>0);
floorVal = min(nz(:));
Bc(B==0) = floorVal;
end

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end

function Hb = topRowMagnitudeMask(H, L)
N = size(H,1);
Hb = zeros(size(H));
[~, ord] = sort(abs(H), 2, 'descend');
keep = ord(:, 1:min(L, N));
for r = 1:N
    Hb(r, keep(r,:)) = H(r, keep(r,:));
end
end
