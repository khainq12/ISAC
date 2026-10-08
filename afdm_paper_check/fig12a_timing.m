% Reproduction of Fig. 12a in Bemani, Ksairi, Kountouris (arXiv:2204.12798).
% AFDM, N=256, embedded-pilot channel estimation, FRACTIONAL Doppler case,
% BER vs SNRd for Ideal CSI and SNRp = 30/35/40 dB.
clear; clc; close all;
rng(6);

N = 256; P = 3; lmax = 2; alpha_max = 2; xi_nu = 1; % fractional case
c1 = (2*(alpha_max+xi_nu)+1)/(2*N);
c2 = 1/(2*N^2);
Q = (lmax+1)*(2*(alpha_max+xi_nu)+1)-1;
NQ = N-Q-1;
dataIdx = (Q+2):(NQ+1);
nData = numel(dataIdx);

U = daft(N, c1, c2)';

SNRd_dB = [10];
SNRp_list = [30 35 40];
numMonte = 20;

BER_ideal = zeros(1,numel(SNRd_dB));
BER_est = zeros(numel(SNRp_list), numel(SNRd_dB));

for isnr = 1:numel(SNRd_dB)
    snrd_lin = 10^(SNRd_dB(isnr)/10);
    N0 = 1/snrd_lin;
    errI = 0; totI = 0;
    errE = zeros(numel(SNRp_list),1); totE = zeros(numel(SNRp_list),1);
    for imc = 1:numMonte
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        delays = randperm(lmax+1, P) - 1;
        theta = -pi + 2*pi*rand(1,P);
        dopplers = alpha_max * cos(theta); % fractional, Jakes

        H = Gen_channel_mtx(N, P, gains, delays, dopplers/N, c1);
        Heff = U' * H * U;

        bits = randi([0 1], nData*2, 1);
        symData = qammod(bits, 4, 'gray', 'InputType','bit', 'UnitAveragePower', true);

        xVec = zeros(N,1);
        xVec(dataIdx) = symData;
        noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));

        % ---- Ideal CSI ----
        s = U*xVec;
        r = H*s + noise;
        z = U'*r;
        zData = z(dataIdx);
        HeffData = Heff(dataIdx, dataIdx);
        symHat = (HeffData'*HeffData + N0*eye(nData)) \ (HeffData' * zData);
        bitsHat = qamdemod(symHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
        errI = errI + sum(bitsHat ~= bits); totI = totI + numel(bits);

        % ---- Estimated CSI, for each SNRp ----
        for ip = 1:numel(SNRp_list)
            SNRp = SNRp_list(ip);
            xpilot = sqrt(N0 * 10^(SNRp/10));
            xVecP = xVec; xVecP(1) = xpilot;
            sP = U*xVecP;
            rP = H*sP + noise;
            zP = U'*rP;

            [gHat, lHat, nuHat] = afdm_pilot_estimate(zP, xpilot, N, P, lmax, alpha_max, xi_nu, c1, U', true);
            Hhat = Gen_channel_mtx(N, P, gHat, lHat, nuHat/N, c1);
            HeffHat = U' * Hhat * U;
            HeffHatData = HeffHat(dataIdx, dataIdx);
            zDataP = zP(dataIdx);
            symHatE = (HeffHatData'*HeffHatData + N0*eye(nData)) \ (HeffHatData' * zDataP);
            bitsHatE = qamdemod(symHatE, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
            errE(ip) = errE(ip) + sum(bitsHatE ~= bits); totE(ip) = totE(ip) + numel(bits);
        end
    end
    BER_ideal(isnr) = errI/totI;
    BER_est(:,isnr) = errE ./ totE;
    fprintf('SNRd=%2d dB: Ideal=%.4e  SNRp30=%.4e  SNRp35=%.4e  SNRp40=%.4e\n', ...
        SNRd_dB(isnr), BER_ideal(isnr), BER_est(1,isnr), BER_est(2,isnr), BER_est(3,isnr));
end

save('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12a_results.mat', 'BER_ideal','BER_est','SNRd_dB','SNRp_list','numMonte');

figure('Color','w');
semilogy(SNRd_dB, BER_ideal, '-ko', 'DisplayName','Ideal'); hold on;
mk={'-s','-d','-^'};
for ip=1:numel(SNRp_list)
    semilogy(SNRd_dB, BER_est(ip,:), mk{ip}, 'DisplayName', sprintf('SNRp = %d dB', SNRp_list(ip)));
end
grid on; xlabel('SNR in dB'); ylabel('BER');
title(sprintf('Fig 12a repro: AFDM channel estimation (fractional), N=%d, numMonte=%d', N, numMonte));
legend('show'); ylim([1e-4 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12a_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
