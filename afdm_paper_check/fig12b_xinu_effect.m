% Reproduction of Fig. 12b in Bemani, Ksairi, Kountouris (arXiv:2204.12798).
% AFDM, N=256, 4QAM, effect of xi_nu (1,4,8) on BER vs SNRd, comparing
% Ideal CSI vs estimated CSI (SNRp = 40 dB), fractional Doppler.
clear; clc; close all;
rng(7);

N = 256; P = 3; lmax = 2; alpha_max = 2;
SNRd_dB = 0:3:18;
SNRp = 40;
numMonte = 300;
xiList = [1 4 8];

BER_ideal = zeros(numel(xiList), numel(SNRd_dB));
BER_est = zeros(numel(xiList), numel(SNRd_dB));

for ix = 1:numel(xiList)
    xi_nu = xiList(ix);
    c1 = (2*(alpha_max+xi_nu)+1)/(2*N);
    c2 = 1/(2*N^2);
    Q = (lmax+1)*(2*(alpha_max+xi_nu)+1)-1;
    NQ = N-Q-1;
    dataIdx = (Q+2):(NQ+1);
    nData = numel(dataIdx);
    U = daft(N, c1, c2)';

    for isnr = 1:numel(SNRd_dB)
        snrd_lin = 10^(SNRd_dB(isnr)/10);
        N0 = 1/snrd_lin;
        errI = 0; totI = 0; errE = 0; totE = 0;
        for imc = 1:numMonte
            gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
            delays = randperm(lmax+1, P) - 1;
            theta = -pi + 2*pi*rand(1,P);
            dopplers = alpha_max * cos(theta);

            H = Gen_channel_mtx(N, P, gains, delays, dopplers/N, c1);
            Heff = U' * H * U;

            bits = randi([0 1], nData*2, 1);
            symData = qammod(bits, 4, 'gray', 'InputType','bit', 'UnitAveragePower', true);
            xVec = zeros(N,1); xVec(dataIdx) = symData;
            noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));

            % Ideal - use the FULL received vector z (all N observations)
            % and the rectangular column-selected Heff(:,dataIdx), since
            % z = Heff(:,dataIdx)*x_data + noise EXACTLY (guard/pilot
            % positions transmit 0). Restricting to the square submatrix
            % Heff(dataIdx,dataIdx) as before discards received energy
            % that leaked (via AFDM's non-diagonal Heff) into the guard
            % region, which grows with xi_nu - this was the bug causing
            % the xi_nu trend to invert relative to the paper.
            s = U*xVec; r = H*s + noise; z = U'*r;
            HeffData = Heff(:, dataIdx);
            symHat = (HeffData'*HeffData + N0*eye(nData)) \ (HeffData' * z);
            bitsHat = qamdemod(symHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
            errI = errI + sum(bitsHat ~= bits); totI = totI + numel(bits);

            % Estimated, SNRp=40
            xpilot = sqrt(N0 * 10^(SNRp/10));
            xVecP = xVec; xVecP(1) = xpilot;
            sP = U*xVecP; rP = H*sP + noise; zP = U'*rP;
            [gHat, lHat, nuHat] = afdm_pilot_estimate(zP, xpilot, N, P, lmax, alpha_max, xi_nu, c1, U', true);
            Hhat = Gen_channel_mtx(N, P, gHat, lHat, nuHat/N, c1);
            HeffHat = U' * Hhat * U;
            HeffHatData = HeffHat(:, dataIdx);
            % Subtract known pilot leakage before data detection (see
            % fig11_channel_estimation.m note - uncancelled
            % HeffHat(:,1)*xpilot otherwise biases all Est curves high).
            zData = zP - HeffHat(:,1) * xpilot;
            symHatE = (HeffHatData'*HeffHatData + N0*eye(nData)) \ (HeffHatData' * zData);
            bitsHatE = qamdemod(symHatE, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
            errE = errE + sum(bitsHatE ~= bits); totE = totE + numel(bits);
        end
        BER_ideal(ix,isnr) = errI/totI;
        BER_est(ix,isnr) = errE/totE;
        fprintf('xi_nu=%d SNRd=%2d dB: Ideal=%.4e  Est(SNRp=40)=%.4e\n', xi_nu, SNRd_dB(isnr), BER_ideal(ix,isnr), BER_est(ix,isnr));
    end
end

save('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12b_results.mat', 'BER_ideal','BER_est','SNRd_dB','xiList','numMonte');

figure('Color','w');
mk = {'-o','-s','-^'};
clr = {[0 0.45 0.74],[0.85 0.33 0.1],[0.47 0.67 0.19]};
for ix=1:numel(xiList)
    semilogy(SNRd_dB, BER_est(ix,:), mk{ix}, 'Color', clr{ix}, 'DisplayName', sprintf('Est, \\xi_\\nu=%d', xiList(ix))); hold on;
    semilogy(SNRd_dB, BER_ideal(ix,:), mk{ix}, 'LineStyle', '--', 'Color', clr{ix}, 'DisplayName', sprintf('Ideal, \\xi_\\nu=%d', xiList(ix)));
end
grid on; xlabel('SNR in dB'); ylabel('BER');
title(sprintf('Fig 12b repro: effect of \\xi_\\nu, N=%d, SNRp=40dB, numMonte=%d', N, numMonte));
legend('show','FontSize',8); ylim([1e-4 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12b_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
