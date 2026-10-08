% Diagnostic: replicate fig12a's EXACT per-trial setup (pilot AT position 1
% PLUS data symbols at dataIdx, both present when estimating), to see if
% data-symbol interference is what breaks monotonicity in SNRp for the
% pilot-based delay/Doppler/gain estimator (the pilot-only diagnostic
% diag_fig12a_snrp.m showed clean monotonic improvement with no data).
clear; clc; close all;
rng(11);

N = 256; P = 3; lmax = 2; alpha_max = 2; xi_nu = 1;
c1 = (2*(alpha_max+xi_nu)+1)/(2*N);
c2 = 1/(2*N^2);
Q = (lmax+1)*(2*(alpha_max+xi_nu)+1)-1;
NQ = N-Q-1;
dataIdx = (Q+2):(NQ+1);
nData = numel(dataIdx);
U = daft(N, c1, c2)';

SNRd_dB = 9; % fixed, mid-range, matches where reversal was clear in fig12a
snrd_lin = 10^(SNRd_dB/10); N0 = 1/snrd_lin;
SNRp_list = [30 35 40];
numMonte = 300;

gErr = zeros(numel(SNRp_list),1); lErr = zeros(numel(SNRp_list),1); nuErr = zeros(numel(SNRp_list),1);

for ip = 1:numel(SNRp_list)
    SNRp = SNRp_list(ip);
    xpilot = sqrt(N0 * 10^(SNRp/10));
    gAcc=0; lAcc=0; nuAcc=0;
    for imc = 1:numMonte
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        delays = randperm(lmax+1, P) - 1;
        theta = -pi + 2*pi*rand(1,P);
        dopplers = alpha_max * cos(theta);

        H = Gen_channel_mtx(N, P, gains, delays, dopplers/N, c1);

        bits = randi([0 1], nData*2, 1);
        symData = qammod(bits, 4, 'gray', 'InputType','bit', 'UnitAveragePower', true);
        xVecP = zeros(N,1); xVecP(dataIdx) = symData; xVecP(1) = xpilot;

        sP = U*xVecP;
        noise = sqrt(N0/2)*(randn(N,1)+1j*randn(N,1));
        rP = H*sP + noise;
        zP = U'*rP;

        [gHat, lHat, nuHat] = afdm_pilot_estimate(zP, xpilot, N, P, lmax, alpha_max, xi_nu, c1, U', true);

        used = false(1,P);
        for k=1:P
            best=-1; bd=inf;
            for t=1:P
                if used(t), continue; end
                d = abs(delays(t)-lHat(k)) + abs(dopplers(t)-nuHat(k));
                if d<bd, bd=d; best=t; end
            end
            used(best)=true;
            gAcc = gAcc + abs(gains(best)-gHat(k))^2;
            lAcc = lAcc + abs(delays(best)-lHat(k));
            nuAcc = nuAcc + abs(dopplers(best)-nuHat(k));
        end
    end
    gErr(ip) = gAcc/(numMonte*P); lErr(ip) = lAcc/(numMonte*P); nuErr(ip) = nuAcc/(numMonte*P);
    fprintf('SNRp=%2d dB (with data interference): gainMSE=%.4e  delayMAE=%.4f  dopplerMAE=%.4f\n', SNRp, gErr(ip), lErr(ip), nuErr(ip));
end

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
