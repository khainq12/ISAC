% Diagnostic: does afdm_pilot_estimate's accuracy degrade at high SNRp
% due to matching-pursuit cancellation leakage under a coarse aFine grid?
clear; clc; close all;
rng(42);

N = 256; P = 3; lmax = 2; alpha_max = 2; xi_nu = 1;
c1 = (2*(alpha_max+xi_nu)+1)/(2*N);
c2 = 1/(2*N^2);
U = dftmtx_daft(N, c1, c2)';

SNRp_list = [20 25 30 35 40 50 70];
numMonte = 200;
N0 = 1; % fix data-noise scale; only pilot amplitude varies with SNRp

gErr = zeros(numel(SNRp_list),1);
lErr = zeros(numel(SNRp_list),1);
nuErr = zeros(numel(SNRp_list),1);

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
        xVecP = zeros(N,1); xVecP(1) = xpilot;
        sP = U*xVecP;
        noise = sqrt(N0/2)*(randn(N,1)+1j*randn(N,1));
        rP = H*sP + noise;
        zP = U'*rP;

        [gHat, lHat, nuHat] = afdm_pilot_estimate(zP, xpilot, N, P, lmax, alpha_max, xi_nu, c1, U', true);

        % match estimated paths to true paths by nearest delay+doppler
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
    gErr(ip) = gAcc/(numMonte*P);
    lErr(ip) = lAcc/(numMonte*P);
    nuErr(ip) = nuAcc/(numMonte*P);
    fprintf('SNRp=%2d dB: gainMSE=%.4e  delayMAE=%.4f  dopplerMAE=%.4f\n', SNRp, gErr(ip), lErr(ip), nuErr(ip));
end

function A = dftmtx_daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
