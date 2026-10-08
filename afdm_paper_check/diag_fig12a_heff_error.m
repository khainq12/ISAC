% Diagnostic: measure ||HeffHat - HeffTrue||_F / ||HeffTrue||_F (the
% quantity that actually drives MMSE detection quality) across SNRp, to
% see if it behaves differently from the per-path parameter MAE (which
% improved monotonically) - e.g. due to a gain SIGN/PHASE issue that
% parameter MAE (measuring |gain error|) wouldn't catch but which could
% still corrupt the reconstructed Heff.
clear; clc; close all;
rng(11);

N = 256; P = 3; lmax = 2; alpha_max = 2; xi_nu = 1;
c1 = (2*(alpha_max+xi_nu)+1)/(2*N);
c2 = 1/(2*N^2);
Q = (lmax+1)*(2*(alpha_max+xi_nu)+1)-1;
NQ = N-Q-1;
dataIdx = (Q+2):(NQ+1);
U = daft(N, c1, c2)';

SNRd_dB = 9; snrd_lin = 10^(SNRd_dB/10); N0 = 1/snrd_lin;
SNRp_list = [30 35 40];
numMonte = 300;

for ip = 1:numel(SNRp_list)
    SNRp = SNRp_list(ip);
    xpilot = sqrt(N0 * 10^(SNRp/10));
    relErrAcc = 0; relErrDataAcc = 0;
    for imc = 1:numMonte
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        delays = randperm(lmax+1, P) - 1;
        theta = -pi + 2*pi*rand(1,P);
        dopplers = alpha_max * cos(theta);
        H = Gen_channel_mtx(N, P, gains, delays, dopplers/N, c1);
        HeffTrue = U' * H * U;

        bits = randi([0 1], numel(dataIdx)*2, 1);
        symData = qammod(bits, 4, 'gray', 'InputType','bit', 'UnitAveragePower', true);
        xVecP = zeros(N,1); xVecP(dataIdx) = symData; xVecP(1) = xpilot;
        sP = U*xVecP;
        noise = sqrt(N0/2)*(randn(N,1)+1j*randn(N,1));
        rP = H*sP + noise;
        zP = U'*rP;

        [gHat, lHat, nuHat] = afdm_pilot_estimate(zP, xpilot, N, P, lmax, alpha_max, xi_nu, c1, U', true);
        Hhat = Gen_channel_mtx(N, P, gHat, lHat, nuHat/N, c1);
        HeffHat = U' * Hhat * U;

        relErrAcc = relErrAcc + norm(HeffHat - HeffTrue, 'fro') / norm(HeffTrue, 'fro');
        relErrDataAcc = relErrDataAcc + norm(HeffHat(:,dataIdx) - HeffTrue(:,dataIdx), 'fro') / norm(HeffTrue(:,dataIdx), 'fro');
    end
    fprintf('SNRp=%2d dB: relFrobErr(full)=%.4f  relFrobErr(dataCols)=%.4f\n', SNRp, relErrAcc/numMonte, relErrDataAcc/numMonte);
end

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
