% Diagnostic: exact fig12a detection pipeline (pilot+data, MMSE on
% HeffHatData) at ONE fixed SNRd, tracking BER AND diagnostics for
% catastrophic failures (duplicate estimated (delay,doppler) pairs,
% condition number of HeffHatData) to see if failure rate grows with SNRp.
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

SNRd_dB = 9;
snrd_lin = 10^(SNRd_dB/10); N0 = 1/snrd_lin;
SNRp_list = [30 35 40];
numMonte = 600;

for ip = 1:numel(SNRp_list)
    SNRp = SNRp_list(ip);
    xpilot = sqrt(N0 * 10^(SNRp/10));
    errE = 0; totE = 0; nDup = 0; badCond = 0; condVals = zeros(numMonte,1);
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

        % detect duplicate (delay,doppler) estimates (rounded doppler to 0.1)
        key = lHat*100 + round(nuHat*10);
        if numel(unique(key)) < P
            nDup = nDup + 1;
        end

        Hhat = Gen_channel_mtx(N, P, gHat, lHat, nuHat/N, c1);
        HeffHat = U' * Hhat * U;
        HeffHatData = HeffHat(:, dataIdx);
        cnum = cond(HeffHatData'*HeffHatData + N0*eye(nData));
        condVals(imc) = cnum;
        if cnum > 1e6, badCond = badCond + 1; end

        symHatE = (HeffHatData'*HeffHatData + N0*eye(nData)) \ (HeffHatData' * zP);
        bitsHatE = qamdemod(symHatE, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
        errE = errE + sum(bitsHatE ~= bits); totE = totE + numel(bits);
    end
    fprintf('SNRp=%2d dB: BER=%.4e  dupRate=%.3f  badCondRate=%.3f  medianCond=%.3e  maxCond=%.3e\n', ...
        SNRp, errE/totE, nDup/numMonte, badCond/numMonte, median(condVals), max(condVals));
end

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
