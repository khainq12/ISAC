clear; clc;
N = 256; P = 3; lmax = 2; alpha_max = 2;
SNRd_dB = 10;
numMonte = 300;
xiList = [1 8];

for ix = 1:numel(xiList)
    xi_nu = xiList(ix);
    c1 = (2*(alpha_max+xi_nu)+1)/(2*N);
    c2 = 1/(2*N^2);
    Q = (lmax+1)*(2*(alpha_max+xi_nu)+1)-1;
    NQ = N-Q-1;
    dataIdx = (Q+2):(NQ+1);
    nData = numel(dataIdx);
    U = daft(N, c1, c2)';

    snrd_lin = 10^(SNRd_dB/10);
    N0 = 1/snrd_lin;
    rng(7);
    errI = 0; totI = 0;
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
        s = U*xVec; r = H*s + noise; z = U'*r;
        HeffData = Heff(:, dataIdx);
        symHat = (HeffData'*HeffData + N0*eye(nData)) \ (HeffData' * z);
        bitsHat = qamdemod(symHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
        errI = errI + sum(bitsHat ~= bits); totI = totI + numel(bits);
    end
    fprintf('xi_nu=%d (Q=%d, nData=%d): Ideal BER = %.4e\n', xi_nu, Q, nData, errI/totI);
end
fprintf('=== DONE ===\n');

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
