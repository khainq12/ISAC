clear; clc;
N = 16;
alpha_max = 1;
c1 = (2*alpha_max+1)/(2*N);
c2 = 1/(N^2);
P = 3; lmax = 2;

nHyp = 2^N;
S_hyp = zeros(N, nHyp);
for idx = 0:nHyp-1
    bits = bitget(idx, N:-1:1);
    S_hyp(:,idx+1) = 1 - 2*bits;
end
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
X_hyp = A' * S_hyp;

rng(1);
numTest = 200;
snr_lin = 10^(10/10); N0 = 1/snr_lin;
tic;
for imc = 1:numTest
    gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
    delays = randperm(lmax+1, P) - 1;
    dopplers_int = randi([-alpha_max alpha_max], 1, P);
    H = Gen_channel_mtx(N, P, gains, delays, dopplers_int/N, c1);
    xBits = randi([0 1], N, 1);
    xSym = 1 - 2*xBits;
    s = A' * xSym;
    noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));
    r = H*s + noise;
    HX = H * X_hyp;
    dist2 = sum(abs(r - HX).^2, 1);
    [~, kBest] = min(dist2);
end
elapsed = toc;
fprintf('=== %d trials took %.3f sec => %.5f sec/trial ===\n', numTest, elapsed, elapsed/numTest);
