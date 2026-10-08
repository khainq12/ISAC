clear; clc;
N = 256; P = 3; lmax = 2; alpha_max = 2; xi_nu = 1;
c1 = (2*(alpha_max+xi_nu)+1)/(2*N);
c2 = 1/(2*N^2);
U = daft(N, c1, c2)';

rng(11);
gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
delays = [0 1 2];
trueA = [0.3, -0.2, 0.45];
trueAlphaInt = [1, -2, 0];
trueNu = trueAlphaInt + trueA;

H = Gen_channel_mtx(N, P, gains, delays, trueNu/N, c1);

xpilot = 10; % very strong pilot, noiseless test
xVec = zeros(N,1); xVec(1) = xpilot;
s = U*xVec;
r = H*s; % NO NOISE
z = U'*r;

[gHat, lHat, nuHat] = afdm_pilot_estimate(z, xpilot, N, P, lmax, alpha_max, xi_nu, c1, U', true);

fprintf('True delays  : %s\n', num2str(sort(delays)));
fprintf('Est  delays  : %s\n', num2str(sort(lHat)));
fprintf('True nu      : %s\n', num2str(sort(trueNu)));
fprintf('Est  nu      : %s\n', num2str(sort(nuHat)));
fprintf('True gains   : %s\n', num2str(gains));
fprintf('Est  gains   : %s\n', num2str(gHat));

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
