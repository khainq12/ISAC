clear; clc; close all;
rng(5);
N = 256; P = 3; lmax = 2; alpha_max = 2; xi_nu = 1;
c1_afdm = (2*(alpha_max+xi_nu)+1)/(2*N);
U_afdm = daft(N, c1_afdm, 1/(2*N^2))';
U_ofdm = dftmtx(N)'/sqrt(N);

SNR_dB = 20; snr_lin = 10^(SNR_dB/10); N0 = 1/snr_lin;
numMonte = 400;
Lvals = [9 15 21 27 33 45];

errAfdm = zeros(numel(Lvals),1); errOfdm = zeros(numel(Lvals),1);
tot = 0;
for imc = 1:numMonte
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
        for iL = 1:numel(Lvals)
            HeffB = topRowMagnitudeMask(Heff, Lvals(iL));
            xHat = (HeffB'*HeffB + N0*eye(N)) \ (HeffB' * z);
            bitsHat = qamdemod(xHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
            if w==1
                errAfdm(iL) = errAfdm(iL) + sum(bitsHat ~= bits);
            else
                errOfdm(iL) = errOfdm(iL) + sum(bitsHat ~= bits);
            end
        end
    end
    tot = tot + numel(bits);
end
for iL=1:numel(Lvals)
    fprintf('L=%3d: AFDM-lowC BER=%.4e  OFDM-lowC BER=%.4e\n', Lvals(iL), errAfdm(iL)/tot, errOfdm(iL)/tot);
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
