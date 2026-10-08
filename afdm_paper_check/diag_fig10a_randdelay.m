% Diagnostic: does randomizing delays per trial (instead of fixed 0:lmax)
% produce OFDM's expected pronounced floor at high SNR (paper shows OFDM
% floor ~1.2e-4 at SNR=30dB, integer Doppler), while AFDM keeps dropping?
clear; clc; close all;
rng(31);

N = 256; P = 3; lmax = 2; alpha_max = 2;
c1_afdm = (2*alpha_max+1)/(2*N);
U_afdm = daft(N, c1_afdm, 1/(2*N^2))';
U_ofdm = dftmtx(N)'/sqrt(N);

SNR_dB_list = [20 25 30];
numMonte = 8000;

for isnr = 1:numel(SNR_dB_list)
    SNR_dB = SNR_dB_list(isnr);
    snr_lin = 10^(SNR_dB/10); N0 = 1/snr_lin;
    errAfdm = 0; errOfdm = 0; tot = 0;
    for imc = 1:numMonte
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        delays = randperm(lmax+1, min(P,lmax+1)); delays = delays-1;
        if numel(delays) < P
            delays = [delays, randi([0 lmax], 1, P-numel(delays))];
        end
        cand = -alpha_max:alpha_max;
        doppler_norm = cand(randperm(numel(cand), P));

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
            xHat = (Heff'*Heff + N0*eye(N)) \ (Heff'*z);
            bitsHat = qamdemod(xHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
            if w==1
                errAfdm = errAfdm + sum(bitsHat ~= bits);
            else
                errOfdm = errOfdm + sum(bitsHat ~= bits);
            end
        end
        tot = tot + numel(bits);
    end
    fprintf('SNR=%2d dB: AFDM BER=%.4e  OFDM BER=%.4e  (ratio=%.1fx)\n', SNR_dB, errAfdm/tot, errOfdm/tot, (errOfdm/tot)/(errAfdm/tot));
end

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
