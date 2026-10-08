% Fix: at high SNR (>=20dB, integer case), BER counts are so low (~4-8
% errors out of numMonte=1500) that curves show non-monotonic "upticks"
% purely from Monte Carlo noise (verified: e.g. AFDM MMSE went from
% 5.2e-6 at SNR=27.5 to 1.04e-5 at SNR=30, an increase, driven by only
% ~4 vs ~8 counted error events - well within Poisson noise). Recompute
% ONLY the high-SNR tail points with much higher numMonte to smooth this,
% keeping low/mid-SNR points (already well-resolved) untouched.
clear; clc; close all;
rng(21);

load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig10_results.mat');

N = 256; P = 3; lmax = 2; alpha_max = 2;
c1_afdm = (2*alpha_max+1)/(2*N);
U_afdm = daft(N, c1_afdm, 1/(2*N^2))';
U_ofdm = dftmtx(N)'/sqrt(N);

SNR_dB = results.integer.SNR_dB;
tailIdx = find(SNR_dB >= 20); % 20,22.5,25,27.5,30
numMonteTail = 30000;

for ii = 1:numel(tailIdx)
    isnr = tailIdx(ii);
    snr_lin = 10^(SNR_dB(isnr)/10);
    N0 = 1/snr_lin;
    errB = zeros(3,2); totB = zeros(3,2);
    for imc = 1:numMonteTail
        gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
        delays = 0:lmax;
        theta = -pi + 2*pi*rand(1,P);
        doppler_norm = round(alpha_max * cos(theta)); % Jakes, rounded (matches fig10_detectors.m)

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

            xHat1 = (Heff'*Heff + N0*eye(N)) \ (Heff'*z);
            Rmat = Heff'*Heff + N0*eye(N);
            bvec = Heff'*z;
            Dg = diag(diag(Rmat)); Lg = tril(Rmat,-1); Ug = triu(Rmat,1);
            xIter = zeros(N,1);
            for it=1:15
                xIter = (Dg+Lg) \ (bvec - Ug*xIter);
            end
            xHat2 = xIter;
            HeffB = topRowMagnitudeMask(Heff, P); % integer case: exact, L=P
            xHat3 = (HeffB'*HeffB + N0*eye(N)) \ (HeffB'*z);

            for d = 1:3
                switch d
                    case 1, xHat = xHat1; case 2, xHat = xHat2; case 3, xHat = xHat3;
                end
                bitsHat = qamdemod(xHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
                errB(d,w) = errB(d,w) + sum(bitsHat ~= bits);
                totB(d,w) = totB(d,w) + numel(bits);
            end
        end
    end
    newBER = errB ./ totB;
    results.integer.BER(:,:,isnr) = newBER;
    fprintf('SNR=%5.1f (numMonte=%d): AFDM(MMSE/DFE/lowC)=%.3e/%.3e/%.3e  OFDM(MMSE/DFE/lowC)=%.3e/%.3e/%.3e\n', ...
        SNR_dB(isnr), numMonteTail, newBER(1,1),newBER(2,1),newBER(3,1),newBER(1,2),newBER(2,2),newBER(3,2));
end

save('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig10_results.mat', 'results','numMonte');
fprintf('=== DONE ===\n');

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
