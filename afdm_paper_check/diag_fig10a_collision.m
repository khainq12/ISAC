clear; clc;
N = 256; P = 3; lmax = 2; alpha_max = 2;
c1_afdm = (2*alpha_max+1)/(2*N);
U_afdm = daft(N, c1_afdm, 1/(2*N^2))';
U_ofdm = dftmtx(N)'/sqrt(N);
SNRd = 15; N0 = 1/10^(SNRd/10);
numMonte = 500;
rng(3);

errAFDM=0; errOFDM=0; tot=0;
for imc=1:numMonte
    gains = sqrt(1/P) * (randn(1,P) + 1j*randn(1,P))/sqrt(2);
    delays = 0:lmax;
    cand = -alpha_max:alpha_max;
    doppler_norm = cand(randperm(numel(cand), P)); % DISTINCT now
    bits = randi([0 1], N*2, 1);
    sym = qammod(bits, 4, 'gray', 'InputType','bit', 'UnitAveragePower', true);
    noise = sqrt(N0/2) * (randn(N,1) + 1j*randn(N,1));
    for w=1:2
        if w==1, Uw=U_afdm; c1w=c1_afdm; else, Uw=U_ofdm; c1w=0; end
        Hw = Gen_channel_mtx(N, P, gains, delays, doppler_norm/N, c1w);
        s = Uw*sym; r = Hw*s+noise; Heff = Uw'*Hw*Uw; z = Uw'*r;
        xHat = (Heff'*Heff + N0*eye(N)) \ (Heff'*z);
        bitsHat = qamdemod(xHat, 4, 'gray', 'OutputType','bit', 'UnitAveragePower', true);
        if w==1, errAFDM = errAFDM + sum(bitsHat~=bits); else, errOFDM = errOFDM + sum(bitsHat~=bits); end
    end
    tot = tot + numel(bits);
end
fprintf('DISTINCT Doppler: AFDM BER=%.4e  OFDM BER=%.4e (should be close now)\n', errAFDM/tot, errOFDM/tot);

function A = daft(N, c1, c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
