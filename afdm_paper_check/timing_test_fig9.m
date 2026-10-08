clear; clc;
N = 256; NOTFS=16; MOTFS=16; P=3; lmax=2; alpha_max=2; xi_nu=1;
c1_afdm = (2*(alpha_max+xi_nu)+1)/(2*N);
c2_afdm = 1/(2*N^2);
c1_ocdm = 1/(2*N); c2_ocdm = 1/(2*N);
F = dftmtx(N)/sqrt(N);
A_ofdm = F';
A_ocdm = daft(N,c1_ocdm,c2_ocdm);
A_afdm = daft(N,c1_afdm,c2_afdm);
FN = dftmtx(NOTFS)/sqrt(NOTFS);
A_otfs = kron(FN', eye(MOTFS));
A_list = {A_ofdm, A_ocdm, A_otfs, A_afdm};
rng(1);
N0 = 1/10^(1);
numTest = 200;
tic;
for imc=1:numTest
    gains = sqrt(1/P)*(randn(1,P)+1j*randn(1,P))/sqrt(2);
    delays = randperm(lmax+1,P)-1;
    theta = -pi+2*pi*rand(1,P);
    doppler_norm = alpha_max*cos(theta);
    bits = randi([0 1], N*2,1);
    sym = qammod(bits,4,'gray','InputType','bit','UnitAveragePower',true);
    noise = sqrt(N0/2)*(randn(N,1)+1j*randn(N,1));
    for w=1:4
        switch w
            case 1, c1w=0; case 2, c1w=c1_ocdm; case 3, c1w=0; case 4, c1w=c1_afdm;
        end
        Hw = Gen_channel_mtx(N,P,gains,delays,doppler_norm/N,c1w);
        Aw = A_list{w};
        s = Aw*sym;
        r = Hw*s+noise;
        Heff = Aw'*Hw*Aw;
        symHat = (Heff'*Heff+N0*eye(N))\(Heff'*(Aw'*r));
        bitsHat = qamdemod(symHat,4,'gray','OutputType','bit','UnitAveragePower',true);
    end
end
elapsed=toc;
fprintf('=== %d trials (4 waveforms each) took %.3f sec => %.5f sec/trial ===\n', numTest, elapsed, elapsed/numTest);

function A = daft(N,c1,c2)
F = dftmtx(N)/sqrt(N);
L1 = diag(exp(-1i*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1i*2*pi*c2*((0:N-1).^2)));
A = L2*F*L1;
end
