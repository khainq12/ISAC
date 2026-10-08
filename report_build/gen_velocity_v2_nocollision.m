% Muc 14 velocity sweep - v2, FRESH RUN, SAME LMMSE/fixed-channel methodology as
% Muc 9.7/11/13.1-13.2/13.3 (get_waveform_operator-style unitary operators, LMMSE
% equalizer, one swept "target" path + fixed background multipath), but with the
% background Dopplers deliberately placed OUTSIDE the swept target-Doppler range so
% the target path's Doppler can never numerically collide with a background path's
% Doppler. This isolates whether BER-vs-velocity is monotonic once the (already-proven-
% real) Doppler-collision resonance mechanism is removed by construction, using the
% same waveform/channel/equalizer building blocks as every other new-methodology section.
clc; clear;
L = 256; M = 16; N = 16;
c = 3e8; fc = 28e9; fs = 500e3;

% background multipath: delays same as Muc 9.7, but Dopplers moved OUTSIDE the
% swept target Doppler range (target sweeps roughly nu in [0.5, 6.7] for v in
% [20,250] km/h at this fc/fs -- see printed nu values below). Background nu chosen
% far outside that band so no collision is possible at any tested velocity.
bgDelays   = [3, 11, 17];
bgDopplers = [10, -9, 13];
bgGains    = [0.75*exp(1j*0.4), 0.55*exp(-1j*0.9), 0.38*exp(1j*1.2)];
targetDelay = 0;
targetGain  = 1.00*exp(1j*0.0);

velocities = [20 40 60 80 100 120 140 160 180 200 220 250]; % km/h
dopplerOfL = @(v_kmh) L * (2*(v_kmh/3.6)*fc/c) / fs;

snrList = [10];
numMonte = 8000; % raised from 1500 -> cut Monte Carlo noise (~2.3x) since true dynamic range without collision is small

rmSeed = 314159;
rng(rmSeed);
permRM = randperm(L);
PermMat = eye(L); PermMat = PermMat(permRM,:);
F = exp(-1j*2*pi*(0:L-1)'*(0:L-1)/L)/sqrt(L);

c1afdm = 5/L; c2afdm = 1/(L^2);
nvec = (0:L-1).';
Lambda1 = diag(exp(1j*pi*c1afdm*nvec.^2));
Lambda2 = diag(exp(1j*pi*c2afdm*nvec.^2));

FN = zeros(N); for n=0:N-1; for k=0:N-1; FN(n+1,k+1)=exp(-1j*2*pi*n*k/N)/sqrt(N); end; end

ops = struct();
ops.OFDM = F';
ops.OTFS = kron(FN', eye(M));
ops.AFDM = Lambda1 * F' * Lambda2;
ops.RM   = PermMat * F;
wfNames = {'OFDM','OTFS','AFDM','RM'};

fprintf('--- Target Doppler nu at each tested velocity (background fixed at 10,-9,13; no overlap expected) ---\n');
for vi = 1:numel(velocities)
    fprintf('v=%dkm/h -> nu_target=%.4f\n', velocities(vi), dopplerOfL(velocities(vi)));
end

BER = zeros(numel(wfNames), numel(velocities), numel(snrList));

for si = 1:numel(snrList)
    snrDb = snrList(si);
    noiseVar = 10^(-snrDb/10);
    for vi = 1:numel(velocities)
        targetDoppler = dopplerOfL(velocities(vi));
        delaysAll = [targetDelay, bgDelays];
        dopplersAll = [targetDoppler, bgDopplers];
        gainsAll = [targetGain, bgGains];
        H_comm = build_delay_doppler_channel_local(L, delaysAll, dopplersAll, gainsAll);
        for wi = 1:numel(wfNames)
            U = ops.(wfNames{wi}); W = U';
            Heff = W * H_comm * U;
            A = Heff'*Heff + noiseVar*eye(L);
            bitErr = 0; bitTot = 0;
            rng(2000*si+vi);
            for t = 1:numMonte
                bitTx = randi([0 1], 2*L, 1);
                symTx = qpsk_mod_local(bitTx);
                x = U*symTx;
                noise = sqrt(noiseVar/2)*(randn(L,1)+1j*randn(L,1));
                y = H_comm*x + noise;
                z = W*y;
                b = Heff'*z;
                symHat = A\b;
                bitHat = qpsk_hard_demod_local(symHat);
                bitErr = bitErr + sum(bitHat~=bitTx);
                bitTot = bitTot + numel(bitTx);
            end
            BER(wi,vi,si) = bitErr/bitTot;
        end
    end
    fprintf('SNR=%ddB done\n', snrDb);
end

save('velocity_v2_nocollision_results.mat','velocities','snrList','BER','wfNames','bgDopplers');
for si=1:numel(snrList)
    fprintf('--- SNR=%d dB ---\n', snrList(si));
    for wi=1:numel(wfNames)
        fprintf('%s: %s\n', wfNames{wi}, mat2str(squeeze(BER(wi,:,si)),4));
    end
end
fprintf('=== DONE ===\n');

function H = build_delay_doppler_channel_local(L, delays, dopplers, gains)
n = (0:L-1).'; I = eye(L); H = zeros(L,L);
for p=1:numel(gains)
    S_tau = circshift(I, delays(p), 1);
    D_nu = diag(exp(1j*2*pi*dopplers(p)*n/L));
    H = H + gains(p)*D_nu*S_tau;
end
end
function s = qpsk_mod_local(b)
b = b(:); I = 1-2*b(1:2:end); Q = 1-2*b(2:2:end); s = (I+1j*Q)/sqrt(2);
end
function b = qpsk_hard_demod_local(s)
I = real(s)<0; Q = imag(s)<0;
b = zeros(2*numel(s),1); b(1:2:end)=I; b(2:2:end)=Q;
end
