% Muc 14 velocity sweep - v6, FRESH CODE. v5 (multipath + smooth estimation-error
% model) restored waveform differentiation (OFDM ~2x worse than OTFS/AFDM/RM) but
% reintroduced a real, non-monotonic ripple in the LOW-speed region (v<=80km/h,
% where estimation error ~=0): because even under EXACT channel knowledge, one
% FIXED background realization's phase relationship with the target's Doppler
% creates a reproducible per-realization collision/ripple pattern (same root cause
% established in v1/v2's genie-CSI runs).
%
% v6 fixes this the same way this report already handles "single fixed channel is
% not a general trend" elsewhere (Muc 12's diversity-order methodology, which
% explicitly requires ensemble averaging over Rayleigh realizations rather than one
% fixed channel): instead of ONE fixed background channel reused every trial, each
% Monte Carlo trial draws a FRESH Rayleigh-faded background (same RMS gain
% magnitudes / delays / Doppler support as the fixed background used in Muc
% 9.7/11/13, but independent complex-Gaussian phase+magnitude per trial). Averaging
% over this ensemble removes the one-realization coincidental phase alignment while
% preserving the genuine average waveform-structural difference (OFDM's lack of a
% Doppler-domain path-separation mechanism, proven in Part I Muc 2.1/4) and the
% growing target-Doppler estimation-error mechanism from v4/v5.
clc; clear;
L = 256; M = 16; N = 16;
c = 3e8; fc = 28e9; fs = 500e3;

bgDelays  = [3, 11, 17];
bgRMS     = [0.75, 0.55, 0.38]; % same magnitudes as the fixed background in Muc 9.7/11/13
bgDopplers = [1, -2, 3];        % same Doppler support as Muc 9.7/11/13
targetDelay = 0;
targetGain  = 1.00*exp(1j*0.0);

velocities = [20 40 60 80 100 120 140 160 180 200 220 250]; % km/h
dopplerOfL = @(v_kmh) L * (2*(v_kmh/3.6)*fc/c) / fs;
designMaxNu = 1.5;
estSlope = 0.02;

snrList = [10];
numMonte = 4000;

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

n = (0:L-1).';
I = eye(L);
S_tau = cell(1,numel(bgDelays)); D_nu_bg = cell(1,numel(bgDelays));
for p = 1:numel(bgDelays)
    S_tau{p} = circshift(I, bgDelays(p), 1);
    D_nu_bg{p} = diag(exp(1j*2*pi*bgDopplers(p)*n/L));
end

fprintf('--- Target Doppler nu / estimation-error std at each tested velocity ---\n');
for vi = 1:numel(velocities)
    nuT = dopplerOfL(velocities(vi));
    sig = estSlope * max(0, nuT - designMaxNu);
    fprintf('v=%dkm/h -> nu_target=%.4f, std(nu_hat)=%.4f\n', velocities(vi), nuT, sig);
end

BER = zeros(numel(wfNames), numel(velocities), numel(snrList));

for si = 1:numel(snrList)
    snrDb = snrList(si);
    noiseVar = 10^(-snrDb/10);
    for vi = 1:numel(velocities)
        nuTrue = dopplerOfL(velocities(vi));
        sigNu = estSlope * max(0, nuTrue - designMaxNu);
        bitErr = zeros(numel(wfNames),1); bitTot = 0;
        rng(6000*si); % SAME seed at every velocity (common random numbers): pairs the
        % channel/data/noise draws across velocities so only nuTrue differs trial-by-
        % trial -- removes velocity-uncorrelated noise, leaving the smooth true trend.
        for t = 1:numMonte
            % fresh Rayleigh-faded background realization this trial
            Hbg = zeros(L,L);
            for p = 1:numel(bgDelays)
                g = bgRMS(p)/sqrt(2)*(randn()+1j*randn());
                Hbg = Hbg + g*D_nu_bg{p}*S_tau{p};
            end

            nuHat = nuTrue + sigNu*randn();
            Dtrue = exp(1j*2*pi*nuTrue*n/L);
            Dhat  = exp(1j*2*pi*nuHat*n/L);
            Htrue = Hbg + targetGain*diag(Dtrue);
            Hhat  = Hbg + targetGain*diag(Dhat);

            bitTx = randi([0 1], 2*L, 1);
            symTx = qpsk_mod_local(bitTx);
            dataNoise = sqrt(noiseVar/2)*(randn(L,1)+1j*randn(L,1));
            for wi = 1:numel(wfNames)
                U = ops.(wfNames{wi}); W = U';
                x = U*symTx;
                y = Htrue*x + dataNoise;
                z = W*y;
                Heff_est = W*Hhat*U;
                A = Heff_est'*Heff_est + noiseVar*eye(L);
                b = Heff_est'*z;
                symHat = A\b;
                bitHat = qpsk_hard_demod_local(symHat);
                bitErr(wi) = bitErr(wi) + sum(bitHat~=bitTx);
            end
            bitTot = bitTot + numel(bitTx);
        end
        BER(:,vi,si) = bitErr/bitTot;
    end
    fprintf('SNR=%ddB done\n', snrDb);
end

save('velocity_v6_ensemble_results.mat','velocities','snrList','BER','wfNames','designMaxNu','estSlope');
for si=1:numel(snrList)
    fprintf('--- SNR=%d dB ---\n', snrList(si));
    for wi=1:numel(wfNames)
        fprintf('%s: %s\n', wfNames{wi}, mat2str(squeeze(BER(wi,:,si)),4));
    end
end
fprintf('=== DONE ===\n');

function s = qpsk_mod_local(b)
b = b(:); I = 1-2*b(1:2:end); Q = 1-2*b(2:2:end); s = (I+1j*Q)/sqrt(2);
end
function b = qpsk_hard_demod_local(s)
I = real(s)<0; Q = imag(s)<0;
b = zeros(2*numel(s),1); b(1:2:end)=I; b(2:2:end)=Q;
end
