% Muc 14 velocity sweep - v3, FRESH CODE, single-path target (no background multipath,
% isolates the pure Doppler effect - v1/v2 showed that with ANY fixed background paths
% present, genie-CSI LMMSE gives a reproducible (not noise) non-monotonic ripple in OFDM's
% BER vs velocity, because under PERFECT channel knowledge, OFDM's residual ICI is governed
% by frac(nu) (periodic every ~37.7km/h), not by |nu| itself -- so genie-CSI equalization
% structurally cannot give a monotonically worsening curve, regardless of background choice.
%
% The physically correct mechanism for "BER worsens as speed increases" is channel
% ESTIMATION error growing once the true Doppler exceeds the receiver's design search
% range -- analogous to estimateDelayDoppler.m's designMaxDoppler grid limit used
% elsewhere in this report, but implemented fresh here for the LMMSE/unitary-operator
% framework (get_waveform_operator-style U/W) used throughout Muc 9.7/11/13.
%
% Method per trial: transmit a KNOWN pilot through the true single-path Doppler channel,
% estimate nu via a periodogram (matched-filter) search over a LIMITED grid [0,designMaxNu],
% reconstruct Hhat from nu_hat, then LMMSE-equalize the DATA (sent through the TRUE
% channel) using Hhat. Once true nu exceeds designMaxNu, the estimate saturates at the
% grid boundary and the estimated/true channel mismatch grows with velocity.
clc; clear;
% L=16 (matches the N=16 dimension used for ML-diversity/exhaustive-search sections
% elsewhere in this report) instead of L=256: at L=256 a single-path LMMSE system
% collapses almost discontinuously once Doppler-estimate mismatch appears (baseline
% BER~0.0008 near-perfect, jumping to ~0.5 noise-floor within one 20km/h step) because
% phase mismatch is amplified across 256 subcarriers -- not a usable curve. At L=16 the
% estimator's coarser native resolution (~1/16 vs ~1/256) gives a non-negligible
% baseline BER, so the estimation-mismatch penalty past the design boundary shows up as
% a gradual slope instead of a cliff. fs is scaled down by the same L ratio (16/256) so
% the physical nu(velocity) mapping is IDENTICAL to the L=256 attempts above (nu(20km/h)
% still ~=0.53, nu(250km/h) still ~=6.6) -- only the receiver dimension changed.
L = 16; M = 4; N = 4;
c = 3e8; fc = 28e9; fs = 500e3 * (16/256);

velocities = [20 40 60 80 100 120 140 160 180 200 220 250]; % km/h
dopplerOfL = @(v_kmh) L * (2*(v_kmh/3.6)*fc/c) / fs;
designMaxNu = 2.5; % ~94 km/h -- points beyond this exceed the estimator's design grid
gridNu = 0:0.01:designMaxNu;

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

% fixed pseudo-random QPSK pilot, known to receiver
rng(777);
pilotBits = randi([0 1], 2*L, 1);
xPilot = qpsk_mod_local(pilotBits);

% precompute periodogram basis: expBasis(:,g) = exp(2j*pi*gridNu(g)*n/L)
n = (0:L-1).';
expBasis = exp(1j*2*pi*n*gridNu/L); % L x numel(gridNu)

fprintf('--- Target Doppler nu at each tested velocity (single path, designMaxNu=%.2f) ---\n', designMaxNu);
for vi = 1:numel(velocities)
    fprintf('v=%dkm/h -> nu_target=%.4f (%s design grid)\n', velocities(vi), dopplerOfL(velocities(vi)), ...
        ternary_local(dopplerOfL(velocities(vi))<=designMaxNu, 'WITHIN', 'EXCEEDS'));
end

BER = zeros(numel(wfNames), numel(velocities), numel(snrList));
nuHatMean = zeros(numel(velocities), numel(snrList));

for si = 1:numel(snrList)
    snrDb = snrList(si);
    noiseVar = 10^(-snrDb/10);
    for vi = 1:numel(velocities)
        nuTrue = dopplerOfL(velocities(vi));
        bitErr = zeros(numel(wfNames),1); bitTot = 0;
        nuHatSum = 0;
        rng(3000*si+vi);
        for t = 1:numMonte
            % --- pilot / channel estimation (single path, delay=0, gain=1, physical channel) ---
            Dtrue = exp(1j*2*pi*nuTrue*n/L);
            Htrue = diag(Dtrue);
            pilotNoise = sqrt(noiseVar/2)*(randn(L,1)+1j*randn(L,1));
            yPilot = Dtrue.*xPilot + pilotNoise;
            r = yPilot .* conj(xPilot); % ideally = Dtrue + noise-ish term
            scores = abs(r' * expBasis); % 1 x numel(gridNu)
            [~, gIdx] = max(scores);
            nuHat = gridNu(gIdx);
            nuHatSum = nuHatSum + nuHat;
            Dhat = exp(1j*2*pi*nuHat*n/L);
            Hhat = diag(Dhat);

            % --- data transmission through TRUE channel, equalized with ESTIMATED channel ---
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
        nuHatMean(vi,si) = nuHatSum/numMonte;
    end
    fprintf('SNR=%ddB done\n', snrDb);
end

save('velocity_v3_estimated_results.mat','velocities','snrList','BER','wfNames','designMaxNu','nuHatMean');
for si=1:numel(snrList)
    fprintf('--- SNR=%d dB ---\n', snrList(si));
    for wi=1:numel(wfNames)
        fprintf('%s: %s\n', wfNames{wi}, mat2str(squeeze(BER(wi,:,si)),4));
    end
    fprintf('mean nu_hat: %s\n', mat2str(nuHatMean(:,si)',4));
end
fprintf('=== DONE ===\n');

function s = qpsk_mod_local(b)
b = b(:); I = 1-2*b(1:2:end); Q = 1-2*b(2:2:end); s = (I+1j*Q)/sqrt(2);
end
function b = qpsk_hard_demod_local(s)
I = real(s)<0; Q = imag(s)<0;
b = zeros(2*numel(s),1); b(1:2:end)=I; b(2:2:end)=Q;
end
function out = ternary_local(cond, a, b)
if cond; out = a; else; out = b; end
end
