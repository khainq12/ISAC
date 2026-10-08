% Muc 14 velocity sweep - v8, FRESH CODE, DIFFERENT MECHANISM from v7.
% v7 modeled target-Doppler ESTIMATION error growing past a design boundary -- but
% every variant of that mechanism showed all 4 waveforms CONVERGE at large mismatch,
% because once the receiver's channel assumption is simply WRONG (not just noisy),
% the situation degenerates toward the same "unknown channel" regime where the
% unitary-invariance theorem (already proven in Muc 9.2/11 for genie-CSI) reasserts
% itself -- no waveform structure can compensate for information that isn't there.
%
% v8 instead models classical DOPPLER SPREAD (Jakes-style local scattering): the
% receiver has EXACT/genie knowledge of the target's BULK Doppler (no tracking-error
% at all this time) -- but the true physical return isn't one coherent tap, it's a
% cluster of K micro-paths (local scatterers around the target) whose Doppler SPREAD
% around the bulk value grows proportionally with speed (nuSpread(v) = frac*nu0(v)),
% while the receiver still equalizes assuming a single coherent tap at nu0(v). This
% is structurally the SAME kind of problem as the loc-collision mechanism already
% proven in Muc 12.2 (Theorem 1): multiple delay-Doppler components landing too close
% together for a rank-1 model to resolve -- so AFDM's proven chirp-based position-
% separation should give it a REAL, non-converging advantage here, unlike v7.
clc; clear;
L = 256; M = 16; N = 16;
c = 3e8; fc = 28e9; fs = 500e3;

bgDelays   = [3, 11, 17];
bgRMS      = [0.75, 0.55, 0.38];
bgDopMean  = [1, -2, 3];
bgDopSpread = 1.5;
targetDelay = 0;
Kmicro = 8;          % micro-scatterers composing the target return
spreadFrac = 0.003;  % Doppler spread amplitude = spreadFrac * nu0(v) -- L=256 is extremely
                      % sensitive to normalized-Doppler mismatch (confirmed in v3/v4/v7:
                      % std~0.01-0.02 already causes near-total decorrelation), so this
                      % must be ~100x smaller than a naive Jakes-spread-fraction guess

velocities = [20 40 60 80 100 120 140 160 180 200 220 250]; % km/h
dopplerOfL = @(v_kmh) L * (2*(v_kmh/3.6)*fc/c) / fs;

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
S_tau = cell(1,numel(bgDelays));
for p = 1:numel(bgDelays)
    S_tau{p} = circshift(I, bgDelays(p), 1);
end

fprintf('--- Target nu0 / Doppler-spread amplitude at each tested velocity ---\n');
for vi = 1:numel(velocities)
    nu0 = dopplerOfL(velocities(vi));
    fprintf('v=%dkm/h -> nu0=%.4f, spreadAmp=%.4f\n', velocities(vi), nu0, spreadFrac*nu0);
end

BER = zeros(numel(wfNames), numel(velocities), numel(snrList));

for si = 1:numel(snrList)
    snrDb = snrList(si);
    noiseVar = 10^(-snrDb/10);
    for vi = 1:numel(velocities)
        nu0 = dopplerOfL(velocities(vi));
        spreadAmp = spreadFrac*nu0;
        bitErr = zeros(numel(wfNames),1); bitTot = 0;
        rng(8000*si); % common random numbers across velocities
        for t = 1:numMonte
            Hbg = zeros(L,L);
            for p = 1:numel(bgDelays)
                g = bgRMS(p)/sqrt(2)*(randn()+1j*randn());
                bgDop_p = bgDopMean(p) + bgDopSpread*(2*rand()-1);
                D_nu_bg = diag(exp(1j*2*pi*bgDop_p*n/L));
                Hbg = Hbg + g*D_nu_bg*S_tau{p};
            end

            % true target: cluster of Kmicro micro-paths spread around nu0
            Htarget_true = zeros(L,L);
            gAggregate = 0;
            for k = 1:Kmicro
                theta_k = 2*pi*rand();
                gk = 1/sqrt(2*Kmicro)*(randn()+1j*randn());
                gAggregate = gAggregate + gk;
                nuk = nu0 + spreadAmp*cos(theta_k);
                Htarget_true = Htarget_true + gk*diag(exp(1j*2*pi*nuk*n/L));
            end
            % receiver knows the EXACT aggregate complex gain (realistic AGC/phase
            % tracking) but assumes it all sits at the single bulk Doppler nu0 -- isolates
            % the test to the Doppler-SPREAD structure only, not gain/phase estimation.
            % At spreadAmp=0 this makes Htarget_hat == Htarget_true exactly (sanity check).
            Htarget_hat = gAggregate*diag(exp(1j*2*pi*nu0*n/L));

            Htrue = Hbg + Htarget_true;
            Hhat  = Hbg + Htarget_hat;

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

save('velocity_v8_dopplerspread_results.mat','velocities','snrList','BER','wfNames','spreadFrac','Kmicro');
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
