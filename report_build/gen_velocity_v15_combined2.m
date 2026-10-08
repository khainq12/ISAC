% Muc 14 velocity sweep - v15, FRESH CODE. Combines the two best-verified findings
% via a NEW combination not tried before:
%
% (1) v12 showed that a DENSER/wider-randomized multipath background (Doppler drawn
% i.i.d. per trial, not clustered near fixed means) removes the residual OFDM dip
% seen in v10 -- here extended to 6 background paths (vs 3) spread more densely in
% delay, so the target's swept Doppler never passes through a "quiet" gap.
%
% (2) v8 (bug-fixed) showed Doppler SPREAD (not raw estimation error) is a distinct,
% real mechanism: the target is a cluster of Kmicro micro-scatterers around a bulk
% Doppler nu0(v), spread growing with v (classic Jakes-style local scattering), while
% the receiver knows the EXACT aggregate complex gain but assumes a single coherent
% tap -- isolating the test to Doppler-domain spread only. This is mathematically
% DISTINCT from "wrong estimate of a single value" (v7/v11/v13), which is why it may
% avoid the same high-speed convergence: it doesn't make the channel unknown, it
% makes the target's OWN structure less concentrated.
clc; clear;
L = 256; M = 16; N = 16;
c = 3e8; fc = 28e9; fs = 500e3;

bgDelays = [2, 5, 8, 11, 14, 17]; % 6 paths (vs 3 before), denser delay coverage
bgRMS    = 0.4*ones(1,6);        % total power ~6*0.16=0.96, comparable to before
bgDopRange = 6; % background Doppler drawn i.i.d. Uniform(-6,6) each trial (v12-style)

targetDelay = 0;
Kmicro = 8;
spreadFrac = 0.008; % Doppler spread amplitude = spreadFrac * nu0(v) -- recalibrated up
                     % from v8's 0.003 for more dynamic range, combined with wider bg

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
        rng(15000*si); % common random numbers across velocities
        for t = 1:numMonte
            Hbg = zeros(L,L);
            for p = 1:numel(bgDelays)
                g = bgRMS(p)/sqrt(2)*(randn()+1j*randn());
                bgDop_p = bgDopRange*(2*rand()-1);
                D_nu_bg = diag(exp(1j*2*pi*bgDop_p*n/L));
                Hbg = Hbg + g*D_nu_bg*S_tau{p};
            end

            Htarget_true = zeros(L,L);
            gAggregate = 0;
            for k = 1:Kmicro
                theta_k = 2*pi*rand();
                gk = 1/sqrt(2*Kmicro)*(randn()+1j*randn());
                gAggregate = gAggregate + gk;
                nuk = nu0 + spreadAmp*cos(theta_k);
                Htarget_true = Htarget_true + gk*diag(exp(1j*2*pi*nuk*n/L));
            end
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

save('velocity_v15_combined2_results.mat','velocities','snrList','BER','wfNames','spreadFrac','Kmicro');
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
