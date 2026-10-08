clc; clear; close all;

params = struct();
params.M = 16; params.N = 16; params.L = params.M*params.N;
L = params.L;

params.target.delays   = [5, 18];
params.target.dopplers = [3, -4];
params.target.gains    = [0.68*exp(1j*0.2), 0.33*exp(1j*0.8)];

params.user.delays   = [0, 3, 11, 17];
params.user.dopplers = [0, 1, -2, 3];
params.user.gains    = [1.00*exp(1j*0.0), 0.75*exp(1j*0.4), 0.55*exp(-1j*0.9), 0.38*exp(1j*1.2)];

params.delayGrid   = 0:24;
params.dopplerGrid = -8:8;
params.afdm.c1 = 5/L; params.afdm.c2 = 1/(L^2);
params.rmPermutationSeed = 314159;

waveformList = {'AFDM'};
numWave = 1;
ops = {get_waveform_operator('AFDM', params)};

guardDelay=1; guardDoppler=1;
c=3e8; fs=20e6; fc=24e9; lambda=c/fc; Ts=1/fs;
rangeRes = c*Ts/2; velRes = lambda*fs/(2*L);

H_comm = build_delay_doppler_channel(L, params.user.delays, params.user.dopplers, params.user.gains);
H_radar = build_delay_doppler_channel(L, params.target.delays, params.target.dopplers, params.target.gains);
trueDelays = params.target.delays; trueDopplers = params.target.dopplers;
[trueRangesSorted, idxTrueSort] = sort(trueDelays*rangeRes, 'ascend');
trueVelsSorted = trueDopplers*velRes; trueVelsSorted = trueVelsSorted(idxTrueSort);

snrVec = -10:5:40;
numMonte = 3000;

fprintf('======================================================\n');
fprintf('Extended SNR sweep (-10 to 40 dB), AFDM only, numMonte=%d\n', numMonte);
fprintf('======================================================\n');

for snrDb = snrVec
    rangeNMSEAcc = 0; velNMSEAcc = 0;
    rng(1,'twister');
    for imc = 1:numMonte
        bitTx = randi([0,1], 2*L, 1);
        symTx = qpsk_mod(bitTx);
        commRes = simulate_comm_link(symTx, bitTx, ops{1}, H_comm, snrDb);
        radarRes = simulate_radar_link(commRes.txSamples, H_radar, snrDb, params.delayGrid, params.dopplerGrid);
        CAFabs2 = abs(radarRes.CAF).^2;
        CAFwork = CAFabs2;
        estDelays = zeros(1,2); estDopplers = zeros(1,2);
        for it=1:2
            [~,idxMax] = max(CAFwork(:));
            [rowMax,colMax] = ind2sub(size(CAFwork), idxMax);
            [delayRef, dopplerRef] = refine_caf_peak_2d(CAFabs2, rowMax, colMax, params.dopplerGrid, params.delayGrid);
            estDopplers(it) = dopplerRef; estDelays(it) = delayRef;
            r1=max(1,rowMax-guardDoppler); r2=min(size(CAFwork,1),rowMax+guardDoppler);
            c1=max(1,colMax-guardDelay);   c2=min(size(CAFwork,2),colMax+guardDelay);
            CAFwork(r1:r2,c1:c2) = -inf;
        end
        estRanges = estDelays*rangeRes; estVels = estDopplers*velRes;
        [estRangesSorted, idxE] = sort(estRanges,'ascend'); estVelsSorted = estVels(idxE);
        rangeNMSEAcc = rangeNMSEAcc + norm(estRangesSorted-trueRangesSorted)^2/(norm(trueRangesSorted)^2+eps);
        velNMSEAcc   = velNMSEAcc   + norm(estVelsSorted-trueVelsSorted)^2/(norm(trueVelsSorted)^2+eps);
    end
    fprintf('SNR=%4d dB : RangeNMSE=%.4e  VelNMSE=%.4e\n', snrDb, rangeNMSEAcc/numMonte, velNMSEAcc/numMonte);
end

fprintf('=== DONE ===\n');
