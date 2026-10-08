clc; clear; close all;

params = struct();
params.M = 16;
params.N = 16;
params.L = params.M * params.N;

params.target.delays   = [5, 18];
params.target.dopplers = [3, -4];
params.target.gains    = [0.68 * exp(1j*0.2),0.33 * exp(1j*0.8)];

params.user.delays   = [0, 3, 11, 17];
params.user.dopplers = [0, 1, -2, 3];
params.user.gains    = [1.00 * exp(1j*0.0), 0.75 * exp(1j*0.4),0.55 * exp(-1j*0.9),0.38 * exp(1j*1.2)];

params.delayGrid   = 0:24;
params.dopplerGrid = -8:8;

params.afdm.c1 = 5 / params.L;
params.afdm.c2 = 1 / (params.L^2);

params.rmPermutationSeed = 314159;

L = params.M * params.N;
waveformList = {'RM','OFDM','OTFS','AFDM'};
numWave = numel(waveformList);

snrDb = 3.0;
numMonte = 2000;
guardDelay = 1;
guardDoppler = 1;

c  = 3e8;
fs = 20e6;
fc = 24e9;
lambda = c / fc;
Ts = 1 / fs;

rangeRes = c * Ts / 2;
velRes   = lambda * fs / (2 * L);

trueDelays   = params.target.delays(:).';
trueDopplers = params.target.dopplers(:).';
numTargets   = numel(trueDelays);

trueRanges = trueDelays * rangeRes;
trueVels   = trueDopplers * velRes;
[trueRangesSorted, idxTrueSort] = sort(trueRanges, 'ascend');
trueVelsSorted = trueVels(idxTrueSort);

H_comm = build_delay_doppler_channel(L, params.user.delays, params.user.dopplers, params.user.gains);
H_radar = build_delay_doppler_channel(L, params.target.delays, params.target.dopplers, params.target.gains);

ops = cell(numWave,1);
for k = 1:numWave
    ops{k} = get_waveform_operator(waveformList{k}, params);
end

tic;
rng(1,"twister");
for imc = 1:numMonte
    bitTx = randi([0,1], 2*L, 1);
    symTx = qpsk_mod(bitTx);
    for k = 1:numWave
        op = ops{k};
        commRes = simulate_comm_link(symTx, bitTx, op, H_comm, snrDb);
        radarRes = simulate_radar_link(commRes.txSamples, H_radar, snrDb, params.delayGrid, params.dopplerGrid);
        CAFabs = abs(radarRes.CAF);
        CAFwork = CAFabs;
        CAFabs2 = CAFabs.^2;
        for it = 1:numTargets
            [~, idxMax] = max(CAFwork(:));
            [rowMax, colMax] = ind2sub(size(CAFwork), idxMax);
            [delayRef, dopplerRef] = refine_caf_peak_2d(CAFabs2, rowMax, colMax, params.dopplerGrid, params.delayGrid);
            r1 = max(1, rowMax - guardDoppler);
            r2 = min(size(CAFwork,1), rowMax + guardDoppler);
            c1 = max(1, colMax - guardDelay);
            c2 = min(size(CAFwork,2), colMax + guardDelay);
            CAFwork(r1:r2, c1:c2) = -inf;
        end
    end
end
elapsed = toc;
fprintf('=== TIMING: %d trials x %d waveforms took %.3f sec => %.5f sec/trial-all-waveforms ===\n', numMonte, numWave, elapsed, elapsed/numMonte);
fprintf('=== DONE ===\n');
