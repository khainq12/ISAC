%% ========================================================================
% main_isac_compare_refined_finegrid.m
% Fine SNR grid (1.0:0.2:4.0 dB) with much higher numMonte (24000) to
% shrink the finite-trial "outlier count" staircase visible in the
% detection-threshold region, on top of the unmodified reference
% compute_caf.m + sub-grid parabolic refinement (refine_caf_peak_2d.m).
% ========================================================================

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
M = params.M;
N = params.N;

waveformList = {'RM','OFDM','OTFS','AFDM'};
numWave = numel(waveformList);

snrVec = 1.0:0.2:4.0;
numSNR = numel(snrVec);

numMonte = 16000;
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

BER_mat        = zeros(numWave, numSNR);
RangeNMSE_mat  = zeros(numWave, numSNR);
VelNMSE_mat    = zeros(numWave, numSNR);

fprintf('================ Monte Carlo (fine grid, high numMonte) ================\n');
fprintf('L = %d, M = %d, N = %d\n', L, M, N);
fprintf('SNR sweep = [%s] dB\n', num2str(snrVec));
fprintf('numMonte = %d\n', numMonte);
fprintf('==========================================================================\n\n');

for iSNR = 1:numSNR
    snrDb = snrVec(iSNR);
    fprintf('---- SNR = %.2f dB (point %d/%d) ----\n', snrDb, iSNR, numSNR);

    bitErrCnt      = zeros(numWave,1);
    bitTotalCnt    = zeros(numWave,1);
    rangeNMSEAcc   = zeros(numWave,1);
    velNMSEAcc     = zeros(numWave,1);

    rng(1,"twister");

    for imc = 1:numMonte
        bitTx = randi([0,1], 2*L, 1);
        symTx = qpsk_mod(bitTx);

        for k = 1:numWave
            op = ops{k};

            commRes = simulate_comm_link(symTx, bitTx, op, H_comm, snrDb);
            bitErrCnt(k)   = bitErrCnt(k) + sum(commRes.bitHat ~= bitTx);
            bitTotalCnt(k) = bitTotalCnt(k) + numel(bitTx);

            radarRes = simulate_radar_link(commRes.txSamples, H_radar, snrDb, params.delayGrid, params.dopplerGrid);

            CAFabs = abs(radarRes.CAF);
            CAFwork = CAFabs;
            CAFabs2 = CAFabs.^2;

            estDelays   = zeros(1, numTargets);
            estDopplers = zeros(1, numTargets);

            for it = 1:numTargets
                [~, idxMax] = max(CAFwork(:));
                [rowMax, colMax] = ind2sub(size(CAFwork), idxMax);

                [delayRef, dopplerRef] = refine_caf_peak_2d(CAFabs2, rowMax, colMax, params.dopplerGrid, params.delayGrid);
                estDopplers(it) = dopplerRef;
                estDelays(it)   = delayRef;

                r1 = max(1, rowMax - guardDoppler);
                r2 = min(size(CAFwork,1), rowMax + guardDoppler);
                c1 = max(1, colMax - guardDelay);
                c2 = min(size(CAFwork,2), colMax + guardDelay);
                CAFwork(r1:r2, c1:c2) = -inf;
            end

            estRanges = estDelays * rangeRes;
            estVels   = estDopplers * velRes;

            if numTargets == 2
                estRanges_case1 = estRanges([1 2]);
                estVels_case1   = estVels([1 2]);
                estRanges_case2 = estRanges([2 1]);
                estVels_case2   = estVels([2 1]);

                cost1 = norm(estRanges_case1 - trueRangesSorted)^2 / (norm(trueRangesSorted)^2 + eps) + ...
                        norm(estVels_case1   - trueVelsSorted)^2   / (norm(trueVelsSorted)^2 + eps);
                cost2 = norm(estRanges_case2 - trueRangesSorted)^2 / (norm(trueRangesSorted)^2 + eps) + ...
                        norm(estVels_case2   - trueVelsSorted)^2   / (norm(trueVelsSorted)^2 + eps);

                if cost1 <= cost2
                    estRangesMatched = estRanges_case1;
                    estVelsMatched   = estVels_case1;
                else
                    estRangesMatched = estRanges_case2;
                    estVelsMatched   = estVels_case2;
                end
            else
                [estRangesMatched, idxEstSort] = sort(estRanges, 'ascend');
                estVelsMatched = estVels(idxEstSort);
            end

            rangeNMSE_this = norm(estRangesMatched - trueRangesSorted)^2 / (norm(trueRangesSorted)^2 + eps);
            velNMSE_this   = norm(estVelsMatched   - trueVelsSorted)^2   / (norm(trueVelsSorted)^2 + eps);

            rangeNMSEAcc(k) = rangeNMSEAcc(k) + rangeNMSE_this;
            velNMSEAcc(k)   = velNMSEAcc(k)   + velNMSE_this;
        end
    end

    BER_mat(:, iSNR)       = bitErrCnt ./ bitTotalCnt;
    RangeNMSE_mat(:, iSNR) = rangeNMSEAcc / numMonte;
    VelNMSE_mat(:, iSNR)   = velNMSEAcc   / numMonte;

    for k = 1:numWave
        fprintf('%-5s : BER = %.4e, Range-NMSE = %.4e, Vel-NMSE = %.4e\n', ...
            waveformList{k}, BER_mat(k,iSNR), RangeNMSE_mat(k,iSNR), VelNMSE_mat(k,iSNR));
    end
    fprintf('\n');
end

results_snr = struct();
results_snr.waveformList   = waveformList;
results_snr.snrVec         = snrVec;
results_snr.numMonte       = numMonte;
results_snr.BER_mat        = BER_mat;
results_snr.RangeNMSE_mat  = RangeNMSE_mat;
results_snr.VelNMSE_mat    = VelNMSE_mat;
results_snr.trueRanges     = trueRangesSorted;
results_snr.trueVels       = trueVelsSorted;
results_snr.rangeRes       = rangeRes;
results_snr.velRes         = velRes;

assignin('base', 'results_isac_demo_snr', results_snr);

fprintf('================ DONE ================\n');
