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

waveformList = {'OFDM','OTFS','AFDM'};
numWave = numel(waveformList);
ops = cell(numWave,1);
for k=1:numWave, ops{k} = get_waveform_operator(waveformList{k}, params); end

guardDelay=1; guardDoppler=1;
c=3e8; fs=20e6; fc=24e9; lambda=c/fc; Ts=1/fs;
rangeRes = c*Ts/2; velRes = lambda*fs/(2*L);

H_comm = build_delay_doppler_channel(L, params.user.delays, params.user.dopplers, params.user.gains);

bitTx0 = []; symTx0 = [];
rng(1,'twister');
bitTx0 = randi([0,1], 2*L, 1);
symTx0 = qpsk_mod(bitTx0);

fprintf('======================================================\n');
fprintf('TEST: fractional-DELAY true target (windowed-sinc), noise = 0\n');
fprintf('(uses build_delay_doppler_channel_frac.m, borrowed technique from\n');
fprintf(' github.com/eric-hs-rou/fractional_AF_generator - our OWN addition,\n');
fprintf(' not part of the original reference repo)\n');
fprintf('======================================================\n');
fracDelays = [5.37, 18.41];
H_radar_frac = build_delay_doppler_channel_frac(L, fracDelays, params.target.dopplers, params.target.gains);
trueRangesFrac = fracDelays*rangeRes;
[trueRangesFracSorted, idxTrueSort] = sort(trueRangesFrac, 'ascend');
trueVelsSorted = params.target.dopplers*velRes; trueVelsSorted = trueVelsSorted(idxTrueSort);

for k = 1:numWave
    op = ops{k};
    commRes = simulate_comm_link(symTx0, bitTx0, op, H_comm, 100);
    txSamples = commRes.txSamples;
    rxRadar = H_radar_frac * txSamples;  % NO noise
    CAF = compute_caf(rxRadar, txSamples, params.delayGrid, params.dopplerGrid);
    CAFabs2 = abs(CAF).^2;
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
    rangeNMSE = norm(estRangesSorted-trueRangesFracSorted)^2/(norm(trueRangesFracSorted)^2+eps);
    velNMSE   = norm(estVelsSorted-trueVelsSorted)^2/(norm(trueVelsSorted)^2+eps);
    fprintf('%-5s : estDelays=[%s] (true frac=[%s])  RangeNMSE=%.3e   VelNMSE(check,should be tiny)=%.3e\n', ...
        waveformList{k}, num2str(sort(estDelays)), num2str(sort(fracDelays)), rangeNMSE, velNMSE);
end

fprintf('\n=== Reference comparison (from earlier sessions) ===\n');
fprintf('On-grid delay, noise=0        : RangeNMSE ~1e-7 to 1e-10\n');
fprintf('Fractional Doppler (3.37/-4.41): VelNMSE ~1e-3 to 1e-2 (earlier test)\n');
fprintf('=== DONE ===\n');
