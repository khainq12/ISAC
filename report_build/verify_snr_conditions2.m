addpath('D:\doubly-dispersive-channel-simulation-main\rm_isac_demo');

params = struct();
params.M = 16; params.N = 16; params.L = params.M*params.N;
params.target.delays   = [5, 18];
params.target.dopplers = [3, -4];
params.target.gains    = [0.68*exp(1j*0.2), 0.33*exp(1j*0.8)];
params.afdm.c1 = 5/params.L; params.afdm.c2 = 1/(params.L^2);
params.rmPermutationSeed = 314159;
L = params.L;

waveformList = {'RM','OFDM','OTFS','AFDM'};
ops = cell(4,1);
for k=1:4
    ops{k} = get_waveform_operator(waveformList{k}, params);
end

H_radar = build_delay_doppler_channel(L, params.target.delays, params.target.dopplers, params.target.gains);
delayGrid = 0:24; dopplerGrid = -8:8;
numTrialsVerify = 300;

snrTestList = [-3, 0];

for si = 1:numel(snrTestList)
    snrDb_test = snrTestList(si);
    noiseVar = 10^(-snrDb_test/10);
    fprintf('=== SNR nhan = %d dB ===\n', snrDb_test);
    for k = 1:4
        peakPowerAcc = 0; noisePowerAcc = 0;
        rng(100+k);
        for t = 1:numTrialsVerify
            bitTx = randi([0 1], 2*L, 1);
            symTx = qpsk_mod(bitTx);
            txSamples = ops{k}.tx * symTx;
            noise = sqrt(noiseVar/2)*(randn(L,1)+1j*randn(L,1));
            rxRadar = H_radar*txSamples + noise;
            CAF = compute_caf(rxRadar, txSamples, delayGrid, dopplerGrid);
            CAFabs2 = abs(CAF).^2;
            [~, idop] = min(abs(dopplerGrid-3)); [~, idel] = min(abs(delayGrid-5));
            peakPowerAcc = peakPowerAcc + CAFabs2(idop, idel);
            bg = CAFabs2(dopplerGrid<=-6, delayGrid>=15);
            noisePowerAcc = noisePowerAcc + mean(bg(:));
        end
        peakPower = peakPowerAcc/numTrialsVerify;
        noisePower = noisePowerAcc/numTrialsVerify;
        effSNR_dB = 10*log10(peakPower/noisePower);
        fprintf('%-5s : SNR hieu dung tai dinh CAF = %.2f dB (do loi = %.2f dB)\n', ...
            waveformList{k}, effSNR_dB, effSNR_dB - snrDb_test);
    end
    fprintf('\n');
end
fprintf('=== DONE ===\n');
