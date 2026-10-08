% Kiểm chứng trực tiếp bằng số các giả định về điều kiện SNR trong mô
% phỏng cảm biến (Muc 12), thay vì chỉ suy luận qua công thức.
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

fprintf('--- 1) Cong suat trung binh moi mau cua tin hieu phat txSamples (ky vong = 1?) ---\n');
rng(7);
numTest = 200;
for k = 1:4
    p = 0;
    for t = 1:numTest
        bitTx = randi([0 1], 2*L, 1);
        symTx = qpsk_mod(bitTx);
        txSamples = ops{k}.tx * symTx;
        p = p + mean(abs(txSamples).^2);
    end
    fprintf('%-5s : cong suat trung binh/mau = %.4f\n', waveformList{k}, p/numTest);
end

fprintf('\n--- 2) Nang luong echo muc tieu (H_radar*tx) so voi nang luong tin hieu phat ---\n');
bitTx = randi([0 1], 2*L, 1);
symTx = qpsk_mod(bitTx);
txSamples = ops{2}.tx * symTx; % OFDM lam vi du
echo = H_radar * txSamples;
fprintf('cong suat tx = %.4f, cong suat echo (khong nhieu) = %.4f (ty le = %.4f, %.2f dB)\n', ...
    mean(abs(txSamples).^2), mean(abs(echo).^2), mean(abs(echo).^2)/mean(abs(txSamples).^2), ...
    10*log10(mean(abs(echo).^2)/mean(abs(txSamples).^2)));

fprintf('\n--- 3) Do loi xu ly thuc te cua CAF: so sanh dinh CAF (co tin hieu) voi nen CAF (chi nhieu) ---\n');
snrDb_test = 0; % nhan SNR moi mau danh gia
noiseVar = 10^(-snrDb_test/10);
delayGrid = 0:24; dopplerGrid = -8:8;
numTrialsVerify = 200;
peakPowerAcc = 0; noisePowerAcc = 0; noiseCount = 0;
rng(11);
for t = 1:numTrialsVerify
    bitTx = randi([0 1], 2*L, 1);
    symTx = qpsk_mod(bitTx);
    txSamples = ops{2}.tx * symTx; % OFDM
    noise = sqrt(noiseVar/2)*(randn(L,1)+1j*randn(L,1));
    rxRadar = H_radar*txSamples + noise;
    CAF = compute_caf(rxRadar, txSamples, delayGrid, dopplerGrid);
    CAFabs2 = abs(CAF).^2;
    % dinh that (biet truoc vi tri muc tieu manh nhat: delay=5,doppler=3)
    [~, idop] = min(abs(dopplerGrid-3)); [~, idel] = min(abs(delayGrid-5));
    peakPowerAcc = peakPowerAcc + CAFabs2(idop, idel);
    % nen: lay trung binh cac o CACH XA moi muc tieu (goc phai tren, vd delay>=15,doppler<=-6)
    bg = CAFabs2(dopplerGrid<=-6, delayGrid>=15);
    noisePowerAcc = noisePowerAcc + mean(bg(:));
    noiseCount = noiseCount + 1;
end
peakPower = peakPowerAcc/numTrialsVerify;
noisePower = noisePowerAcc/noiseCount;
effSNR_dB = 10*log10(peakPower/noisePower);
fprintf('SNR nhan (moi mau) = %d dB\n', snrDb_test);
fprintf('Cong suat dinh CAF (muc tieu manh, trung binh %d trial) = %.3f\n', numTrialsVerify, peakPower);
fprintf('Cong suat nen CAF (chi nhieu, xa muc tieu) = %.3f\n', noisePower);
fprintf('SNR hieu dung tai dinh CAF (do truc tiep, KHONG qua cong thuc) = %.2f dB\n', effSNR_dB);
fprintf('Do loi xu ly do duoc = %.2f dB (so voi cong thuc 10log10(L)=%.2f dB)\n', effSNR_dB - snrDb_test, 10*log10(L));

fprintf('\n=== DONE ===\n');
