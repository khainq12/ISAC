velocities = [30 60 90 120 150];
snrdB = [5 15 25];
N = 16; maxDelay = 2; guardWidth = 1; designMaxDoppler = 0.08;
carrierFreq = 28e9; samplingFreq = 2e6;
numTrials = 3000;
waveforms = {'OFDM','OTFS','AFDM'};

fprintf('=== Velocity/Doppler ML sweep ===\n');
BER_vel = zeros(3, length(velocities), length(snrdB));
for i = 1:3
    tic;
    BER_vel(i,:,:) = simDiversityMLVelocity(waveforms{i}, velocities, snrdB, N, maxDelay, designMaxDoppler, guardWidth, carrierFreq, samplingFreq, numTrials);
    t = toc;
    fprintf('%s done in %.1f s\n', waveforms{i}, t);
    disp(squeeze(BER_vel(i,:,:)));
end

save('tables_velocity.mat', ...
    'velocities', 'snrdB', 'BER_vel', 'waveforms', 'N', 'maxDelay', 'guardWidth', 'designMaxDoppler', 'carrierFreq', 'samplingFreq', 'numTrials');
fprintf('Velocity sweep saved.\n');
