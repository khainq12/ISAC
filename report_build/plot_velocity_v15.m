outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
load([outdir 'velocity_v15_combined2_results.mat']); % velocities, BER, wfNames

colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.5 0.5 0.5];
names = {'OFDM','OTFS','AFDM','RM (tham chiếu)'};
mk = {'s-','d-','^-','o-'};

fig = figure('Position',[100 100 680 420],'Visible','off');
hold on;
for i=1:4
    plot(velocities, squeeze(BER(i,:,1)), mk{i}, 'Color', colors(i,:), 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', names{i});
end
hold off; grid on; box on;
xlabel('Tốc độ mục tiêu [km/h]'); ylabel('BER (thu LMMSE, kênh nền 6 đường + Doppler mục tiêu trải rộng)');
title('BER theo tốc độ mục tiêu, SNR=10dB');
legend('Location','northwest','FontSize',9);
saveas(fig, [outdir 'repo_velocity_ber_v15.png']);
fprintf('Saved.\n=== DONE ===\n');
