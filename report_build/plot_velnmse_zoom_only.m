outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
A = load([outdir 'refined_snr_sweep_results.mat']);
B = load([outdir 'refined_extra_snr_sweep_results.mat']);

wf = A.results_snr.waveformList;
snrVec = [A.results_snr.snrVec, B.results_snr.snrVec];
RangeNMSE = [A.results_snr.RangeNMSE_mat, B.results_snr.RangeNMSE_mat];
VelNMSE   = [A.results_snr.VelNMSE_mat,   B.results_snr.VelNMSE_mat];
[snrVec, sortIdx] = sort(snrVec);
RangeNMSE = RangeNMSE(:, sortIdx);
VelNMSE   = VelNMSE(:, sortIdx);

save([outdir 'merged_snr_sweep_results.mat'], 'snrVec', 'wf', 'RangeNMSE', 'VelNMSE');

order = [find(strcmp(wf,'AFDM')) find(strcmp(wf,'OFDM')) find(strcmp(wf,'OTFS')) find(strcmp(wf,'RM'))];
names = {'AFDM','OFDM','OTFS','RM (tham chiếu)'};
colors = [0.466 0.674 0.188; 0 0.447 0.741; 0.85 0.325 0.098; 0.5 0.5 0.5];
mk = {'^-','s-','d-','o-'};

zoomIdx = snrVec >= -1 & snrVec <= 4;

fig = figure('Position',[100 100 700 440],'Visible','off');
hold on;
for i = 1:4
    y = VelNMSE(order(i), zoomIdx);
    plot(snrVec(zoomIdx), y, mk{i}, 'Color', colors(i,:), 'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', names{i});
end
set(gca, 'YScale', 'log');
hold off; grid on; box on;
xlabel('SNR [dB]'); ylabel('Velocity NMSE (chuẩn hóa)');
title('Vùng ngưỡng phóng to (Velocity)');
legend('Location','northeast','FontSize',9);
exportgraphics(fig, [outdir 'repo_velnmse_zoom_v1.png'], 'Resolution', 150);
fprintf('=== DONE ===\n');
