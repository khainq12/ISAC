outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
A = load([outdir 'refined_snr_sweep_results.mat']);
B = load([outdir 'refined_extra_snr_sweep_results.mat']);

wf = A.results_snr.waveformList; % same order in both

snrVec = [A.results_snr.snrVec, B.results_snr.snrVec];
RangeNMSE = [A.results_snr.RangeNMSE_mat, B.results_snr.RangeNMSE_mat];
VelNMSE   = [A.results_snr.VelNMSE_mat,   B.results_snr.VelNMSE_mat];

[snrVec, sortIdx] = sort(snrVec);
RangeNMSE = RangeNMSE(:, sortIdx);
VelNMSE   = VelNMSE(:, sortIdx);

order = [find(strcmp(wf,'AFDM')) find(strcmp(wf,'OFDM')) find(strcmp(wf,'OTFS')) find(strcmp(wf,'RM'))];
names = {'AFDM','OFDM','OTFS','RM (tham chiếu)'};
colors = [0.466 0.674 0.188; 0 0.447 0.741; 0.85 0.325 0.098; 0.5 0.5 0.5];
mk = {'^-','s-','d-','o-'};

plotOne(snrVec, RangeNMSE, order, names, colors, mk, ...
    'Range NMSE (chuẩn hóa)', 'Sai số ước lượng khoảng cách (NMSE) theo SNR - có nội suy dưới-lưới', outdir, 'repo_rangenmse_merged_v1.png');
plotOne(snrVec, VelNMSE, order, names, colors, mk, ...
    'Velocity NMSE (chuẩn hóa)', 'Sai số ước lượng vận tốc (NMSE) theo SNR - có nội suy dưới-lưới', outdir, 'repo_velnmse_merged_v1.png');

% zoomed view of the threshold region for closer inspection
zoomIdx = snrVec >= -1 & snrVec <= 4;
plotOne(snrVec(zoomIdx), RangeNMSE(:,zoomIdx), order, names, colors, mk, ...
    'Range NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Range)', outdir, 'repo_rangenmse_zoom_v1.png');
plotOne(snrVec(zoomIdx), VelNMSE(:,zoomIdx), order, names, colors, mk, ...
    'Velocity NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Velocity)', outdir, 'repo_velnmse_zoom_v1.png');

save([outdir 'merged_snr_sweep_results.mat'], 'snrVec', 'wf', 'RangeNMSE', 'VelNMSE');
fprintf('=== DONE ===\n');

function plotOne(snrVec, dataMat, order, names, colors, mk, yLbl, titleStr, outdir, fname)
    fig = figure('Position',[100 100 700 440],'Visible','off');
    hold on;
    for i = 1:4
        y = dataMat(order(i), :);
        plot(snrVec, y, mk{i}, 'Color', colors(i,:), 'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', names{i});
    end
    set(gca, 'YScale', 'log');
    hold off; grid on; box on;
    xlabel('SNR [dB]'); ylabel(yLbl);
    title(titleStr);
    legend('Location','northeast','FontSize',9);
    saveas(fig, [outdir fname]);
end
