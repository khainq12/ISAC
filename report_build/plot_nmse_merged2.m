outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
load([outdir 'merged_snr_sweep_results.mat']);

order = [find(strcmp(wf,'AFDM')) find(strcmp(wf,'OFDM')) find(strcmp(wf,'OTFS')) find(strcmp(wf,'RM'))];
names = {'AFDM','OFDM','OTFS','RM (tham chiếu)'};
colors = [0.466 0.674 0.188; 0 0.447 0.741; 0.85 0.325 0.098; 0.5 0.5 0.5];
mk = {'^-','s-','d-','o-'};

plotOne(snrVec, RangeNMSE, order, names, colors, mk, ...
    'Range NMSE (chuẩn hóa)', 'Sai số ước lượng khoảng cách (NMSE) theo SNR - có nội suy dưới-lưới', outdir, 'repo_rangenmse_merged_v1.png');
plotOne(snrVec, VelNMSE, order, names, colors, mk, ...
    'Velocity NMSE (chuẩn hóa)', 'Sai số ước lượng vận tốc (NMSE) theo SNR - có nội suy dưới-lưới', outdir, 'repo_velnmse_merged_v1.png');

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
    exportgraphics(fig, [outdir fname], 'Resolution', 150);
    close(fig);
end
