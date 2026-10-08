outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
load([outdir 'merged_snr_sweep_results.mat']);

order = [find(strcmp(wf,'AFDM')) find(strcmp(wf,'OFDM')) find(strcmp(wf,'OTFS')) find(strcmp(wf,'RM'))];
names = {'AFDM','OFDM','OTFS','RM (tham chiếu)'};
colors = [0.466 0.674 0.188; 0 0.447 0.741; 0.85 0.325 0.098; 0.5 0.5 0.5];

plotOne(snrVec, RangeNMSE, order, names, colors, ...
    'Range NMSE (chuẩn hóa)', 'Sai số ước lượng khoảng cách (NMSE) theo SNR', outdir, 'repo_rangenmse_final_v1.png', []);
plotOne(snrVec, VelNMSE, order, names, colors, ...
    'Velocity NMSE (chuẩn hóa)', 'Sai số ước lượng vận tốc (NMSE) theo SNR', outdir, 'repo_velnmse_final_v1.png', []);

zoomIdx = snrVec >= -1 & snrVec <= 4;
plotOne(snrVec(zoomIdx), RangeNMSE(:,zoomIdx), order, names, colors, ...
    'Range NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Range)', outdir, 'repo_rangenmse_zoomfinal_v1.png', []);
plotOne(snrVec(zoomIdx), VelNMSE(:,zoomIdx), order, names, colors, ...
    'Velocity NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Velocity)', outdir, 'repo_velnmse_zoomfinal_v1.png', []);

fprintf('=== DONE ===\n');

function plotOne(snrVec, dataMat, order, names, colors, yLbl, titleStr, outdir, fname, ~)
    fig = figure('Position',[100 100 1100 680],'Visible','off');
    hold on;
    xq = linspace(min(snrVec), max(snrVec), 400);
    for i = 1:4
        y = dataMat(order(i), :);
        % shape-preserving smooth interpolation through the real measured
        % points (pchip) - purely a visual line-smoothing between actual
        % data, does not alter or fabricate any measured value
        logy = log10(y);
        logyq = pchip(snrVec, logy, xq);
        plot(xq, 10.^logyq, '-', 'Color', colors(i,:), 'LineWidth', 2.6, 'DisplayName', names{i});
    end
    set(gca, 'YScale', 'log', 'FontSize', 13);
    hold off; grid on; box on;
    xlabel('SNR [dB]', 'FontSize', 14); ylabel(yLbl, 'FontSize', 14);
    title(titleStr, 'FontSize', 15);
    legend('Location','northeast','FontSize',12);
    exportgraphics(fig, [outdir fname], 'Resolution', 200);
    close(fig);
end
