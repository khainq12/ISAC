outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
load([outdir 'merged_snr_sweep_results_v2.mat']);  % snrVec, wf, RangeNMSE, VelNMSE (from v2)

order = [find(strcmp(wf,'AFDM')) find(strcmp(wf,'OFDM')) find(strcmp(wf,'OTFS')) find(strcmp(wf,'RM'))];
names = {'AFDM','OFDM','OTFS','RM (tham chiếu)'};
colors = [0.302 0.686 0.290; 0.216 0.494 0.722; 0.894 0.102 0.110; 0.5 0.5 0.5];

% Light Gaussian-weighted local smoothing IN LOG10(NMSE), window = 0.3 dB,
% applied only to average out the single-trial-outlier-count staircase in
% the finite-numMonte threshold region. Window is much narrower than any
% genuine cross-waveform feature scale (>=0.3-0.5dB), so real ordering/
% crossover differences between waveforms are preserved, not erased -
% unlike the earlier rejected smoothingspline attempt (SmoothingParam
% 0.35/0.9 on much sparser data, which washed out OTFS's early transition).
windowDb = 0.3;

plotOne(snrVec, RangeNMSE, order, names, colors, windowDb, ...
    'Range NMSE (chuẩn hóa)', 'Sai số ước lượng khoảng cách (NMSE) theo SNR', outdir, 'repo_rangenmse_polished_v3.png');
plotOne(snrVec, VelNMSE, order, names, colors, windowDb, ...
    'Velocity NMSE (chuẩn hóa)', 'Sai số ước lượng vận tốc (NMSE) theo SNR', outdir, 'repo_velnmse_polished_v3.png');

zoomIdx = snrVec >= -1 & snrVec <= 4;
plotOne(snrVec(zoomIdx), RangeNMSE(:,zoomIdx), order, names, colors, windowDb, ...
    'Range NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Range)', outdir, 'repo_rangenmse_zoompolished_v3.png');
plotOne(snrVec(zoomIdx), VelNMSE(:,zoomIdx), order, names, colors, windowDb, ...
    'Velocity NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Velocity)', outdir, 'repo_velnmse_zoompolished_v3.png');

fprintf('=== DONE ===\n');

function plotOne(snrVec, dataMat, order, names, colors, windowDb, yLbl, titleStr, outdir, fname)
    fig = figure('Position',[100 100 1150 700],'Visible','off','Color','w');
    ax = axes(fig); hold(ax, 'on');
    xq = linspace(min(snrVec), max(snrVec), 600);
    for i = 1:4
        y = dataMat(order(i), :);
        logy = log10(y);
        logySmooth = smoothdata(logy, 'gaussian', windowDb, 'SamplePoints', snrVec);
        logyq = pchip(snrVec, logySmooth, xq);
        plot(ax, xq, 10.^logyq, '-', 'Color', colors(i,:), 'LineWidth', 3.2, 'DisplayName', names{i});
    end
    set(ax, 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1.1, 'TickDir', 'out');
    ax.GridAlpha = 0.15; ax.MinorGridAlpha = 0.08;
    grid(ax, 'on'); grid(ax, 'minor'); box(ax, 'on');
    xlabel(ax, 'SNR [dB]', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel(ax, yLbl, 'FontSize', 16, 'FontWeight', 'bold');
    title(ax, titleStr, 'FontSize', 17, 'FontWeight', 'bold');
    lg = legend(ax, 'Location','northeast','FontSize',13);
    lg.Box = 'on'; lg.EdgeColor = [0.7 0.7 0.7];
    xlim(ax, [min(snrVec) max(snrVec)]);
    exportgraphics(fig, [outdir fname], 'Resolution', 220);
    close(fig);
end
