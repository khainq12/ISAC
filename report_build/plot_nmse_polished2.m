outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
A  = load([outdir 'refined_snr_sweep_results.mat']);
B  = load([outdir 'refined_extra_snr_sweep_results.mat']);
C  = load([outdir 'refined_extra2_snr_sweep_results.mat']);
Fg = load([outdir 'refined_finegrid_snr_sweep_results.mat']);

wf = A.results_snr.waveformList;

snrA = A.results_snr.snrVec; RA = A.results_snr.RangeNMSE_mat; VA = A.results_snr.VelNMSE_mat;
snrB = B.results_snr.snrVec; RB = B.results_snr.RangeNMSE_mat; VB = B.results_snr.VelNMSE_mat;
snrC = C.results_snr.snrVec; RC = C.results_snr.RangeNMSE_mat; VC = C.results_snr.VelNMSE_mat;
snrF = Fg.results_snr.snrVec; RF = Fg.results_snr.RangeNMSE_mat; VF = Fg.results_snr.VelNMSE_mat;

% drop old low-numMonte points inside [1.0, 4.0] dB - replaced by the new
% high-numMonte (16000) fine grid (main_isac_compare_refined_finegrid.m)
keepA = ~(snrA >= 1.0 & snrA <= 4.0);
keepB = ~(snrB >= 1.0 & snrB <= 4.0);
keepC = ~(snrC >= 1.0 & snrC <= 4.0);

snrVec    = [snrA(keepA), snrB(keepB), snrC(keepC), snrF];
RangeNMSE = [RA(:,keepA), RB(:,keepB), RC(:,keepC), RF];
VelNMSE   = [VA(:,keepA), VB(:,keepB), VC(:,keepC), VF];

[snrVec, sortIdx] = sort(snrVec);
RangeNMSE = RangeNMSE(:, sortIdx);
VelNMSE   = VelNMSE(:, sortIdx);

save([outdir 'merged_snr_sweep_results_v2.mat'], 'snrVec', 'wf', 'RangeNMSE', 'VelNMSE');

order = [find(strcmp(wf,'AFDM')) find(strcmp(wf,'OFDM')) find(strcmp(wf,'OTFS')) find(strcmp(wf,'RM'))];
names = {'AFDM','OFDM','OTFS','RM (tham chiếu)'};
colors = [0.302 0.686 0.290; 0.216 0.494 0.722; 0.894 0.102 0.110; 0.5 0.5 0.5];

plotOne(snrVec, RangeNMSE, order, names, colors, ...
    'Range NMSE (chuẩn hóa)', 'Sai số ước lượng khoảng cách (NMSE) theo SNR', outdir, 'repo_rangenmse_polished_v2.png', []);
plotOne(snrVec, VelNMSE, order, names, colors, ...
    'Velocity NMSE (chuẩn hóa)', 'Sai số ước lượng vận tốc (NMSE) theo SNR', outdir, 'repo_velnmse_polished_v2.png', []);

zoomIdx = snrVec >= -1 & snrVec <= 4;
plotOne(snrVec(zoomIdx), RangeNMSE(:,zoomIdx), order, names, colors, ...
    'Range NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Range)', outdir, 'repo_rangenmse_zoompolished_v2.png', []);
plotOne(snrVec(zoomIdx), VelNMSE(:,zoomIdx), order, names, colors, ...
    'Velocity NMSE (chuẩn hóa)', 'Vùng ngưỡng phóng to (Velocity)', outdir, 'repo_velnmse_zoompolished_v2.png', []);

fprintf('=== DONE ===\n');

function plotOne(snrVec, dataMat, order, names, colors, yLbl, titleStr, outdir, fname, ~)
    fig = figure('Position',[100 100 1150 700],'Visible','off','Color','w');
    ax = axes(fig); hold(ax, 'on');
    xq = linspace(min(snrVec), max(snrVec), 600);
    for i = 1:4
        y = dataMat(order(i), :);
        logy = log10(y);
        logyq = pchip(snrVec, logy, xq);
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
