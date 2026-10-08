outdir = 'D:\doubly-dispersive-channel-simulation-main\report_build\';
load([outdir 'smooth_snr_sweep_results.mat']);

snrVec = results_snr.snrVec;
wf = results_snr.waveformList; % RM, OFDM, OTFS, AFDM
RangeNMSE = results_snr.RangeNMSE_mat;
VelNMSE   = results_snr.VelNMSE_mat;

order = [find(strcmp(wf,'AFDM')) find(strcmp(wf,'OFDM')) find(strcmp(wf,'OTFS')) find(strcmp(wf,'RM'))];
names = {'AFDM','OFDM','OTFS','RM (tham chiếu)'};
colors = [0.466 0.674 0.188; 0 0.447 0.741; 0.85 0.325 0.098; 0.5 0.5 0.5];
mk = {'^-','s-','d-','o-'};

floorVal = 1e-5;

plotOne(snrVec, RangeNMSE, order, names, colors, mk, floorVal, ...
    'Range NMSE (chuẩn hóa)', 'Sai số ước lượng khoảng cách (NMSE) theo SNR - 8000 lần MC', outdir, 'repo_rangenmse_smooth_v1.png');
plotOne(snrVec, VelNMSE, order, names, colors, mk, floorVal, ...
    'Velocity NMSE (chuẩn hóa)', 'Sai số ước lượng vận tốc (NMSE) theo SNR - 8000 lần MC', outdir, 'repo_velnmse_smooth_v1.png');

fprintf('=== DONE ===\n');

function plotOne(snrVec, dataMat, order, names, colors, mk, floorVal, yLbl, titleStr, outdir, fname)
    fig = figure('Position',[100 100 700 440],'Visible','off');
    hold on;
    for i = 1:4
        y = dataMat(order(i), :);
        yPlot = max(y, floorVal);
        plot(snrVec, yPlot, mk{i}, 'Color', colors(i,:), 'LineWidth', 1.6, 'MarkerSize', 5, 'DisplayName', names{i});
    end
    set(gca, 'YScale', 'log');
    hold off; grid on; box on;
    xlabel('SNR [dB]'); ylabel(yLbl);
    title(titleStr);
    legend('Location','northeast','FontSize',9);
    ylim([floorVal/2, 1]);
    saveas(fig, [outdir fname]);
end
