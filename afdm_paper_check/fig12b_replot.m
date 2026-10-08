clear; clc; close all;
load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12b_results.mat');

figure('Color','w');
mk = {'-o','-s','-^'};
clr = {[0 0.45 0.74],[0.85 0.33 0.1],[0.47 0.67 0.19]};
for ix=1:numel(xiList)
    semilogy(SNRd_dB, BER_est(ix,:), mk{ix}, 'Color', clr{ix}, 'DisplayName', sprintf('Est, \\xi_\\nu=%d', xiList(ix))); hold on;
    semilogy(SNRd_dB, BER_ideal(ix,:), mk{ix}, 'LineStyle', '--', 'Color', clr{ix}, 'DisplayName', sprintf('Ideal, \\xi_\\nu=%d', xiList(ix)));
end
grid on; xlabel('SNR in dB'); ylabel('BER');
title(sprintf('Fig 12b repro: effect of \\xi_\\nu, N=256, SNRp=40dB, numMonte=%d', numMonte));
legend('show','FontSize',8); ylim([1e-2 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12b_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');
