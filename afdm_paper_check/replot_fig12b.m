clear; clc; close all;
load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12b_results.mat');
BER_ideal = floorZeros(BER_ideal); BER_est = floorZeros(BER_est);
figure('Color','w','Position',[100 100 800 600]);
mk = {'-o','-s','-^'};
clr = {[0.05 0.35 0.85],[0.85 0.35 0.05],[0.10 0.55 0.15]};
for ix=1:numel(xiList)
    semilogy(SNRd_dB, BER_est(ix,:), mk{ix}, 'Color', clr{ix}, 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', sprintf('Est, \\xi_\\nu=%d', xiList(ix))); hold on;
    semilogy(SNRd_dB, BER_ideal(ix,:), mk{ix}, 'LineStyle', '--', 'Color', clr{ix}, 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', sprintf('Ideal, \\xi_\\nu=%d', xiList(ix)));
end
grid on; box on; set(gca,'FontSize',12,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
legend('show','FontSize',9); ylim([1e-4 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig12b_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function Bc = floorZeros(B)
Bc = B; nz = B(B>0); floorVal = min(nz(:)); Bc(B==0) = floorVal;
end
