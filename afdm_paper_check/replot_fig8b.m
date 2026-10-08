clear; clc; close all;
load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8b_results.mat');
BER = floorZeros(BER);
figure('Color','w','Position',[100 100 750 560]);
mk = {':d','--v','-o','-^'};
clr = {[0.85 0.1 0.1],[0.9 0.5 0.05],[0.05 0.35 0.85],[0.10 0.55 0.15]};
for w=1:4
    semilogy(SNR_dB, BER(w,:), mk{w}, 'Color', clr{w}, 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', waveNames{w}); hold on;
end
grid on; box on; set(gca,'FontSize',12,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
legend('show','FontSize',11); ylim([1e-5 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8b_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function Bc = floorZeros(B)
Bc = B; nz = B(B>0); floorVal = min(nz(:)); Bc(B==0) = floorVal;
end
