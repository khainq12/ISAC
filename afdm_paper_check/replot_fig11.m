clear; clc; close all;
load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig11_results.mat');
BER_ideal = floorZeros(BER_ideal); BER_est = floorZeros(BER_est);
figure('Color','w','Position',[100 100 750 560]);
semilogy(SNRd_dB, BER_ideal, '-o', 'Color', [0.05 0.6 0.05], 'LineWidth', 2.2, 'MarkerSize', 6, 'DisplayName','Ideal'); hold on;
mk={':d','--v','-.^'};
clr = {[0.85 0.1 0.1],[0.9 0.55 0.05],[0.05 0.35 0.85]};
for ip=1:numel(SNRp_list)
    semilogy(SNRd_dB, BER_est(ip,:), mk{ip}, 'Color', clr{ip}, 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', sprintf('SNRp = %d dB', SNRp_list(ip)));
end
grid on; box on; set(gca,'FontSize',12,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
legend('show','FontSize',11); ylim([1e-4 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig11_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function Bc = floorZeros(B)
Bc = B; nz = B(B>0); floorVal = min(nz(:)); Bc(B==0) = floorVal;
end
