clear; clc; close all;
load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8a_results.mat');
BER = floorZeros(BER);
figure('Color','w','Position',[100 100 750 560]);
mk = {'-o','-s','-^'};
clr = {[0.75 0.15 0.05],[0 0.30 0.65],[0.10 0.55 0.15]};
for ip=1:numel(Pvals)
    semilogy(SNR_dB, BER(ip,:), mk{ip}, 'Color', clr{ip}, 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', sprintf('P = %d', Pvals(ip))); hold on;
end
grid on; box on; set(gca,'FontSize',12,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
legend('show','FontSize',11); ylim([1e-6 1]);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig8a_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function Bc = floorZeros(B)
Bc = B; nz = B(B>0); floorVal = min(nz(:)); Bc(B==0) = floorVal;
end
