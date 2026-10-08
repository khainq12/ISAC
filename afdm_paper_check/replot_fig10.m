clear; clc; close all;
load('D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig10_results.mat');
detNames = {'MMSE','MRC-DFE','low-complexity MMSE'};
waveNames = {'AFDM','OFDM'};
mk = {'-o','-s','-d'};
clr = {[0 0.30 0.65],[0.75 0.20 0.05]};
LW = 2; MS = 5;
figure('Color','w','Position',[100 100 1300 520]);
subplot(1,2,1);
BERa = floorZeros(results.integer.BER); SNRa = results.integer.SNR_dB;
for w=1:2
    for d=1:3
        semilogy(SNRa, squeeze(BERa(d,w,:)), mk{d}, 'Color', clr{w}, 'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', sprintf('%s, %s', waveNames{w}, detNames{d})); hold on;
    end
end
grid on; box on; set(gca,'FontSize',11,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
legend('show','FontSize',8);
subplot(1,2,2);
BERb = floorZeros(results.fractional.BER); SNRb = results.fractional.SNR_dB;
for w=1:2
    for d=1:3
        semilogy(SNRb, squeeze(BERb(d,w,:)), mk{d}, 'Color', clr{w}, 'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', sprintf('%s, %s', waveNames{w}, detNames{d})); hold on;
    end
end
grid on; box on; set(gca,'FontSize',11,'LineWidth',1);
xlabel('SNR in dB','FontWeight','bold'); ylabel('BER','FontWeight','bold');
legend('show','FontSize',8);
exportgraphics(gcf, 'D:\doubly-dispersive-channel-simulation-main\afdm_paper_check\fig10_reproduction.png', 'Resolution', 180);
fprintf('=== DONE ===\n');

function Bc = floorZeros(B)
Bc = B; nz = B(B>0); floorVal = min(nz(:)); Bc(B==0) = floorVal;
end
