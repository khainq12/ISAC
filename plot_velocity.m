load('tables_velocity.mat');
% BER_vel: [waveform x velocity x snr], waveforms={'OFDM','OTFS','AFDM'}
outdir = '';
styles = {'-o','-s','-^'};
for k = 1:length(snrdB)
    fig = figure('Position',[100 100 900 400],'Visible','off');
    hold on;
    for i = 1:3
        y = max(squeeze(BER_vel(i,:,k)), 1e-5);
        plot(velocities, y, styles{i}, 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', waveforms{i});
    end
    hold off;
    set(gca,'YScale','log');
    grid on; box on;
    xlabel('Tốc độ mục tiêu [km/h]'); ylabel('BER (thu ML)');
    title(sprintf('BER theo tốc độ mục tiêu, SNR = %d dB (fc = 28 GHz)', snrdB(k)));
    legend('Location','best');
    xlim([min(velocities) max(velocities)]);
    saveas(fig, sprintf('%svel_ber_snr%d.png', outdir, snrdB(k)));
end
fprintf('Plots saved.\n');
