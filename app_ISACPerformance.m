function app_ISACPerformance()
% app_ISACPerformance() Interactive comm./sensing performance dashboard for the
% Doubly-Dispersive Channel Communications Package (AFDM/OTFS/OFDM, AWGN/Rayleigh).
% Left panel: system parameters, matching the visual style of app_DDChannel.mlapp.
% Right panel: BER/SER vs SNR, and range/velocity estimation RMSE vs SNR against the CRB.

fig = uifigure('Name', 'ISAC Performance Simulator', 'Position', [80 40 1320 950]);
fig.AutoResizeChildren = 'off';

mainGrid = uigridlayout(fig, [1 2]);
mainGrid.ColumnWidth = {440, '1x'};
mainGrid.RowHeight = {'1x'};
mainGrid.ColumnSpacing = 0;
mainGrid.RowSpacing = 0;
mainGrid.Padding = [0 0 0 0];

leftPanel = uipanel(mainGrid);
leftPanel.Layout.Row = 1; leftPanel.Layout.Column = 1;

rightPanel = uipanel(mainGrid);
rightPanel.Layout.Row = 1; rightPanel.Layout.Column = 2;

%% --- Left panel: controls -------------------------------------------------

uilabel(leftPanel, 'FontSize', 14, 'FontWeight', 'bold', ...
    'Position', [15 905 200 22], 'Text', 'System Parameters:');

% Waveform button group
waveformGroup = uibuttongroup(leftPanel, 'Title', 'Waveform', 'Position', [15 785 130 106]);
uitogglebutton(waveformGroup, 'Text', 'OFDM', 'Position', [11 52 100 23], 'Value', true);
uitogglebutton(waveformGroup, 'Text', 'OTFS', 'Position', [11 31 100 23]);
uitogglebutton(waveformGroup, 'Text', 'AFDM', 'Position', [11 10 100 23]);

% Channel button group
channelGroup = uibuttongroup(leftPanel, 'Title', 'Channel', 'Position', [160 785 120 106]);
uitogglebutton(channelGroup, 'Text', 'AWGN', 'Position', [10 52 100 23], 'Value', true);
uitogglebutton(channelGroup, 'Text', 'Rayleigh', 'Position', [10 20 100 23]);

% Grid size knob (matches app_DDChannel.mlapp's presets)
uilabel(leftPanel, 'HorizontalAlignment', 'center', 'Position', [300 835 100 22], 'Text', 'Grid Size');
gridSizeKnob = uiknob(leftPanel, 'discrete');
gridSizeKnob.Items = {'16', '36', '64', '100', '144', '256'};
gridSizeKnob.Value = '16';
gridSizeKnob.Position = [325 870 60 60];

% Bits per QAM symbol
uilabel(leftPanel, 'HorizontalAlignment', 'left', 'Position', [300 795 130 22], 'Text', 'QAM (bits/sym)');
bpsDropdown = uidropdown(leftPanel, 'Items', {'2 (4-QAM)', '4 (16-QAM)', '6 (64-QAM)'}, ...
    'ItemsData', [2 4 6], 'Value', 4, 'Position', [300 770 120 22]);

uilabel(leftPanel, 'FontSize', 14, 'FontWeight', 'bold', ...
    'Position', [15 735 200 22], 'Text', 'Channel Parameters:');

numPathsSpinner   = addLabeledSpinner(leftPanel, 'Number of Paths (Rayleigh)', 705, 1, [1 Inf], 1);
maxDelaySpinner   = addLabeledSpinner(leftPanel, 'Max Norm. Integer Delay',    675, 4, [1 Inf], 1);
maxDopplerSpinner = addLabeledSpinner(leftPanel, 'Max Norm. Digital Doppler',  645, 2, [0 Inf], 1);
guardWidthSpinner = addLabeledSpinner(leftPanel, 'AFDM Guard Width',           615, 2, [0 Inf], 1);
numPathsSpinner.Enable = 'off';
guardWidthSpinner.Enable = 'off';

uilabel(leftPanel, 'FontSize', 14, 'FontWeight', 'bold', ...
    'Position', [15 575 200 22], 'Text', 'System / SNR Parameters:');

carrierFreqSpinner = addLabeledSpinner(leftPanel, 'Carrier Frequency [Hz]', 545, 4e9, [1 Inf], 1e8);
bandwidthSpinner   = addLabeledSpinner(leftPanel, 'Bandwidth [Hz]',         515, 1e6, [1 Inf], 1e5);
snrMinSpinner       = addLabeledSpinner(leftPanel, 'SNR Min [dB]',   485, -5,  [-Inf Inf], 1);
snrMaxSpinner       = addLabeledSpinner(leftPanel, 'SNR Max [dB]',   455, 20,  [-Inf Inf], 1);
snrStepSpinner       = addLabeledSpinner(leftPanel, 'SNR Step [dB]', 425, 5,   [0.1 Inf], 0.5);
numTrialsSpinner     = addLabeledSpinner(leftPanel, 'Monte Carlo Trials/SNR', 395, 200, [1 Inf], 50);

runButton = uibutton(leftPanel, 'push', 'Text', 'Run Simulation', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Position', [15 345 200 34]);
statusLabel = uilabel(leftPanel, 'Position', [15 320 400 22], 'Text', 'Ready.');

%% --- Diversity order demo (exhaustive ML, small N/BPSK) ---------------------

uilabel(leftPanel, 'FontSize', 14, 'FontWeight', 'bold', ...
    'Position', [15 280 260 22], 'Text', 'Diversity Order Demo (ML, BPSK):');
uilabel(leftPanel, 'FontSize', 10, 'FontColor', [0.4 0.4 0.4], ...
    'Position', [15 258 400 20], 'Text', 'Exact ML detection reveals diversity order that ZF cannot.');

mlNDropdown = addLabeledDropdown(leftPanel, 'Grid Size N (small, exact ML)', 225, {'4','16'}, '16');
mlTrialsSpinner = addLabeledSpinner(leftPanel, 'Trials/SNR (ML, slow)', 195, 300, [10 Inf], 50);

runMLButton = uibutton(leftPanel, 'push', 'Text', 'Run Diversity Demo (ML)', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Position', [15 145 250 30]);
mlStatusLabel = uilabel(leftPanel, 'Position', [15 120 400 22], 'Text', 'Ready.');

%% --- Right panel: plots -----------------------------------------------------

plotGrid = uigridlayout(rightPanel, [3 2]);
plotGrid.RowHeight = {'1x', '1x', '1x'};
plotGrid.ColumnWidth = {'1x', '1x'};

axesBER = uiaxes(plotGrid);
axesBER.Layout.Row = 1; axesBER.Layout.Column = [1 2];
title(axesBER, 'Communication: BER/SER vs SNR');
xlabel(axesBER, 'SNR [dB]'); ylabel(axesBER, 'Error Rate');
grid(axesBER, 'on');

axesRange = uiaxes(plotGrid);
axesRange.Layout.Row = 2; axesRange.Layout.Column = 1;
title(axesRange, 'Sensing: Range RMSE vs CRB');
xlabel(axesRange, 'SNR [dB]'); ylabel(axesRange, 'Range error [m]');
grid(axesRange, 'on');

axesVel = uiaxes(plotGrid);
axesVel.Layout.Row = 2; axesVel.Layout.Column = 2;
title(axesVel, 'Sensing: Velocity RMSE vs CRB');
xlabel(axesVel, 'SNR [dB]'); ylabel(axesVel, 'Velocity error [m/s]');
grid(axesVel, 'on');

axesDiversity = uiaxes(plotGrid);
axesDiversity.Layout.Row = 3; axesDiversity.Layout.Column = [1 2];
title(axesDiversity, 'Diversity Order Demo: BER vs SNR (exact ML, BPSK, equal-strength 2-path)');
xlabel(axesDiversity, 'SNR [dB]'); ylabel(axesDiversity, 'Error Rate');
grid(axesDiversity, 'on');

%% --- Callbacks ---------------------------------------------------------

waveformGroup.SelectionChangedFcn = @(~,~) onWaveformChanged();
channelGroup.SelectionChangedFcn = @(~,~) onChannelChanged();
runButton.ButtonPushedFcn = @(~,~) onRun();
runMLButton.ButtonPushedFcn = @(~,~) onRunML();

    function onWaveformChanged()
        guardWidthSpinner.Enable = strcmp(waveformGroup.SelectedObject.Text, 'AFDM');
    end

    function onChannelChanged()
        isRayleigh = strcmp(channelGroup.SelectedObject.Text, 'Rayleigh');
        numPathsSpinner.Enable = isRayleigh;
        if ~isRayleigh
            numPathsSpinner.Value = 1;
        end
    end

    function onRun()
        runButton.Enable = 'off';
        drawnow;
        try
            waveform = waveformGroup.SelectedObject.Text;
            channelType = channelGroup.SelectedObject.Text;
            N = str2double(gridSizeKnob.Value);
            bps = bpsDropdown.Value;
            P = round(numPathsSpinner.Value);
            maxDelay = round(maxDelaySpinner.Value);
            maxDoppler = maxDopplerSpinner.Value;
            guardWidth = round(guardWidthSpinner.Value);
            carrierFreq = carrierFreqSpinner.Value;
            bandwidth = bandwidthSpinner.Value;
            samplingFreq = 2*bandwidth;
            snrdB = snrMinSpinner.Value:snrStepSpinner.Value:snrMaxSpinner.Value;
            numTrials = round(numTrialsSpinner.Value);

            statusLabel.Text = 'Running BER/SER sweep...';
            [BER, SER] = simBERSER(waveform, channelType, snrdB, N, bps, P, maxDelay, maxDoppler, ...
                guardWidth, samplingFreq, numTrials, 'ideal', N, 'genie', makeProgressFcn('BER/SER'));

            cla(axesBER, 'reset');
            semilogy(axesBER, snrdB, max(BER, eps), '-o', snrdB, max(SER, eps), '-s', 'LineWidth', 1.5);
            legend(axesBER, {'BER', 'SER'}, 'Location', 'southwest');
            title(axesBER, sprintf('Communication: BER/SER vs SNR (%s, %s)', waveform, channelType));
            xlabel(axesBER, 'SNR [dB]'); ylabel(axesBER, 'Error Rate');
            grid(axesBER, 'on');

            statusLabel.Text = 'Running sensing sweep...';
            [rangeRMSE, velRMSE, crbRangeStd, crbVelStd] = simSensing(waveform, channelType, snrdB, N, bps, ...
                P, maxDelay, maxDoppler, guardWidth, samplingFreq, carrierFreq, numTrials, makeProgressFcn('Sensing'));

            % floor at a fraction of the smallest CRB point so exact-zero RMSE (e.g. the delay
            % estimator classifying a small integer-delay grid perfectly) still renders on the log axis
            rangeFloor = max(min(crbRangeStd)*0.1, 1e-9);
            velFloor = max(min(crbVelStd)*0.1, 1e-9);

            cla(axesRange, 'reset');
            semilogy(axesRange, snrdB, max(rangeRMSE, rangeFloor), '-o', snrdB, crbRangeStd, '--', 'LineWidth', 1.5);
            legend(axesRange, {'Empirical RMSE', 'CRB (\sigma)'}, 'Location', 'northeast');
            title(axesRange, 'Sensing: Range RMSE vs CRB');
            xlabel(axesRange, 'SNR [dB]'); ylabel(axesRange, 'Range error [m]');
            grid(axesRange, 'on');

            cla(axesVel, 'reset');
            semilogy(axesVel, snrdB, max(velRMSE, velFloor), '-o', snrdB, crbVelStd, '--', 'LineWidth', 1.5);
            legend(axesVel, {'Empirical RMSE', 'CRB (\sigma)'}, 'Location', 'northeast');
            title(axesVel, 'Sensing: Velocity RMSE vs CRB');
            xlabel(axesVel, 'SNR [dB]'); ylabel(axesVel, 'Velocity error [m/s]');
            grid(axesVel, 'on');

            statusLabel.Text = 'Done.';
        catch ME
            statusLabel.Text = ['Error: ' ME.message];
        end
        runButton.Enable = 'on';
    end

    function onRunML()
        runMLButton.Enable = 'off';
        drawnow;
        try
            N = str2double(mlNDropdown.Value);
            % fixed small-scale parameters satisfying AFDM's orthogonality condition at
            % this N (independent of the main panel's larger-N spinners above)
            if N == 4
                maxDelay = 1; maxDoppler = 0; guardWidth = 0;
            else
                maxDelay = 1; maxDoppler = 1; guardWidth = 1;
            end
            snrdB = snrMinSpinner.Value:snrStepSpinner.Value:min(snrMaxSpinner.Value, 20);
            numTrials = round(mlTrialsSpinner.Value);

            waveforms = {'OFDM', 'OTFS', 'AFDM'};
            styles = {'-o', '-s', '-^'};
            BERs = cell(1,3);
            for i = 1:3
                mlStatusLabel.Text = sprintf('Running %s...', waveforms{i});
                drawnow;
                BERs{i} = simDiversityML(waveforms{i}, snrdB, N, maxDelay, maxDoppler, guardWidth, ...
                    numTrials, 'genie', makeMLProgressFcn(waveforms{i}));
            end

            cla(axesDiversity, 'reset');
            hold(axesDiversity, 'on');
            for i = 1:3
                semilogy(axesDiversity, snrdB, max(BERs{i}, eps), styles{i}, 'LineWidth', 1.5, 'DisplayName', waveforms{i});
            end
            hold(axesDiversity, 'off');
            axesDiversity.YScale = 'log';    % explicit, in case a later call resets it
            legend(axesDiversity, 'show', 'Location', 'southwest');
            title(axesDiversity, sprintf('Diversity Order Demo: BER vs SNR (exact ML, BPSK, N=%d, equal-strength 2-path)', N));
            xlabel(axesDiversity, 'SNR [dB]'); ylabel(axesDiversity, 'Error Rate');
            grid(axesDiversity, 'on');

            mlStatusLabel.Text = 'Done.';
        catch ME
            mlStatusLabel.Text = ['Error: ' ME.message];
        end
        runMLButton.Enable = 'on';
    end

    function progressFcn = makeProgressFcn(prefix)
        lastPct = -1;
        function cb(frac)
            pct = floor(frac*50)*2;
            if pct ~= lastPct
                statusLabel.Text = sprintf('%s: %d%%', prefix, pct);
                drawnow limitrate;
                lastPct = pct;
            end
        end
        progressFcn = @cb;
    end

    function progressFcn = makeMLProgressFcn(prefix)
        lastPct = -1;
        function cb(frac)
            pct = floor(frac*20)*5;
            if pct ~= lastPct
                mlStatusLabel.Text = sprintf('%s: %d%%', prefix, pct);
                drawnow limitrate;
                lastPct = pct;
            end
        end
        progressFcn = @cb;
    end

end

function spinner = addLabeledSpinner(parent, labelText, y, value, limits, step)
uilabel(parent, 'HorizontalAlignment', 'left', 'Position', [15 y 220 22], 'Text', labelText);
spinner = uispinner(parent, 'Limits', limits, 'Value', value, 'Step', step, 'Position', [245 y 175 22]);
end

function dd = addLabeledDropdown(parent, labelText, y, items, value)
uilabel(parent, 'HorizontalAlignment', 'left', 'Position', [15 y 220 22], 'Text', labelText);
dd = uidropdown(parent, 'Items', items, 'Value', value, 'Position', [245 y 175 22]);
end
