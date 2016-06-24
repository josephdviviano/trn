%% Periodic Signal + Noise Simulation 

clear all; clc;

%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% All strengths are relative to opt.cycle1Strength = 1

opt.TR             = 1; % TR in seconds

% HRF period (rotating hemifield - carrier wave)
opt.cycle1Number   = 10;   % number of cycles per scan 1
opt.cycle1Length   = 18;   % length of cycle 1

% modulation period (flicker)
opt.cycle2Number   = 6;    % number of cycles per scan 2
opt.cycle2Length   = 30;   % length of cycle 2
opt.cycle2Amp      = 0.2;  % amp of cycle 2 in multiples of cycle1Strength = 1

opt.runNumber      = 40; % number of runs collected
opt.iterations     = 15;  % number of simulations (multiverse)

opt.noiseAmp       = 10;   % amp of noise in multiples of cycle1Strength = 1
        
% Smoothing Options
opt.smoothMethod   = 3;    % 0 = Nothing, 1 = MA, 2 = LOWESS, 3 = Cubic Spline
opt.smoothMA       = 5;    % 1) number of TRs for averaging
opt.smoothLO       = 5;    % 2) number of TRs for linear fitting
opt.smoothSP       = 0;    % 3) fit error (0 < x < 1, 0 for default)

opt.printName = ['sim1.jpg']

%% check settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if rem(opt.cycle1Number*opt.cycle1Length/opt.TR, 1) ~= 0;
    error('TR does not divide into run');
end

if rem(opt.cycle2Number*opt.cycle2Length/opt.TR, 1) ~= 0;
    error('TR does not divide into run');
end

if rem(opt.cycle1Number*opt.cycle1Length/opt.cycle2Length, 1) ~= 0;
    error('Cycles are not integral');
end

%if rem(opt.cycle1Length/opt.TR, 1) ~= 0;
%    error('TR does not divide into cycle 1')
%end

if rem(opt.cycle2Length/opt.TR, 1) ~= 0;
    error('TR does not divide into cycle 2')
end

harm.cycle1 = [1:2] * opt.cycle1Number;
harm.cycle2 = [1:2] * opt.cycle2Number;
harm.cycle1 = perms(harm.cycle1);
harm.dim = size(harm.cycle1);
harm.cycle2 = repmat(harm.cycle2, [harm.dim(1) 1]);
harm.test1 = harm.cycle1 ./ harm.cycle2;
harm.test2 = harm.cycle2 ./ harm.cycle1;
harm.test0 = ones(harm.dim(1), harm.dim(2));
harm.testIDX1 = find(rem(harm.test1, harm.test0) == 0);
harm.testIDX2 = find(rem(harm.test2, harm.test0) == 0);
harm.out = cat(1, harm.testIDX1, harm.testIDX2);

if isempty(harm.out) == 0;
    error('Cycles are harmonic');
end

%% initialize some things %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% colormap scaling
cmap.back = bone(round(opt.iterations*3/2));
cmap.frnt1 = autumn(round(opt.iterations*2));
cmap.frnt2 = summer(round(opt.iterations*2));
cmap.axis = [0.3, 0.3, 0.3];

% resample everything to TRs
opt.cycle1Length = opt.cycle1Length/opt.TR;
opt.cycle2Length = opt.cycle2Length/opt.TR;

% init some arrays
rArray = zeros(opt.iterations, 2);
datMeanArray = zeros(opt.iterations, opt.cycle1Length*opt.cycle1Number);

% heamodynamic response fxn = SPM double gamma w/ standard settings
[HRF, p] = spm_hrf(opt.TR);

figure('Position', [0 1200 1200 1200]);
subplot(4,2,[1 2]);
title(sprintf(...
     ['Phase-encoded simulation: timeseries, spectra, and r values'          ... 
      '\n 1st cycle: %g cycles x %g TRs, 2nd cycle: %g cycles x %g TRs'      ...
      '\n TR of %g s, with %g runs at %gx noise and %g smooth method'],      ...
      opt.cycle1Number, opt.cycle1Length, opt.cycle2Number, opt.cycle2Length,...
      opt.TR, opt.runNumber, opt.noiseAmp, opt.smoothMethod),                ...
      'Color', cmap.axis, 'FontSize', 14);

%% Run Simulation n times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:opt.iterations;
    
    %% Generate Data and Analyze
    % generate multiple noisy voxels and take the demeaned, detrended, average
    for run = 1:opt.runNumber;    
        timeRun = 1:opt.cycle1Length * opt.cycle1Number;
        timeCyc = 1:opt.cycle2Length;
        
        % build HRF time series (kludge to allow for odd periods)
        
        flipper = 0;
        counter = 1;
        signal1 = [];
        
        while counter <= opt.cycle1Number;
            
            if flipper == 0;
            
                signalNewCyc = [ones(1, floor(opt.cycle1Length/2)) ...
                              zeros(1, ceil(opt.cycle1Length/2))];
                signal1 = [signal1, signalNewCyc];     
                counter = counter + 1;
                flipper = 1;
            
            end
            
            if flipper == 1;
                signalNewCyc = [ones(1, ceil(opt.cycle1Length/2)) ...
                               zeros(1, floor(opt.cycle1Length/2))];
                signal1 = [signal1, signalNewCyc];     
                counter = counter + 1;
                flipper = 0;
            
            end
        end
                
        signal1 = conv(signal1, HRF, 'same');
        %signal1 = sin(timeRun*2*pi / opt.cycle1Length);
        
        % build modulatory time series and noise
        signal2 = sin(timeRun*2*pi / opt.cycle2Length) * opt.cycle2Amp;
        noise = unifrnd(-1, 1, [1 opt.cycle1Length*opt.cycle1Number]) * opt.noiseAmp;
        
        % now create the final time series, demean, and detrend
        TS = signal1+signal2+noise;
        TSmean = mean(TS);
        TS = TS - TSmean;
        TS = detrend(TS);
        
        % smooth the data
        if opt.smoothMethod == 1;
            TSlength = size(TS);
            TS = repmat(TS, [1 3]);
            TS = smooth(TS, opt.smoothMA, 'moving');
            TS = TS(TSlength(2)+1:TSlength(2)*2);
        end
        
        if opt.smoothMethod == 2;
            TS = smooth(TS, opt.smoothLO, 'rlowess');
        end
        
        if opt.smoothMethod == 3;
            if opt.smoothSP == 0;
                opt.smoothSP = 1 / (1 + (opt.TR^3) / 6);
            end
            TS = csaps(timeRun, TS, opt.smoothSP, timeRun);
        end
        
        % load time series into data array
        datRun(run,:) = TS;
        TS = reshape(TS, opt.cycle2Length, opt.cycle2Number);
        TS = mean(TS, 2);
        datCyc(run,:) = TS;
        
    end
    
    datMeanRun = mean(datRun, 1);
    datMeanCyc = mean(datCyc, 1);
    datMeanRunArray(iter, :) = datMeanRun;
    datMeanCycArray(iter, :) = datMeanCyc;
    
    % fft analysis of data
    datFT = fft(datMeanRun);
    scaledAmp = abs(datFT(:,1:length(timeRun)/2));
    tsPower = sqrt(sum((scaledAmp(:, 1:length(timeRun)/2).^2)));

    r1 = scaledAmp(:, opt.cycle1Number+1) ./tsPower;
    r1(isnan(r1)) = 0;

    r2 = scaledAmp(:, opt.cycle2Number+1) ./tsPower;
    r2(isnan(r2)) = 0;

    % store found correlations in an array (:,1) = r1, (:,2) = r2
    rArray(iter, 1) = r1;
    rArray(iter, 2) = r2;
    
    %% plot results
    subplot(4,2,[1 2]);
    hold all;
    plot(timeRun, datMeanRun, 'LineWidth', 2, 'color', cmap.back(iter, :));
    set(gcf,'color','w');
    set(gca,'FontSize', 12);
    set(gca,'xColor', cmap.axis);
    set(gca,'yColor', cmap.axis);
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlim([1 length(timeRun)]);
    ylim([-5 5]);
    xlabel('Time (s)', 'Color', cmap.axis);
    ylabel('% signal change', 'Color', cmap.axis);

    subplot(4,2,[3 4]);
    hold all;
    plot(timeCyc, datMeanCyc, 'LineWidth', 2, 'color', cmap.back(iter, :));
    set(gcf,'color','w');
    set(gca,'FontSize', 12);
    set(gca,'xColor', cmap.axis);
    set(gca,'yColor', cmap.axis);
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlim([1 length(timeCyc)]);
    ylim([-opt.cycle2Amp*5 opt.cycle2Amp*5]);
    xlabel('Time (s)', 'Color', cmap.axis);
    ylabel('% signal change', 'Color', cmap.axis);
    
    subplot(4,2,[5 6]);
    hold all;
    stem(1:length(timeRun)/2, scaledAmp, ...
        'LineWidth', 1.5, 'color', cmap.back(iter, :));
    stem(opt.cycle1Number+1, scaledAmp(opt.cycle1Number+1), ...
        'LineWidth', 2, 'color', cmap.frnt1(iter, :));
    stem(opt.cycle2Number+1, scaledAmp(opt.cycle2Number+1), ...
        'LineWidth', 2, 'color', cmap.frnt2(iter, :));
    set(gcf,'color','w');
    set(gca,'FontSize', 12);
    set(gca,'xColor', cmap.axis);
    set(gca,'yColor', cmap.axis);
    xlim([1 length(timeRun)/2]);
    xlabel('Cycles per run', 'Color', cmap.axis);
    ylabel('Power', 'Color', cmap.axis);

end

%% Finish Plots and Print Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
datMeanRunGrand = mean(datMeanRunArray);
subplot(4,2,[1 2]);
hold all;
plot(timeRun, datMeanRunGrand, 'LineWidth', 2, 'color', 'red');

datMeanCycGrand = mean(datMeanCycArray);
subplot(4,2,[3 4]);
hold all;
plot(timeCyc, datMeanCycGrand, 'LineWidth', 2, 'color', 'green');

subplot(4,2,7);
hold all;
title('r values at 1st cycle frequency ', 'Color', cmap.axis,'FontSize', 12);
hist(rArray(:,1), 20);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', 'r', 'EdgeColor', 'w')
set(gcf,'color','w');
set(gca,'FontSize', 12);
set(gca,'xColor', cmap.axis);  
set(gca,'yColor', cmap.axis);
set(gca, 'xlim', [0 1]);
xlabel('correlation (r)', 'Color', cmap.axis);
ylabel('Count', 'Color', cmap.axis);

subplot(4,2,8);
hold all;
title('r values at 2nd cycle frequency ', 'Color', cmap.axis, 'FontSize', 12);
hist(rArray(:,2), 20);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', 'g', 'EdgeColor', 'w')
set(gcf,'color','w');
set(gca,'FontSize', 12);
set(gca,'xColor', cmap.axis);
set(gca,'yColor', cmap.axis);
set(gca, 'xlim', [0 1]);
xlabel('correlation (r)', 'Color', cmap.axis);


saveas(gcf, opt.printName, 'fig');
outName = eval('str2mat(opt.printName)');
outCmd = strcat(eval('str2mat(opt.printName)'), [' -native -nocrop -a1']);
export_fig(outName, '-native', '-nocrop', '-a1');

%% JDV Feb 6th 2013 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%