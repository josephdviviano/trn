%% Mean Run Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filters data, rewinds to account for BOLD lag, and censoring head motion.
% Req: datum, mask, /SESS/params_motX.1D
% This program will eat your yummy RAM (rec: 16 GB)
% Method for computing standard deviation:
% stdev = sqrt( sum_x2/n - mean^2 )

%% %%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%input and output file names
o.dir = '/Users/jdv/data/trn/'; % experiment directory
o.sub = [5, 9, 13];     

% per experiment
o.exp = {'flicker_120_thal'};     % experiment
o.masks = {'anat_THAL_mask'};     % input mask

% data within each experiment
o.datum = {'smoothHead'};

% framewise displacement censoring
o.TRdel = 0;        % Number of TRs to remove from motParams for rewind
o.TRrest= 0;        % Number of TRs to remove from motParams for rest
o.radius = 50;      % Radius of head (mm)
o.threshP = 0.25;   % Framewise displacement (mm) for TR censor (0.25 - 0.4)

% bandpass filter design
o.filter = 1; 	        % 1 = bandpass (2Xpass kaiser, beta=2.5 FIR) 0 = none
o.TR = 1.25;            % length of TRs in seconds
o.bStart = 0.01904761904761904761904761; % lo-end of bandpass
o.bStop = 0.25;                           % hi-end of bpass

% BOLD-lag correction
o.rewind = 4;       % # of TRs to rewind to compensate for BOLD lag (~ 5s)

% cycles to average over
o.cycles = [7];

% statistics
o.takeSERR = 0;    % TR-wise SEM (BROKEN), add STEV(x) / SQRT(n);
o.renormalize = 0  % variance normalize the means

% misc 
o.print = 1;

%% Compute Framewise Displacement, Identify TRs to Censor, Plot Changes %%%%% %%
for sub = o.sub;
for ex = 1:length(o.exp);
for dat = 1:length(o.datum);

directory = [o.dir 's' int2str(sub) '/' o.exp{ex}];

numSess = length(dir([directory '/SESS*']));

aDisp = [];
runCount = 0;

if o.print == 1;
    figure('Position', [100 900 1200, 900]) ;
    subplot(3, 2, [1 2]);
    cmap.elements = bone(numSess*20);
end

numTotRuns = 0;

TR1 = o.TRdel;
if TR1 == 0;
    TR1 = 1;
end

for sess = 1:numSess;
    
    % determine the number of runs for each session
    dirRuns = dir([directory '/SESS' int2str(sess) '/raw/func*']);
    numRuns = length(dirRuns);
    numTotRuns = [numTotRuns, numRuns]; 
    
    for run = 1:numRuns;
        
        % convert deg rotation to mm, find sum(abs(deriv(motionParams)))
        runCount = runCount + 1;
        motParams = dlmread([directory '/SESS' int2str(sess) ...
                            '/params_mot' int2str(run) '.1D']);
        motParams(:,1) = degtorad(motParams(:,1))*o.radius;
        motParams(:,2) = degtorad(motParams(:,2))*o.radius;
        motParams(:,3) = degtorad(motParams(:,3))*o.radius;
        motParams = abs(diff(motParams,1,1));
        motParams = sum(motParams,2);
        motParams = motParams(TR1:length(motParams) - o.TRrest);
        aDisp(:, runCount) = motParams;
        
        if o.print == 1;
            plot(1:length(aDisp(:, runCount)), aDisp(:, runCount), ...
                        'color', cmap.elements(runCount, :), ...
                        'LineWidth', 1);
            hold all
        end
    end
end

numTotRuns = sum(numTotRuns);

% in the odd case that you aren't deleting any TRs at the beginning, insert
% a dummy motion time point for the 1st to to ensure proper alignment
if o.TRdel == 0;
    aDisp(2:length(aDisp)+1, :) = aDisp;
    aDisp(1, :) = 0;
end 

numTR = length(aDisp);

% create TR masks based on threshold
censorIDX = find(aDisp > o.threshP);
censorIDX = union(union(union(censorIDX, censorIDX-1), ...
                              censorIDX+1), censorIDX+2);
keepIDX = find(aDisp <= o.threshP);
keepIDX = setdiff(keepIDX, censorIDX);

% optionally plot the results
if o.print == 1;
cmap.axis = [0.3 0.3 0.3];
line([0 numTR], [o.threshP o.threshP], 'color', 'black', 'LineWidth', 2);
set(gcf,'color','w');
set(gca,'FontSize', 10, 'xColor', cmap.axis, 'yColor', cmap.axis);
set(gca,'XTick', [0:24:numTR]);
hold off
ylabel('Displacement (mm)', 'FontSize', 10, 'Color', cmap.axis);
xlabel('Time (TRs)', 'FontSize', 10, 'Color', cmap.axis);
xlim([1 numTR]);
title(['Removed TRs with movement > ' sprintf('%0.5g', o.threshP) ' mm'], ...
       'FontSize', 12, 'Color', cmap.axis);

aDispPlot = zeros(size(aDisp, 1), size(aDisp, 2));
aDispPlot(censorIDX) = 1;
dims = size(aDispPlot);
aDispCycPlot = reshape(aDispPlot, dims(1)/o.cycles(1), dims(2)*o.cycles(1));

subplot(3, 2, 3);
colormap(flipud(bone));
imagesc(aDispPlot');
set(gcf,'color','w');
set(gca,'FontSize', 10, 'xColor', cmap.axis, 'yColor', cmap.axis);
set(gca,'XTick', [0:24:numTR], 'YTick', [1:numTotRuns]);
box off
ylabel('Run Number', 'FontSize', 10, 'Color', cmap.axis);
xlabel('Time (TRs)', 'FontSize', 10, 'Color', cmap.axis);
ylim([0.5 numTotRuns+0.5]);
xlim([0 numTR]);
title('Raster of removed TRs per run', 'FontSize', 10, 'Color', cmap.axis);

subplot(3, 2, 4);
colormap(flipud(bone));
imagesc(aDispCycPlot');
set(gcf,'color','w');
set(gca,'FontSize', 10, 'xColor', cmap.axis, 'yColor', cmap.axis);
yTicks = [1:numTotRuns:numTotRuns*o.cycles(1)+1];
yTicks(end) = yTicks(end) - 1;
set(gca,'XTick', [0:2:numTR/o.cycles(1)], 'YTick', yTicks);
box off
ylabel('Cycle Number', 'FontSize', 10, 'Color', cmap.axis);
xlabel('Time (TRs)', 'FontSize', 10, 'Color', cmap.axis);
ylim([0.5 numTotRuns*o.cycles(1)+0.5]);
xlim([0 numTR/o.cycles(1)]);
title('Raster of removed TRs per cycle', 'FontSize', 10, 'Color', cmap.axis);

subplot(3, 2, 5);
colormap(flipud(bone));
stairs(1:length(aDisp), sum(aDisp <= o.threshP, 2), ...
      'Color', 'black', 'LineWidth', 2);
set(gcf,'color','w');
set(gca,'FontSize', 10, 'xColor', cmap.axis, 'yColor', cmap.axis);
set(gca,'XTick', [0:24:numTR], 'YTick',[min(sum(aDisp <= o.threshP, 2)):numTotRuns]);
box off
xlabel('Time (TRs)', 'FontSize', 10, 'Color', cmap.axis);
ylabel('Number of Averages Per TR', 'FontSize', 10, 'Color', cmap.axis);
xlim([0 numTR]);
ylim([min(sum(aDisp <= o.threshP, 2))-1 numTotRuns]);
title('Number of averages per TR over run', 'FontSize', 10, 'Color', cmap.axis);

subplot(3, 2, 6);
colormap(flipud(bone));
stairs(1:length(aDisp)/o.cycles(1), sum(aDispCycPlot < 1, 2), ...
      'Color', 'black', 'LineWidth', 2);
set(gcf,'Color','w');
set(gca,'FontSize', 10, 'xColor', cmap.axis, 'yColor', cmap.axis);
set(gca,'XTick', [0:2:numTR/o.cycles(1)], 'YTick', ...
        [min(sum(aDispCycPlot < 1, 2)):max(sum(aDispCycPlot < 1, 2))]);
box off
xlabel('Time (TRs)', 'FontSize', 10, 'Color', cmap.axis);
ylabel('Number of Averages Per TR', 'FontSize', 10, 'Color', cmap.axis);
xlim([0 numTR/o.cycles(1)]);
ylim([min(sum(aDispCycPlot < 1, 2))-1 max(sum(aDispCycPlot < 1, 2))]);
title('Number of averages per TR over cycle', 'FontSize', 10, 'Color', cmap.axis);

saveas(gcf, [directory '/fig_mot.fig'], 'fig');
saveas(gcf, [directory '/fig_mot.eps'], 'epsc');
clf
end

% breakerSwitch = 1;
% while breakerSwitch == 1;
%     continueText = input('Review the plot. Desire to continue? Y/N [N]: ', 's');
%     if isempty(continueText)
%         continueText = 'N';
%     end
%     if continueText == 'N';
%         error('pre_censorMeans cancelled by user')
%     elseif continueText == 'Y';
%         breakerSwitch = 0;
%     else
%         clear continueText
%     end 
% end

disp(['Processing subject ' int2str(sub) ' ' o.datum{dat}]); 
timeStart = GetSecs;
%clf(gcf);

%% Import Data & Intersession Mask, Censor TRs, & Create Mean RUN / CYCLE
mask = load_nifti([directory '/' o.masks{ex} '.nii.gz']);
dims = size(mask.vol);
numVox = dims(1)*dims(2)*dims(3);
mask.vol = reshape(mask.vol, numVox, 1);
IDXmask = find(mask.vol > 0);

% for all runs, load masked voxels into mean and sum of squares data array
aDatMean = zeros(numVox, numTR);     % Stores EPI data mean

if o.takeSERR == 1;
    aDatSD = zeros(numVox, numTR); % Stores EPI data sum of squares
end

arrayTRs = zeros(1, numTR);              % Stores # of runs in each TR
runCount = 0;                             % Run # across sessions

% compute the filter parameters
ny = 1/o.TR/2;
filtOrd = (numTR/3)-1;
filt    = fir1(filtOrd, [o.bStart/ny, o.bStop/ny], kaiser(filtOrd+1, 2.5));
                     
for sess = 1:numSess;
    
    % determine the number of runs for each session
    numRuns = length(dir([directory '/SESS' int2str(sess) '/raw/func*']));
    
    for run = 1:numRuns;
        
        % load each masked functional and add to array, recording censored TRs
        runCount = runCount + 1;
        nii = load_nifti([directory '/SESS' int2str(sess) ... 
                         '/func_' o.datum{dat} int2str(run) '.nii.gz']);
        nii.vol = reshape(nii.vol, numVox, numTR);
        
        censorIDX = find(aDisp(:, runCount) >  o.threshP);
        censorIDX = union(union(union(censorIDX,    censorIDX-1), ...
                                      censorIDX+1), censorIDX+2);
        censorIDX(censorIDX == 0) = [];                       
        keepIDX   = find(aDisp(:, runCount) <= o.threshP);
        keepIDX = setdiff(keepIDX, censorIDX);
        
        % update the TR-wise count
        arrayTRs(1, keepIDX) = arrayTRs(1, keepIDX) + 1;
        
        % bandpass
        if o.filter == 1;
            nii.vol(IDXmask, :) = filtfilt(filt, 1, nii.vol(IDXmask, :)')';
        end
        
        % load sum data into mean array
        aDatMean(IDXmask, keepIDX) = ...
             nii.vol(IDXmask, keepIDX) + aDatMean(IDXmask, keepIDX);
        
        % load sum of square data into standard error array
        if o.takeSERR == 1;
            aDatSD(IDXmask, keepIDX) = ...
             nii.vol(IDXmask, keepIDX).^2 + aDatSD(IDXmask, keepIDX);
        end
        
        % print an update
        disp(fprintf('Run %d processed. \n', run));
        
    end
end

arrayTRs = repmat(arrayTRs, numVox, 1);

txt = fprintf('All data arrays loaded. \n');
disp(txt);

%% Create Mean Run & Cycles, ignoring Censored TRs
if isempty(o.cycles) == 0;
    
    for cyc = o.cycles;
    
    % reshape data and sum over data to find mean cycle
    arrayCycMean = aDatMean;
    arrayCycMean = reshape(arrayCycMean, numVox, numTR/cyc, cyc);
    arrayCycMean = sum(arrayCycMean, 3);
    
    arrayCTR = arrayTRs;
    arrayCTR = reshape(arrayCTR, numVox, numTR/cyc, cyc);
    arrayCTR = sum(arrayCTR, 3);
    
    arrayCycMean = arrayCycMean ./ arrayCTR;
    
    % variance normalize
    if o.renormalize == 1;
        arrayCycVariance = std(arrayCycMean, [], 2);
        arrayCycVariance = repmat(arrayCycVariance, [1 numTR/cyc]);
        arrayCycMean = arrayCycMean ./ arrayCycVariance;
    end
	
	% rewind data (compensate for BOLD lag)
	arrayCycMean = circshift(arrayCycMean, [0 -o.rewind]);
    
    % save
    nii.vol = arrayCycMean;
    nii.vol = reshape(nii.vol, dims(1), dims(2), dims(3), numTR/cyc);
    nii.dim(5) = numTR/cyc;
    save_nifti(nii, [directory '/mean_CYCLE_' o.datum{dat} '.nii.gz']);
    
    % find TR-wise standard error of the mean, save
    if o.takeSERR == 1;
        
        arrayCycSD = aDatSD;
        arrayCycSD = reshape(arrayCycSD, numVox, numTR/cyc, cyc);
        arrayCycSD = sum(arrayCycSD, 3);
        arrayCycSD = sqrt(arrayCycSD ./ arrayCTR - arrayCycMean.^2);
        arrayCycSD = arrayCycTSRD ./ sqrt(arrayCTR);
        
        niiSD = nii;
        niiSD.vol = arrayCycSD;
        niiSD.vol = reshape(niiSD.vol, dims(1), dims(2), dims(3), numTR/cyc);
        niiSD.dim(5) = numTR/cyc;
        save_nifti(niiSD, [directory '/sdev_CYCLE_' o.datum{dat} '.nii.gz']);
        
    end

    txt = fprintf('Mean Cycle %d Complete. \n', cyc); 
    disp(txt);
    
    end
end

%% take mean run
aDatMean = aDatMean ./ arrayTRs; 

%% variance normalize
if o.renormalize == 1;
    
    aDatVariance = std(aDatMean, [], 2);
    aDatVariance = repmat(aDatVariance, [1 numTR]);
    aDatMean = aDatMean ./ aDatVariance;
end

% rewind data (compensate for BOLD lag)
aDatMean = circshift(aDatMean, [0 -o.rewind]);

%% save mean and standard deviation of the full run
nii.vol = aDatMean;
nii.vol = reshape(nii.vol, dims(1), dims(2), dims(3), numTR);
save_nifti(nii, [directory '/mean_RUN_' o.datum{dat} '.nii.gz']);

%% find TR-wise standard deviation
if o.takeSERR == 1;
    
    aDatSD = sqrt(aDatSD ./ arrayTRs - aDatMean.^2);
    
    niiSD = nii;
    niiSD.vol = aDatSD;
    niiSD.vol = reshape(niiSD.vol, dims(1), dims(2), dims(3), numTR);
    save_nifti(niiSD, [directory '/sdev_RUN_' o.datum{dat} '.nii.gz']);
end

% print exit text
disp(fprintf('Mean Run Complete. \n'));
disp(fprintf('Processing time %0.5g minutes \n', (GetSecs - timeStart)/60)); 
disp(['Run ' o.datum{dat} ' complete.']);

end
end
end

%% Jun 30th 2013, JDV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%