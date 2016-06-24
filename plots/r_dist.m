%% Computes Distrubution of r Values Obtained via FFT Analysis %%%%%%%%%%%%%% %%
% Directory structures
opt.dir = '/media/OS/Users/jdv/viviano/Dropbox/data/';
opt.sub = [1, 2, 5, 7, 9, 11];
%opt.exp = {'flicker_120_thal'}; 

% file names (no suffix)
opt.stats = {'stats_FFT_single_multitaper_mean_RUN_allpass_smoothHead'};
opt.mask = 'anat_THAL_mask';
opt.titles ={'Retinotopy', 'Flicker'};

% Mask threshold settings
opt.statVols = [4, 14];          % r value volumes in FFT stats
opt.threshVols = [9, 19];        % r value for FFT stat thresh
opt.threshVals = [0.1, 0.1];
opt.threshCont = 9           % # of contigous voxels req.

%% Import Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100 900 1200, 900]);
set(gcf,'Units','normal')
set(gcf,'Position',[0.05 0.05 0.9 0.9])

cmap.axis = [0.3, 0.3, 0.3];
rVals = [];
binVals = [0:2:50];

figureIter = 1;

maxVec = [];

for sub = opt.sub;
%for ex = 1:length(opt.exp);

    directory = [opt.dir 's' int2str(sub) '/']; % opt.exp{ex}
    % load & reshape stats
    stat = load_nifti([directory '/' opt.stats{1} '.nii.gz']);
    mask = load_nifti([directory '/' opt.mask '.nii.gz']);
    dims = size(stat.vol);
    numVox = dims(1)*dims(2)*dims(3);

    % binarize % remove islands smaller than (opt.threshCont)
    if opt.threshCont > 0;
        iter = 1;
        for vol = opt.threshVols;
            bin = zeros(dims(1), dims(2), dims(3), 1);
            tmp = stat.vol(:, :, :, vol);
            tmp(tmp < opt.threshVal(iter)) = 0;
            tmp(tmp >= opt.threshVal(iter)) = 1;
            tmp = bwareaopen(tmp, opt.threshCont(iter));
            tmp = find(tmp > 0);
            bin(tmp) = 1;
            stat.vol(:, :, :, vol) = bin;
            iter = iter + 1;
        end
    end

    stat.vol = reshape(stat.vol, numVox, 20)';
    mask.vol = reshape(mask.vol, numVox, 1)';

    % find dual-frequency masks with thresholds
    for iter = 1:length(opt.statVols);
        
        statvol = opt.statVols(iter);
        maskvol = opt.threshVols(iter);
        thrshld = opt.threshVals(iter);
        IDXm = find(stat.vol(maskvol, :) > thrshld);
        rVec = stat.vol(statvol, IDXm); 
        n_sig = length(rVec);
        n_tot = length(mask.vol(mask.vol == 1));
        
        % plot normalized histograms showing % voxels
        subplot(length(opt.sub), length(opt.statVols), figureIter);
        binCount = hist(rVec, binVals);
        bar(binVals-0.025, binCount/sum(binCount), 1);

        maxVec = [maxVec, max(binCount / sum(binCount))];

        h = findobj(gca, 'Type', 'patch');
        set(h, 'FaceColor', 'black', 'EdgeColor', 'w');
        set(gcf,'color','w');
        set(gca,'FontSize', 12);
        set(gca,'xColor', cmap.axis);
        set(gca,'yColor', cmap.axis);
        set(gca, 'xlim', [0 50]);
        xlabel('F', 'Color', cmap.axis);
        ylabel('voxels (%)', 'Color', cmap.axis);

        if rem(figureIter, 2) == 0;
            title([opt.titles{1} ' proportion of voxels = ' num2str(n_sig/n_tot)], 'FontSize', 12, 'Color', cmap.axis);
        else;
            title([opt.titles{2} ' proportion of voxels = ' num2str(n_sig/n_tot)], 'FontSize', 12, 'Color', cmap.axis);
        end;

        box off;
        figureIter = figureIter + 1;    
    end
    
%end
%mtit(['r value distributions for subject ' sprintf('%d', sub)], ...
%       'FontSize', 14, 'Color', cmap.axis, 'xoff', 0, 'yoff', 0.05);
end

figureIter = 1;
for sub = opt.sub;
    for iter = 1:length(opt.statVols);
        subplot(length(opt.sub), length(opt.statVols), figureIter);
        set(gca, 'ylim', [0 max(maxVec)]);
        figureIter = figureIter + 1;
    end
end

%% May 17th 2013 JDV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%