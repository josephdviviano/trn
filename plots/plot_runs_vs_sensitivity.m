%% Computes Distrubution of r Values Obtained via FFT Analysis %%%%%%%%%%%%%% %%
% Directory structures
opt.dir = '/media/OS/Users/jdv/viviano/Dropbox/data/';
opt.sub = [1, 2, 5, 7, 9, 11];

opt.runs = [25, 30, 20, 30, 23, 30];

%opt.exp = {'flicker_120_thal'}; 

% file names (no suffix)
opt.stats = {'stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHead'};
opt.mask = 'mask_LGN_TRN_nopul';
opt.titles ={'Retinotopy', 'Flicker'};

% Mask threshold settings
opt.statVols = [14];          % r value volumes in FFT stats
opt.threshVols = [19];        % r value for FFT stat thresh
opt.threshVals = [0.1];
opt.threshCont = 0           % # of contigous voxels req.

%% Import Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = [];
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
            tmp(tmp < opt.threshVals(iter)) = 0;
            tmp(tmp >= opt.threshVals(iter)) = 1;
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

        n = [n n_sig/n_tot];
    end
end

figure;
scatter(opt.runs, n, [], [1 0 0], 'filled');
hold all;
[b, bint, r, rint, stats] = regress(n', [opt.runs' ones(length(opt.runs), 1)]);
x = 18:32;
y = x * b(1) + b(2);
chi = x * bint(1,1) + bint(2,1);
clo = x * bint(1,2) + bint(2,2);

plot(x,y, 'linewidth', 2, 'color', 'black');
plot(x,chi, 'linestyle', '--', 'color', 'black');
plot(x,clo, 'linestyle', '--', 'color', 'black');

xlabel(['Number of Averaged Runs']);
ylabel(['Proportion of Significant Voxels (%)']);

%xlim([18 32]);
%ylim([-0.2 0.3]);

title(['R^2 = ' num2str(stats(2)) ', p = ' num2str(stats(3))]);

%% Oct 2nd 2013 JDV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%