%% create stat mask                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Input statistical map                                                        %
%% %%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
o.dir = '/srv/data/MSC_JNEUROSCI/';             % experiment directory
o.sub = {'01','02','05','07','09','11'};       

% per experiment
%o.exp = {'flicker_120_thal'};     % experiment
o.datum = {'mean_CYCLE_bpass_smoothHead'};
o.mask = 'mask_stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHead';

o.stats = 'stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHead';
o.vol = 1; %% flicker phase

o.labels = 'mask_LGN_TRN_resamp';

o.TRNs = [133, 27, 31; ...
          67,  17, 32; ...
          62,  29, 48; ...
          73,  20, 36; ... 
          123,  23, 40; ...
          131, 17, 39];

o.test = [2, 1; ...
          1, 2; ...
          1, 2; ...
          1, 2; ...
          2, 1; ...
          2, 1]; %ipsi, contra

o.mean = 1; % 1 = plot the mean, not the max anticorr

o.prefix = 'LGNTRN'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o.TRNs = o.TRNs + 1;
labels = {'unknown', 'R LGN', 'L LGN', 'R TRN', 'L TRN', 'R PUL', 'L PUL'};
figure;

R_VEC = [];
GROUP_MEAN=[];
LGN_MEAN=[];
TRN_MEAN=[];
iter = 1;
for s = o.sub;
    
%directory = [o.dir 's' int2str(s) '/' o.exp{1}];
directory = [o.dir 's' s{1}];
for d = 1:length(o.datum);

    %% Import Mask & Data, Compute Dimensions
    data = load_nifti([directory '/' o.datum{d} '.nii.gz']);
    mask = load_nifti([directory '/' o.mask '.nii.gz']);
    stat = load_nifti([directory '/' o.stats '.nii.gz']);
    labs = load_nifti([directory '/' o.labels '.nii.gz']);

    dims = size(data.vol);

    nVox = dims(1)*dims(2)*dims(3);
    nTRs = dims(4);

    seed = data.vol;
    seed = zeros(dims(1), dims(2), dims(3));
    seed(o.TRNs(iter, 1), o.TRNs(iter, 2), o.TRNs(iter, 3)) = 1;

    %OUT.dim(5) = 10;

    data = reshape(data.vol, nVox, nTRs);
    mask = reshape(mask.vol, nVox, 1);
    stat = stat.vol(:, :, :, o.vol);
    stat = reshape(stat, nVox, 1);
    labs = reshape(labs.vol, nVox, 1);
    seed = reshape(seed, nVox, 1);

    idx = find(seed == 1);
    idx_lgns = find(labs == o.test(iter, 1));
    
    if o.mean == 0;
      rVec = corr(data(find(seed == 1), :)', data(idx_lgns, :)');
      idx_anti = find(rVec == min(rVec));
    end

    subplot(3, 2, iter);

    if o.mean == 0;
      plot(1:nTRs, data(idx_lgns(idx_anti), :), 'color', 'r', 'linewidth', 2); 
    end
    
    if o.mean == 1;
      ts = data(idx_lgns, :) ./ repmat(max(abs(data(idx_lgns, :)), [], 2), [1 nTRs]);
      nan_idx = find(isnan(ts(:, 1)) == 0);
      ts = ts(nan_idx, :);
      ts_std = std(ts);
      ts = mean(ts);

      %ts = resample(ts, 301, 24);
      %ts = smooth(ts, 10);
      plot(1:24, ts, 'color', 'r', 'linewidth', 2); hold all;
      plot(1:24, ts+ts_std, 'color', 'r', 'linewidth', 1);
      plot(1:24, ts-ts_std, 'color', 'r', 'linewidth', 1);
    end

    GROUP_MEAN = [GROUP_MEAN; data(idx_lgns, :)];
    LGN_MEAN = [LGN_MEAN; nanmean(data(idx_lgns, :))];
    TRN_MEAN = [TRN_MEAN; data(find(seed == 1), :)];

    hold all;

    %scaled data
    data(isnan(data)) = 1;
    data = data ./ repmat(max(abs(data), [], 2), [1 nTRs]);

    ts = data(find(seed == 1), :);
    ts = ts ./ repmat(max(abs(ts), [], 2), [1 nTRs]);

    plot(1:24, ts, 'color', 'black', 'linewidth', 2);

    xlim([1 24]);
    ylim([-1.1 1.1]);

    lins = linspace(1,24,11);
    ticks = [];
    for x = 1:length(lins)-1;
        ticks = [ticks ((lins(x)+lins(x+1))/2)];
    end;

    set(gca,'XTick', ticks, 'YTick', [-1:1:1]);

    set(gca, 'xTickLabel', [1, 2.5, 5, 7.5, 10, 12, 15, 20, 30, 60]);
    xlabel('Flicker (Hz)', 'FontSize', 12);
    ylabel('Normalized Response (A.U.)', 'FontSize', 12);
    
end

disp([':::Subject ' s ' complete:::']);

iter = iter + 1;

end

    % %% save outputs
    %OUT.vol = reshape(OUT.vol, dims(1), dims(2), dims(3), 10);
    %OUTstr = [directory '/stats_network_' o.prefix '_' o.datum{d} '.nii.gz'];
    %save_nifti(OUT, OUTstr);


%% Sepy 20th 2013, JDV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% stats of difference between TRN and LGN...

H = [];
P = [];

crit = 0.05/24;

for x = 1:24; 
	[h, p] = ttest2(LGN_MEAN(:,x), TRN_MEAN(:,x)); 
    
    if p > crit;
    	h = 0;
    end

    H = [H, h];
    P = [P, p];
end