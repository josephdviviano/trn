%% create stat mask                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Input statistical map                                                        %
%% %%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
o.dir = '/home/jdv/Dropbox/data/';             % experiment directory
o.sub = [1,2,5,7,9,11];       

% per experiment
o.datum = 'mean_CYCLE_bpass_smoothHeadLGNTRN';
o.mask = 'mask_stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHeadLGNTRN';
o.affinity = 'stats_correlation_affinity_LGNTRN';
o.stats = 'stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHeadLGNTRN';
o.vol = 1; %% flicker phase

o.labels = 'mask_LGN_TRN_resamp_split';

o.prefix = 'LGNTRN'

o.freq = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

TRN_LO_GROUP = [];
TRN_HI_GROUP = [];
LGN_GROUP = [];

iter = 1;
for s = o.sub;
    
    %directory = [o.dir 's' int2str(s) '/' o.exp{1}];
    directory = [o.dir 's' int2str(s)];
    
    %% Import Mask & Data, Compute Dimensions
    mask = load_nifti([directory '/' o.mask '.nii.gz']);
    affi = load_nifti([directory '/' o.affinity '.nii.gz']);
    stat = load_nifti([directory '/' o.stats '.nii.gz']);
    labs = load_nifti([directory '/' o.labels '.nii.gz']);

    dims = size(data.vol);

    nVox = dims(1)*dims(2)*dims(3);
    nTRs = dims(4);

    data = reshape(data.vol, nVox, nTRs);
    mask = reshape(mask.vol, nVox, 1);
    stat = stat.vol(:, :, :, o.vol);
    stat = reshape(stat, nVox, 1);
    labs = reshape(labs.vol, nVox, 1);
    rval = affi.vol(:, :, :, 2);
    affi = affi.vol(:, :, :, 1);
    affi = reshape(affi, nVox, 1);
    rval = reshape(rval, nVox, 1);
    
    TRN_LO = [];
    TRN_HI = [];
    LGN = [];

    freqBins = linspace(0, 2*pi, 6); 


    idx_lo = intersect(union(find(labs == 3), find(labs == 4)), find(mask > 0));;

    % labs 3/7 = right, 4/8 = left
    % affi 1 = right, 2 = left

    for x = 1:length(unique(mask))-1;
        idx_mask = find(mask == x);
        for y = 1:length(idx_mask);              
            % LGN
            if labs(idx_mask(y)) == 1;
                LGN = [LGN; data(idx_mask(y), :)];
                LGN_GROUP = [LGN_GROUP; data(idx_mask(y), :)];
            end;
            if labs(idx_mask(y)) == 2;
                LGN = [LGN; data(idx_mask(y), :)];
                LGN_GROUP = [LGN_GROUP; data(idx_mask(y), :)];
            end;
            % TRN
            if stat(idx_mask(y)) > freqBins(o.freq) && stat(idx_mask(y)) <= freqBins(o.freq+1);
                if rval(idx_mask(y)) > percentile(rval(idx_lo), 30);    
                    if labs(idx_mask(y)) == 3 && affi(idx_mask(y)) == 2;
                        TRN_LO = [TRN_LO; data(idx_mask(y), :)];
                        TRN_LO_GROUP = [TRN_LO_GROUP; data(idx_mask(y), :)];
                    end;
                    if labs(idx_mask(y)) == 4 && affi(idx_mask(y)) == 1;
                        TRN_LO = [TRN_LO; data(idx_mask(y), :)];
                        TRN_LO_GROUP = [TRN_LO_GROUP; data(idx_mask(y), :)];
                    end;
                end
            end
        end
    end

    %nan_idx = find(isnan(LGN(:, 1)) == 0);
    %LGN = LGN(nan_idx, :);
    %nan_idx = find(isnan(TRN_LO(:, 1)) == 0);
    %TRN_LO = TRN_LO(nan_idx, :);
    %nan_idx = find(isnan(TRN_HI(:, 1)) == 0);
    %TRN_HI = TRN_HI(nan_idx, :);

    LGN = LGN ./ repmat(max(abs(LGN)')', [1 24]);
    TRN_LO = TRN_LO ./ repmat(max(abs(TRN_LO)')', [1 24]);
    %TRN_HI = TRN_HI ./ repmat(max(abs(TRN_HI)')', [1 24]);

    %LGN = LGN ./ repmat(std(LGN')', [1 24]);
    %TRN_LO = TRN_LO ./ repmat(std(TRN_LO')', [1 24]);

    LGN_MEAN = mean(LGN);
    TRN_LO_MEAN = mean(TRN_LO);
    %TRN_HI_MEAN = mean(TRN_HI);

    LGN_SEM = std(LGN) ./ repmat(sqrt(length(LGN(:, 1))), [1 24]);
    TRN_LO_SEM = std(TRN_LO) ./ repmat(sqrt(length(TRN_LO(:, 1))), [1 24]);
    %if isempty(TRN_HI) == 0; 
    %    TRN_HI_SEM = std(TRN_HI) ./ repmat(sqrt(length(TRN_HI(:, 1))), [1 24]);
    %end;

    subplot(5, 2, iter);

    plot(1:24, LGN_MEAN, 'color', [0.6 0.6 0.6], 'linewidth', 2); hold all;
    plot(1:24, LGN_MEAN+LGN_SEM, 'color', [0.6 0.6 0.6], 'linewidth', 1);
    plot(1:24, LGN_MEAN-LGN_SEM, 'color', [0.6 0.6 0.6], 'linewidth', 1);
    plot(1:24, TRN_LO_MEAN, 'color', [0 0 0], 'linewidth', 2);
    plot(1:24, TRN_LO_MEAN+TRN_LO_SEM, 'color', [0 0 0], 'linewidth', 1);
    plot(1:24, TRN_LO_MEAN-TRN_LO_SEM, 'color', [0 0 0], 'linewidth', 1);

    xlim([1 24]);
    ylim([-1 1]);
    lins = linspace(1,24,11);
    ticks = [];
    for x = 1:length(lins)-1;
        ticks = [ticks ((lins(x)+lins(x+1))/2)];
    end;
    set(gca,'XTick', ticks, 'YTick', [-1:1:1]);
    set(gca, 'xTickLabel', [1, 2.5, 5, 7.5, 10, 12, 15, 20, 30, 60]);
    xlabel('Flicker (Hz)', 'FontSize', 12);
    ylabel('Normalized Response (A.U.)', 'FontSize', 12);

    % subplot(6, 2, iter + 1);

    % plot(1:24, LGN_MEAN, 'color', [0.6 0.6 0.6], 'linewidth', 2); hold all;
    % plot(1:24, LGN_MEAN+LGN_SEM, 'color', [0.6 0.6 0.6], 'linewidth', 1);
    % plot(1:24, LGN_MEAN-LGN_SEM, 'color', [0.6 0.6 0.6], 'linewidth', 1);
    % plot(1:24, TRN_HI_MEAN, 'color', [0 0 0], 'linewidth', 2);
    % plot(1:24, TRN_HI_MEAN+TRN_HI_SEM, 'color', [0 0 0], 'linewidth', 1);
    % plot(1:24, TRN_HI_MEAN-TRN_HI_SEM, 'color', [0 0 0], 'linewidth', 1);

    % xlim([1 24]);
    % ylim([-1.1 1.1]);
    % lins = linspace(1,24,11);
    % ticks = [];
    % for x = 1:length(lins)-1;
    %     ticks = [ticks ((lins(x)+lins(x+1))/2)];
    % end;
    % set(gca,'XTick', ticks, 'YTick', [-1:1:1]);
    % set(gca, 'xTickLabel', [1, 2.5, 5, 7.5, 10, 12, 15, 20, 30, 60]);
    % xlabel('Flicker (Hz)', 'FontSize', 12);
    % ylabel('Normalized Response (A.U.)', 'FontSize', 12);

    disp([':::Subject ' int2str(s) ' complete:::']);

    iter = iter + 1;

end

subplot(5, 2, [7 10])

    LGN_GROUP = LGN_GROUP ./ repmat(max(abs(LGN_GROUP)')', [1 24]);
    TRN_LO_GROUP = TRN_LO_GROUP ./ repmat(max(abs(TRN_LO_GROUP)')', [1 24]);

    %LGN_GROUP = LGN_GROUP ./ repmat(std(LGN_GROUP')', [1 24]);
    %TRN_LO_GROUP = TRN_LO_GROUP ./ repmat(std(TRN_LO_GROUP')', [1 24]);

    LGN_GROUP_MEAN = mean(LGN_GROUP);
    TRN_LO_GROUP_MEAN = mean(TRN_LO_GROUP);

    LGN_GROUP_SEM = std(LGN_GROUP) ./ repmat(sqrt(length(LGN_GROUP(:, 1))), [1 24]);
    TRN_LO_GROUP_SEM = std(TRN_LO_GROUP) ./ repmat(sqrt(length(TRN_LO_GROUP(:, 1))), [1 24]);

    plot(1:24, LGN_GROUP_MEAN, 'color', [0.6 0.6 0.6], 'linewidth', 2); hold all;
    plot(1:24, LGN_GROUP_MEAN+LGN_GROUP_SEM, 'color', [0.6 0.6 0.6], 'linewidth', 1);
    plot(1:24, LGN_GROUP_MEAN-LGN_GROUP_SEM, 'color', [0.6 0.6 0.6], 'linewidth', 1);
    plot(1:24, TRN_LO_GROUP_MEAN, 'color', [0 0 0], 'linewidth', 2);
    plot(1:24, TRN_LO_GROUP_MEAN+TRN_LO_GROUP_SEM, 'color', [0 0 0], 'linewidth', 1);
    plot(1:24, TRN_LO_GROUP_MEAN-TRN_LO_GROUP_SEM, 'color', [0 0 0], 'linewidth', 1);

    xlim([1 24]);
    ylim([-1 1]);
    lins = linspace(1,24,11);
    ticks = [];
    for x = 1:length(lins)-1;
        ticks = [ticks ((lins(x)+lins(x+1))/2)];
    end;
    set(gca,'XTick', ticks, 'YTick', [-1:1:1]);
    set(gca, 'xTickLabel', [1, 2.5, 5, 7.5, 10, 12, 15, 20, 30, 60]);
    xlabel('Flicker (Hz)', 'FontSize', 12);
    ylabel('Normalized Response (A.U.)', 'FontSize', 12);

%% Sepy 20th 2013, JDV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%