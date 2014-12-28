%% descriptive stats on tremotopy

%% create stat mask                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Input statistical map                                                        %
%% %%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
o.dir = '/srv/data/MSC_JNEUROSCI/';             % experiment directory
o.sub = {'01','02','05','07','09','11'};

% per experiment
%o.exp = {'flicker_120_thal'};     % experiment
o.stats = 'stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHead';
o.flicker = [1, 9]; %% phase
o.retino = [11, 19]; %% phase
o.mask = 'mask_LGN_TRN_resamp_split';

o.rotate = 1;
o.ratio = 60; %degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labels = {'R LGN', 'L LGN', 'R vTRN', 'L vTRN', 'R dTRN', 'L dTRN'};
x_labels = {'1', '2.5', '5', '7.5', '10', '12', '15', '20', '30', '60'};

LGN_R_R = [];
LGN_L_R = [];
TRN_R_R = [];
TRN_L_R = [];
PUL_R_R = [];
PUL_L_R = [];

LGN_R_F = [];
LGN_L_F = [];
TRN_R_F = [];
TRN_L_F = [];
PUL_R_F = [];
PUL_L_F = [];

% keeps track of how much each participant contributes to each histogram. 
LGN_COUNTS = [];
vTRN_COUNTS = [];
dTRN_COUNTS = [];

LGN_ROSE_MAX = [];
TRN_ROSE_MAX = [];
PUL_ROSE_MAX = [];
ALL_BINS_MAX = [];

figure;

for s = o.sub;    
	%directory = [o.dir 's' int2str(s) '/' o.exp{1}];
	directory = [o.dir 's' char(s)];
    disp(directory)

    %% Import Mask & Data, Compute Dimensions
    mask = load_nifti([directory '/' o.mask '.nii.gz']);
    stat = load_nifti([directory '/' o.stats '.nii.gz']);
    dims = size(mask.vol);
    nVox = dims(1)*dims(2)*dims(3);

    mask = reshape(mask.vol, nVox, 1);
    retino = stat.vol(:, :, :, o.retino(1));
    r_thresh = stat.vol(:, :, :, o.retino(2));
    retino = reshape(retino, nVox, 1);
    r_thresh = reshape(r_thresh, nVox, 1);

    flicker = stat.vol(:, :, :, o.flicker(1));
    f_thresh = stat.vol(:, :, :, o.flicker(2));
    flicker = reshape(flicker, nVox, 1);
    f_thresh = reshape(f_thresh, nVox, 1);



    count = 1;

    for x = 1:length(unique(mask))-1;
    	THR_R = find(r_thresh == 1);
    	THR_F = find(f_thresh == 1);
    	MSK = find(mask == x);

    	ROI_R = retino(intersect(THR_R, MSK));
    	ROI_F = flicker(intersect(THR_F, MSK));
    	ROI_F = ROI_F / (2*pi);

    	if count == 1;
    		LGN_R_R = [LGN_R_R, ROI_R'];
			LGN_R_F = [LGN_R_F, ROI_F'];
			count_LGN = length(ROI_R);

		elseif count == 2;
    		LGN_L_R = [LGN_L_R, ROI_R'];
			LGN_L_F = [LGN_L_F, ROI_F'];
			count_LGN = count_LGN + length(ROI_R);

		elseif count == 3;
    		TRN_R_R = [TRN_R_R, ROI_R'];
			TRN_R_F = [TRN_R_F, ROI_F'];
			count_vTRN = length(ROI_R);

		elseif count == 4;
    		TRN_L_R = [TRN_L_R, ROI_R'];
			TRN_L_F = [TRN_L_F, ROI_F'];
			count_vTRN = count_vTRN + length(ROI_R);

		elseif count == 7;
    		PUL_R_R = [PUL_R_R, ROI_R'];
			PUL_R_F = [PUL_R_F, ROI_F'];
			count_dTRN = length(ROI_R);

		elseif count == 8;
    		PUL_L_R = [PUL_L_R, ROI_R'];
			PUL_L_F = [PUL_L_F, ROI_F'];
			count_dTRN = count_dTRN + length(ROI_R);
		end

		count = count + 1;
	end

    % add counts to mega thing (for pie charts)
	LGN_COUNTS = [LGN_COUNTS; count_LGN];
    vTRN_COUNTS = [vTRN_COUNTS; count_vTRN];
    dTRN_COUNTS = [dTRN_COUNTS; count_dTRN];

end

count = 1;
for x = 1:length(unique(mask))-1;
	if count == 1;
		[a, b] = rose(LGN_R_R, 36);
		LGN_ROSE_MAX = [LGN_ROSE_MAX, max(b)];
	elseif count == 2;
		[a, b] = rose(LGN_L_R, 36);
		LGN_ROSE_MAX = [LGN_ROSE_MAX, max(b)];
	elseif count == 3;
		[a, b] = rose(TRN_R_R, 36);
		TRN_ROSE_MAX = [TRN_ROSE_MAX, max(b)];
	elseif count == 4;
		[a, b] = rose(TRN_L_R, 36);
		TRN_ROSE_MAX = [TRN_ROSE_MAX, max(b)];
	elseif count == 7;
		[a, b] = rose(PUL_R_R, 36);
		PUL_ROSE_MAX = [PUL_ROSE_MAX, max(b)];
	elseif count == 8;
		[a, b] = rose(PUL_L_R, 36);
		PUL_ROSE_MAX = [PUL_ROSE_MAX, max(b)];
	end
	count = count + 1;
end

figure;
subplot(1,3,1);
pie(LGN_COUNTS, {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'})
colormap('bone')
title('LGN')

subplot(1,3,2);
pie(vTRN_COUNTS, {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'})
colormap('bone')
title('vTRN')

subplot(1,3,3);
pie(dTRN_COUNTS, {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'})
colormap('bone')
title('dTRN')

%% ROSE PLOTS
if o.rotate == 1;
	ratio = o.ratio/360;
end;

% LGN
subplot(2, 3, 1);
polar(max(LGN_ROSE_MAX));
view([-270 90])
hold all;

if o.rotate == 1;
	LGN_R_R = LGN_R_R - (2*pi*ratio);
	LGN_R_R(LGN_R_R < 0) = LGN_R_R(LGN_R_R < 0) + (2*pi);
end
hline = rose(LGN_R_R, 36);
set(hline,'LineWidth', 1, 'Color', [0,0,0]);

if o.rotate == 1;
	LGN_L_R = LGN_L_R - (2*pi*ratio);
	LGN_L_R(LGN_L_R < 0) = LGN_L_R(LGN_L_R < 0) + (2*pi);
end
hline = rose(LGN_L_R, 36);
set(hline,'LineWidth', 1, 'Color', [1,0,0]);
hold off;
title('LGN')

ks_LGN_R = [LGN_R_R, LGN_L_R];
ks_LGN_F = [LGN_R_F, LGN_L_F];

%% TRN ROSE plot
subplot(2, 3, 2);
polar(max(TRN_ROSE_MAX));
view([-270 90])
hold all;

if o.rotate == 1;
	TRN_R_R = TRN_R_R - (2*pi*ratio);
	TRN_R_R(TRN_R_R < 0) = TRN_R_R(TRN_R_R < 0) + (2*pi);
end
hline = rose(TRN_R_R, 36);
set(hline,'LineWidth', 1, 'Color', [0,0,0]);

if o.rotate == 1;
	TRN_L_R = TRN_L_R - (2*pi*ratio);
	TRN_L_R(TRN_L_R < 0) = TRN_L_R(TRN_L_R < 0) + (2*pi);
end
hline = rose(TRN_L_R, 36);
set(hline,'LineWidth', 1, 'Color', [1,0,0]);
hold off;
title('TRN')

ks_TRN_R = [TRN_R_R, TRN_L_R];
ks_TRN_F = [TRN_R_F, TRN_L_F];

%% PUL ROSE plot
subplot(2, 3, 3);
polar(max(PUL_ROSE_MAX));
view([-270 90])
hold all;

if o.rotate == 1;
	PUL_R_R = PUL_R_R - (2*pi*ratio);
	PUL_R_R(PUL_R_R < 0) = PUL_R_R(PUL_R_R < 0) + (2*pi);
end
hline = rose(PUL_R_R, 36);
set(hline,'LineWidth', 1, 'Color', [0,0,0]);

if o.rotate == 1;
	PUL_L_R = PUL_L_R - (2*pi*ratio);
	PUL_L_R(PUL_L_R < 0) = PUL_L_R(PUL_L_R < 0) + (2*pi);
end
hline = rose(PUL_L_R, 36);
set(hline,'LineWidth', 1, 'Color', [1,0,0]);
hold off;
title('PUL')

ks_PUL_R = [PUL_R_R, PUL_L_R];
ks_PUL_F = [PUL_R_F, PUL_L_F];

%% LGN histogram
subplot(2, 3, 4);
bins1 = histc(LGN_R_F, linspace(0, 1, 11));
bins2 = histc(LGN_L_F, linspace(0, 1, 11));
n1 = sum(bins1);
n2 = sum(bins2);
bins1 = bins1 / n1;
bins2 = bins2 / n2;

ALL_BINS_MAX = [ALL_BINS_MAX, max(bins1), max(bins2)];

h = bar(linspace(0, 1, 11) , [bins1', bins2'],'BarWidth', 1);
set(h(1), 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0])
set(h(2), 'FaceColor', [1 0 0], 'EdgeColor', [1 0 0])

title(['n right = ' num2str(n1) ', n left = ' num2str(n2)]);
xlim([-0.1 1])
set(gca, 'xTick', [0:0.1:0.9])
set(gca, 'xTickLabel', x_labels)
xlabel('Flicker (Hz)')
ylabel('Proportion of Total Voxels (%)')

%% TRN histogram
subplot(2, 3, 5);
bins1 = histc(TRN_R_F, linspace(0, 1, 11));
bins2 = histc(TRN_L_F, linspace(0, 1, 11));
n1 = sum(bins1);
n2 = sum(bins2);
bins1 = bins1 / n1;
bins2 = bins2 / n2;

ALL_BINS_MAX = [ALL_BINS_MAX, max(bins1), max(bins2)];

h = bar(linspace(0, 1, 11) , [bins1', bins2'],'BarWidth', 1);
set(h(1), 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0])
set(h(2), 'FaceColor', [1 0 0], 'EdgeColor', [1 0 0])

title(['n right = ' num2str(n1) ', n left = ' num2str(n2)]);
xlim([-0.1 1])
set(gca, 'xTick', [0:0.1:0.9])
set(gca, 'xTickLabel', x_labels)
xlabel('Flicker (Hz)')
ylabel('Proportion of Total Voxels (%)')

%% PUL histogram
subplot(2, 3, 6);
bins1 = histc(PUL_R_F, linspace(0, 1, 11));
bins2 = histc(PUL_L_F, linspace(0, 1, 11));
n1 = sum(bins1);
n2 = sum(bins2);
bins1 = bins1 / n1;
bins2 = bins2 / n2;

ALL_BINS_MAX = [ALL_BINS_MAX, max(bins1), max(bins2)];

h = bar(linspace(0, 1, 11) , [bins1', bins2'],'BarWidth', 1);
set(h(1), 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0])
set(h(2), 'FaceColor', [1 0 0], 'EdgeColor', [1 0 0])

title(['n right = ' num2str(n1) ', n left = ' num2str(n2)]);
xlim([-0.1 1])
set(gca, 'xTick', [0:0.1:0.9])
set(gca, 'xTickLabel', x_labels)
xlabel('Flicker (Hz)')
ylabel('Proportion of Total Voxels (%)')

%% Rescale all histograms
for x = 4:6;
	subplot(2, 3, x);
	ylim([0 max(ALL_BINS_MAX)+0.01])
end

%% KS tests

alpha = 0.05/3;
fid = fopen([o.dir 'KS_tests.txt'], 'w');

[h,p,ks2stat] = kstest2(ks_LGN_R, ks_TRN_R, alpha, 'unequal');
spc = ['retinotopy: LGN vs TRN, p = ' num2str(p) ', stat = ' num2str(ks2stat) '\n'];
fprintf(fid, spc);

[h,p,ks2stat] = kstest2(ks_PUL_R, ks_TRN_R, alpha, 'unequal');
spc = ['retinotopy: TRN vs PUL, p = ' num2str(p) ', stat = ' num2str(ks2stat) '\n'];
fprintf(fid, spc);

[h,p,ks2stat] = kstest2(ks_LGN_R, ks_PUL_R, alpha, 'unequal');
spc = ['retinotopy: LGN vs PUL, p = ' num2str(p) ', stat = ' num2str(ks2stat) '\n'];
fprintf(fid, spc);

[h,p,ks2stat] = kstest2(ks_LGN_F, ks_TRN_F, alpha, 'unequal');
spc = ['flicker: LGN vs TRN, p = ' num2str(p) ', stat = ' num2str(ks2stat) '\n'];
fprintf(fid, spc);

[h,p,ks2stat] = kstest2(ks_PUL_F, ks_TRN_F, alpha, 'unequal');
spc = ['flicker: TRN vs PUL, p = ' num2str(p) ', stat = ' num2str(ks2stat) '\n'];
fprintf(fid, spc);

[h,p,ks2stat] = kstest2(ks_LGN_F, ks_PUL_F, alpha, 'unequal');
spc = ['flicker: LGN vs PUL, p = ' num2str(p) ', stat = ' num2str(ks2stat) '\n'];
fprintf(fid, spc);

fclose(fid);