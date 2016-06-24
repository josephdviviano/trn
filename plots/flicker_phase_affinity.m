% Community Detection of Thalamus%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Options 
% Directory structures: {thalamus, cortex}

opt.sub = [1, 2, 5, 7, 9, 11];

opt.dir = '/home/jdv/Dropbox/data/msc/';
opt.clusters = 'stats_cluster_TRN'; 
opt.mask = 'mask_LGN_TRN_resamp_split';
opt.affinity = 'stats_correlation_affinity_LGNTRN';
opt.stats = 'stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHeadLGNTRN';

%% Import Data, Stats, & Masks
vector_ipsi = [];
vector_contra = [];
bins_max = [];

vector_ipsi = [];
vector_contra = [];

for s = opt.sub;

	directory = [opt.dir 's' int2str(s) '/' ];

	clusters = load_nifti([directory '/' opt.clusters '.nii.gz']);
	affinity = load_nifti([directory '/' opt.affinity '.nii.gz']);
	mask = load_nifti([directory '/' opt.mask '.nii.gz']);
	stats = load_nifti([directory opt.stats '.nii.gz']);

	stats.vol = stats.vol(:,:,:,1);

	dims = size(mask.vol);
	numVox = dims(1)*dims(2)*dims(3);

	clusters = reshape(clusters.vol, [numVox, 1]);
	affinity = reshape(affinity.vol(:,:,:,1), [numVox, 1]);
	mask = reshape(mask.vol, [numVox, 1]);
	stats = reshape(stats.vol, [numVox, 1])';

	% ipsi
	idx_affinity = find(affinity == 1);
	idx_mask = find(mask == 3);
	idx = intersect(idx_affinity, idx_mask);
	vector_ipsi = [vector_ipsi stats(idx)];

	idx_affinity = find(affinity == 2);
	idx_mask = find(mask == 4);
	idx = intersect(idx_affinity, idx_mask);
	vector_ipsi = [vector_ipsi stats(idx)];

	% contra
	idx_affinity = find(affinity == 1);
	idx_mask = find(mask == 4);
	idx = intersect(idx_affinity, idx_mask);
	vector_contra = [vector_contra stats(idx)];

	idx_affinity = find(affinity == 2);
	idx_mask = find(mask == 3);
	idx = intersect(idx_affinity, idx_mask);
	vector_contra = [vector_contra stats(idx)];
end

bins_contra = histc(vector_contra, linspace(0, 2*pi+(2*pi/10), 11));
bins_ipsi = histc(vector_ipsi, linspace(0, 2*pi+(2*pi/10), 11));

n = sum(bins_contra) + sum(bins_ipsi);

bins_contra = bins_contra / n;
bins_ipsi = bins_ipsi / n;

bar(linspace(0, 2*pi+(2*pi/10), 11) , [bins_contra', bins_ipsi'],'BarWidth', 1, 'FaceColor', [0 0 0], 'EdgeColor', [1 1 1]);
title(['n = ' num2str(n)]);
xlim([-1 2*pi+1])
ylim([0 max([bins_ipsi bins_contra])+0.01])
set(gca, 'xTick', linspace(0, 2*pi+(2*pi/10), 11))