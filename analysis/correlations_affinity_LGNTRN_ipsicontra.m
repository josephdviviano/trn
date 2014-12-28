% Community Detection of Thalamus%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Options 
% Directory structures: {thalamus, cortex}
opt.dir = '/home/jdv/Dropbox/data/';
opt.sub = [1,2,5,7,9,11];

% file names (no suffix)
opt.datum = 'mean_RUN_bpass_filterHead';
opt.stats = 'stats_FFT_single_multitaper_LGNTRN_mean_RUN_allpass_smoothHeadLGNTRN';
opt.labels = 'mask_LGN_TRN_resamp';

opt.threshvol = 19;
opt.threshval = 0.1;
opt.polarvol = 11;
opt.flickervol = 1;

opt.rthresh = 0

%% Import Data, Stats, & Masks

VEC_CONTRA = [];
VEC_IPSI = [];

for s = opt.sub;

	directory = [opt.dir 's' int2str(s) '/'];

	data = load_nifti([directory '/' opt.datum '.nii.gz']);
	stat = load_nifti([directory '/' opt.stats '.nii.gz']);
	labs = load_nifti([directory '/' opt.labels '.nii.gz']);
	out  = data;
	out.dim(5) = 1;

	dims = size(data.vol);
	numVox = dims(1)*dims(2)*dims(3);
	numTRs = dims(4);

	out.vol = zeros(numVox, 2)';

	% reshape everything
	data = reshape(data.vol, numVox, numTRs)';
	stat = reshape(stat.vol, numVox, 40)';
	labs = reshape(labs.vol, numVox, 1)';

	%polar = stat(opt.polarvol, :);
	flicker = stat(opt.flickervol, :);
	stat = stat(opt.threshvol, :);
	stat(stat <= opt.threshval) = 0;
	stat(stat > opt.threshval) = 1;

	idx_seeds = intersect(find(stat == 1), union(find(labs == 3), find(labs == 4)));
	idx_test  = union(find(labs == 1), find(labs == 2));

	for x = 1:length(idx_seeds);
		r = corr(data(:, idx_seeds(x)), data(:, idx_test));
		if max(r) > opt.rthresh; 
			idx_r = find(r == max(r));
			out.vol(1, idx_seeds(x)) = labs(:, idx_test(idx_r));
			out.vol(2, idx_seeds(x)) = max(r);
			%out.vol(1, idx_test(idx_r)) = labs(:, idx_test(idx_r));
			%out.vol(2, idx_test(idx_r)) = max(r);
		end
	end

	idx_r_r = intersect(find(out.vol(1, :) == 1), find(labs == 3));
	idx_r_l = intersect(find(out.vol(1, :) == 2), find(labs == 3));
	idx_l_r = intersect(find(out.vol(1, :) == 1), find(labs == 4));
	idx_l_l = intersect(find(out.vol(1, :) == 2), find(labs == 4));
	vec_r_r = flicker(:, idx_r_r);
	vec_r_l = flicker(:, idx_r_l);
	vec_l_r = flicker(:, idx_l_r);
	vec_l_l = flicker(:, idx_l_l);

	VEC_IPSI = [VEC_IPSI vec_r_r vec_l_l];
	VEC_CONTRA = [VEC_CONTRA vec_r_l vec_l_r];

	% figure;
	% subplot(2,2,1);
	% hist(vec_r_r, 10);
	% xlim([0 2*pi]);
	% title('flicker r TRN ipsi LGN')

	% subplot(2,2,2);
	% hist(vec_r_l, 10);
	% xlim([0 2*pi]);
	% title('flicker r TRN contra LGN')

	% subplot(2,2,3);
	% hist(vec_l_r, 10);
	% xlim([0 2*pi]);
	% title('flicker l TRN contra LGN')

	% subplot(2,2,4);
	% hist(vec_l_l, 10);
	% xlim([0 2*pi]);
	% title('flicker l TRN ipsi LGN')

	% Write output
	out.vol = reshape(out.vol', dims(1), dims(2), dims(3), 2);
	outStr = sprintf([directory '/stats_correlation_affinity_LGNTRN.nii.gz']);
	save_nifti(out, outStr);
	disp(['s' int2str(s) ' finished'])

end

figure;
subplot(2,1,1);
hist(VEC_IPSI, 10);
xlim([0 2*pi]);
title('flicker TRN corr with ipsi LGN');

subplot(2,1,2);
hist(VEC_CONTRA, 10);
xlim([0 2*pi]);
title('flicker TRN corr with contra LGN');