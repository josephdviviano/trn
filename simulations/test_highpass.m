%% Options

opt.input = 'mean_ROTFLICK.nii.gz';
opt.cycle1 = 7;
opt.cycle2 = 10;

mask = load_nifti('mask_LGN_FFT.nii.gz');
data = load_nifti(opt.input);

dims = size(data.vol);
numVox = dims(1)*dims(2)*dims(3);
numTRs = dims(4);

mask = reshape(mask.vol, numVox, 1);
data = reshape(data.vol, numVox, numTRs);

IDX = find(mask > 0);

figure;
hold on;

cmap.back = bone(round(length(IDX)*3/2));
cmap.frnt1 = autumn(round(length(IDX)*2));
cmap.frnt2 = summer(round(length(IDX)*2));
cmap.axis = [0.3, 0.3, 0.3];

%% FFT
dataFFT = fft(data, [], 2);
scaledAmp = abs(dataFFT(IDX,1:numTRs/2));
tsPower = sqrt(sum((scaledAmp(:, 1:numTRs/2).^2), 2));

r1 = scaledAmp(:, opt.cycle1+1) ./tsPower;
r1(isnan(r1)) = 0;

r2 = scaledAmp(:, opt.cycle2+1) ./tsPower;
r2(isnan(r2)) = 0;
    
subplot(2,2,[1 2]);
hold all;
 
for VOX = 1:length(IDX);
    
    stem(1:numTRs/2, scaledAmp(VOX, :), 'LineWidth', 1.5, 'color', cmap.back(VOX, :));
    stem(opt.cycle1+1, scaledAmp(VOX, opt.cycle1+1), 'LineWidth', 2, 'color', cmap.frnt1(VOX, :));
    stem(opt.cycle2+1, scaledAmp(VOX, opt.cycle2+1), 'LineWidth', 2, 'color', cmap.frnt2(VOX, :));
    set(gcf,'color','w');
    set(gca,'FontSize', 12);
    set(gca,'xColor', cmap.axis);
    set(gca,'yColor', cmap.axis);
    xlim([1 numTRs/2]);
    xlabel('Cycles per run', 'Color', cmap.axis);
    ylabel('Power', 'Color', cmap.axis);
    
end

    subplot(2,2,3);
    hold all;
    
    title('r values at 1st cycle frequency ', 'Color', cmap.axis,'FontSize', 12);
    hist(r1, 20);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', 'r', 'EdgeColor', 'w')
    set(gcf,'color','w');
    set(gca,'FontSize', 12);
    set(gca,'xColor', cmap.axis);  
    set(gca,'yColor', cmap.axis);
    set(gca, 'xlim', [0 1]);
    xlabel('correlation (r)', 'Color', cmap.axis);
    ylabel('Count', 'Color', cmap.axis);
    
	subplot(2,2,4);
    hold all;

    title('r values at 2nd cycle frequency ', 'Color', cmap.axis, 'FontSize', 12);
    hist(r2, 20);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', 'g', 'EdgeColor', 'w')
    set(gcf,'color','w');
    set(gca,'FontSize', 12);
    set(gca,'xColor', cmap.axis);
    set(gca,'yColor', cmap.axis);
    set(gca, 'xlim', [0 1]);
    xlabel('correlation (r)', 'Color', cmap.axis);

