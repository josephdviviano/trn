%% Multi-Frequency Fourier Analysis w/ FDR thresholding %%%%%%%%%%%%%%%%%%%%%
% Input co-registered & detrended phase-encoded images w/ blank TRs deleted %
%% %%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
o.dir = '/Users/jdv/data/trn/';             % experiment directory
o.sub = [1, 2, 5, 7, 9, 11];       

% per experiment
o.exp = {'flicker_120_thal'};     % experiment
o.masks = {'mask_LGN_TRN_nopul'};     % input mask

% data within each experiment
o.datum = {'mean_RUN_allpass_smoothHeadLGNTRN'};

o.stim = [7, 10];           % frequency bins
o.harm = [14, 20];

o.q = [0.25, 0.2, 0.1, 0.05, 0.01];        %  q[FDR]

% time series properties
o.TR = 1.25;                  % TR (seconds)
o.lp = 0;                  % Lowpass frequency in Hz (0 = off)
o.nw = 1.25

%o.sideband = 1;       % 1 = Include sideband analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = o.sub;
    
directory = [o.dir 's' int2str(s) '/' o.exp{1}];

for d = 1:length(o.datum);

    %% Import Mask & Data, Compute Dimensions
    dat = load_nifti([directory '/' o.datum{d} '.nii.gz']);
    msk = load_nifti([directory '/' o.masks{1} '.nii.gz']);

    dim1 = size(dat.vol);
    nVox1 = dim1(1)*dim1(2)*dim1(3);
    nTRs = dim1(4);
    nCyc = length(o.stim) + length(o.harm);

    % Generate output file
    OUT = dat;
    OUT.vol = zeros(dim1(1), dim1(2), dim1(3), 6+length(o.q) * nCyc);
    OUT.dim(5) = 6+length(o.q) * nCyc;

    % Create concatenated analysis file & masks
    dat.vol = reshape(dat.vol, nVox1, nTRs);
    msk.vol = reshape(msk.vol, nVox1, 1);
    
    msk = find(msk.vol > 0);
    dat = dat.vol;
    dat(isnan(dat)) = 0;

    % low-pass the data, and find frequencies of interest
    ny = 1/o.TR/2;
    if o.lp == 0;
        o.lp = ny;
    end
    lowRatio = o.lp / ny;
    nBin = ceil(nTRs/2*lowRatio);
    bins = 2:nBin;

    %% Fourier Analysis of Data (Coherence / Voxel with o.cycles Frequencies)    
    ft  = fft(dat(msk, :), [], 2);
    pxx = [];
    for x = 1:length(msk);   
        pxx(x, :) = pmtm(dat(msk(x), :)', o.nw, nTRs);
    end
    pxx = pxx(:, 2:length(pxx(1, :)));

    % NB: this protect adjacent bins to prevent spectral leakage
    IDXs = union(union(o.stim-1, o.stim), o.stim+1);
    IDXh = union(union(union((o.stim*2)-1, (o.stim*2)), (o.stim*2)+1), ... 
                 union(union((o.stim*3)-1, (o.stim*3)), (o.stim*3)+1));
    IDXn = setdiff(bins-1, union(IDXs, IDXh));

    %% Write FFT stats for all cycles
    iter = 1;
    cycles = union(o.stim, o.harm);
    for c = cycles;
        
        countIter = 0; 
            
        %% phase
        phase = mod(-angle(ft(:, c+1)), 2 * pi);

        %% amplitude
        amp = 2 * (pxx(:, c)) / nTRs;
        amp(isnan(amp)) = 0;
        %amp = sum(amp, 2);
        
        %% r (engle, 1997), F & p Wei & Craigmile 2010: 
        r = pxx(:, c) ./ sqrt(sum(pxx.^2, 2)); 
        F = (nBin-length(union(IDXs, IDXh))) * (pxx(:, c).^2 ./ sum(pxx(:, IDXn).^2, 2));
        p = 1-fcdf(F, 2, 2*(nBin-length(union(IDXs, IDXh))));
    
        %% FDR thresholds
        FDR = zeros(length(p), length(o.q));
        OUTstats = [];
        loop = 1;
        for q = o.q;
            
            [FDRval, threshP] = fdr_bh(p, q);
            FDR(:, loop) = FDRval;
            loop = loop+1;

            disp(['::: q[' num2str(q, 5) '] p < ' num2str(threshP, 5)]);
            
            tmpStats = [q, threshP];
            OUTstats = vertcat(OUTstats, tmpStats);

        end
        
        OUTstr = [directory '/FDR_single_multitaper_LGNTRN_' o.datum{d} '_' int2str(c) '.txt'];
        dlmwrite(OUTstr, OUTstats, 'delimiter', '\t');    
        
        %% write outputs
        OUTphase = zeros(nVox1, 1);
        OUTphase(msk) = phase;
        OUTphase = reshape(OUTphase, dim1(1), dim1(2), dim1(3), 1);
        OUT.vol(:, :, :, iter) = OUTphase;

        OUTamp = zeros(nVox1, 1);
        OUTamp(msk) = amp;
        OUTamp = reshape(OUTamp, dim1(1), dim1(2), dim1(3), 1);
        OUT.vol(:, :, :, iter + 1) = OUTamp;

        OUTr = zeros(nVox1, 1);
        OUTr(msk) = r;
        OUTr = reshape(OUTr, dim1(1), dim1(2), dim1(3), 1);
        OUT.vol(:, :, :, iter + 2) = OUTr;

        OUTF = zeros(nVox1, 1);
        OUTF(msk) = F;
        OUTF = reshape(OUTF, dim1(1), dim1(2), dim1(3), 1);
        OUT.vol(:, :, :, iter + 3) = OUTF;

        OUTp = zeros(nVox1, 1);
        OUTp(msk) = p;
        OUTp = reshape(OUTp, dim1(1), dim1(2), dim1(3), 1);
        OUT.vol(:, :, :, iter + 4) = OUTp;
        
        OUTFDR = zeros(nVox1, length(o.q));
        OUTFDR(msk, :) = FDR;
        OUTFDR = reshape(OUTFDR, dim1(1), dim1(2), dim1(3), length(o.q));
        OUT.vol(:, :, :, iter + 5:iter + 4 + length(o.q)) = OUTFDR;

        iter = iter + (5+length(o.q));
            
    end

    %% save outputs
    OUTstr = [directory '/stats_FFT_single_multitaper_LGNTRN_' o.datum{d} '.nii.gz'];
    save_nifti(OUT, OUTstr);

    disp([':::Datum ' o.datum{d} ' complete:::'])

end

disp([':::Subject ' int2str(s) ' complete:::']);

end
%% Jul 22nd 2013, JDV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%