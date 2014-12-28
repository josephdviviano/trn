%% options
opt.time = 240;
opt.cycles = 8:20;
opt.TR = 1;
opt.HRF = 2; % 1 = AFNI single gamma, 2 = SPM double gamma

timeRun = 1:opt.time / opt.TR;

%% exclude bad freqs
IDX = find(rem(opt.time, opt.cycles) == 0);
opt.cycles = opt.cycles(IDX);

%% plot surviving periods
figure;
subplot(2,1,1);
cycVec = [];
hold all;

for cycle = opt.cycles;
    
    period = opt.time / cycle / opt.TR;
   
    if rem(period/2, 1) == 0;
        cycVec = [cycVec, cycle];
    end
end

cmap.frnt = winter(length(cycVec));
cmap.axis = [0.3, 0.3, 0.3];

%% plot all HRFs using a single gamma
HRF = hrf_gamma_hrf(opt.TR);

cmap.counter = 1;

for cycle = cycVec;

    period = opt.time / cycle / opt.TR;
    
    TS = [ones(1, period/2) zeros(1, period/2)];
    TS = repmat(TS, [1, cycle]);
    TS = conv(TS, HRF, 'same');
    
    plot(timeRun(40:80), TS(40:80), 'LineWidth', 2, 'color', cmap.frnt(cmap.counter,:));
    cmap.counter = cmap.counter + 1;
    hold all;
    
end

title(sprintf('run length = %g seconds \n single gamma HRF (AFNI)', length(timeRun)), 'FontSize', 18);
set(gcf,'color','w');
set(gca,'FontSize', 14);
set(gca,'xColor', cmap.axis);  
set(gca,'yColor', cmap.axis);
set(gca, 'xlim', [40 80]); 


%% plot all HRFs using a double gamma
subplot(2,1,2);
[HRF, p] = spm_hrf(opt.TR);
 
cmap.counter = 1;

for cycle = cycVec;

    period = opt.time / cycle / opt.TR;
    
    TS = [ones(1, period/2) zeros(1, period/2)];
    TS = repmat(TS, [1, cycle]);
    TS = conv(TS, HRF, 'same');
    
    plot(timeRun(40:80), TS(40:80), 'LineWidth', 2, 'color', cmap.frnt(cmap.counter,:));
    cmap.counter = cmap.counter + 1;
    hold all;
    
end

title(sprintf('double gamma HRF (SPM)', length(timeRun)), 'FontSize', 18);
set(gcf,'color','w');
set(gca,'FontSize', 14);
set(gca,'xColor', cmap.axis);  
set(gca,'yColor', cmap.axis);
set(gca, 'xlim', [40 80]);
box off

legMat = cell(size(cycVec));
    
counter = 1;
for txt = cycVec;
    legMat{counter} = strcat('cycles=', num2str(txt));
    counter = counter + 1;
end
    
columnlegend(cycVec,legMat);
