%% stimulus plot
opt.freqs = [2, 5, 10, 15, 20, 25, 30, 40, 60, 120];
opt.time = 210;
opt.rotPeriod = 21;
opt.flickPeriod = 30;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Generate stimulus vectors
time = 1:opt.time*1000;
rot = sin(2*pi*time/(opt.rotPeriod*1000));

flickCycles = opt.time/opt.flickPeriod;
flickUnits = opt.flickPeriod / length(opt.freqs);

flickCycle = [];
for freq = opt.freqs;
    
    flickUnit = repmat(freq, [1, flickUnits*1000]);
    flickCycle = horzcat(flickCycle, flickUnit);
    
end
 
flick = repmat(flickCycle, [1 flickCycles]);

%% Continuously Varying Sin Wave
figure;
cmapAxis = [0.3 0.3 0.3];
cmapFrnt = repmat((1:length(flickCycle(1,:))), [1 flickCycles]);
cmapBack = bone(10);

surface([time(1, :); time(1, :)], ...
        [rot(1, :); rot(1, :)], ...
        [flick(1, :); flick(1, :)], ...
        [cmapFrnt; cmapFrnt], ...
        'facecol','no', 'edgecol','interp', 'linew',4);

% global settings
box off
colormap(cmapBack);
axis([1 length(time) -1 1]);
set(gcf,'color','w');
set(gca,'FontSize', 14);
set(gca,'xColor', cmapAxis, 'yColor', cmapAxis, 'zColor', cmapAxis);
title('Dual Frequency Stimulus','Color', cmapAxis, 'FontSize', 14);

% y axis
set(gca,'YTick',[-1 1], ...
        'YTickLabel',{'-pi', '+pi'});
    %, ...
    %    'YGrid', 'on', 'YMinorTick', 'on');
ylabel('Rotation (radians)', 'Color', cmapAxis);

% x axis
set(gca,'XTick', [1, (1:10).*(opt.rotPeriod*1000)], ...
        'XTickLabel',[(0:10) .* opt.rotPeriod]);
    %, ...
    %    'XGrid', 'on', 'XMinorTick', 'on');
xlabel('Time (seconds)', 'Color', cmapAxis);

% colorbar
yTck = (1:length(opt.freqs)) .* length(flickUnit);
cbar = colorbar('YTick', yTck - yTck(1)/2 , ...
                'YTickLabel', opt.freqs, ...
                'yColor', cmapAxis, 'xColor', cmapAxis, ...
                'FontSize', 14);
cbarTitle = ylabel(cbar, 'Flicker Frequency (Hz)');
pos = get(cbarTitle,'position'); 
pos(1,1) = pos(1,1)+2;
set(cbarTitle,'position', pos);
set(cbarTitle,'Rotation',270);
set(gcf, 'renderer', 'zbuffer');