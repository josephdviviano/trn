opt.TR = 1.25;
opt.loPassHz = 0.25;
opt.filterOrder = 2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_data = dlmread('params_mot1.1D');
cmap.axis = [0.3 0.3 0.3];
cmap.frnt = lines(5);

%% Normalize Input Data
input_dims = size(input_data);
 input_std = std(input_data);
input_mean = mean(input_data);
 input_std = repmat(input_std, [length(input_data), 1]);
input_mean = repmat(input_mean, [length(input_data), 1]);
input_data = (input_data - input_mean) ./ input_std;

%% Butterworth
nyquistHz = 1/opt.TR/2;
[butterB, butterA] = butter(opt.filterOrder, [opt.loPassHz/nyquistHz],'low');
output_data_butter = filtfilt(butterB, butterA, input_data);

% plot the response function
figure('Position', [100 900 1200, 900]);
subplot(3,2,1);
[h,w] = freqz(butterB, butterA, 128);
p = abs(h);
f = w/(pi*2);
plot(f, p.^2, 'LineWidth', 2, 'Color', 'black');
set(gcf,'color','w');
set(gca,'FontSize', 12);
set(gca,'xColor', cmap.axis);
set(gca,'yColor', cmap.axis);
box off
xlabel('Frequency (Hz)', 'FontSize', 12, 'Color', cmap.axis);
xlim([0 nyquistHz]);
set(gca,'XTick', [0:0.1:nyquistHz]);
ylabel('Power', 'FontSize', 12, 'Color', cmap.axis);
ylim([0 1.2]);
set(gca,'YTick', [0:0.1:1.2]);
title('Frequency Response: Butterworth', 'FontSize', 14, 'Color', cmap.axis);

% plot the impulse response function
subplot(3,2,3);
impz(butterB, butterA, 25);
set(gcf,'color','w');
set(gca,'FontSize', 12);
set(gca,'xColor', cmap.axis);
set(gca,'yColor', cmap.axis);
box off
%xlabel('Frequency (Hz)', 'FontSize', 12, 'Color', cmap.axis);
%xlim([0 nyquistHz]);
%set(gca,'XTick', [0:0.1:nyquistHz]);
%ylabel('Power', 'FontSize', 12, 'Color', cmap.axis);
%ylim([0 1.2]);
%set(gca,'YTick', [0:0.1:1.2]);
title('Impulse Response: Butterworth', 'FontSize', 14, 'Color', cmap.axis);

%% define Chebyshev Type II filter, with bi-directional filtering
nyquistHz = 1/opt.TR/2;
[cheb2B, cheb2A] = cheby2(opt.filterOrder, 20, [opt.loPassHz/nyquistHz],'low');
output_data_chebyshev = filtfilt(cheb2B, cheb2A, input_data);

% plot the response function
subplot(3,2,2);
[h,w] = freqz(cheb2B, cheb2A, 128);
p = abs(h);
f = w/(pi*2);
plot(f, p.^2, 'LineWidth', 2, 'Color', 'black');
set(gcf,'color','w');
set(gca,'FontSize', 12);
set(gca,'xColor', cmap.axis);
set(gca,'yColor', cmap.axis);
box off
xlabel('Frequency (Hz)', 'FontSize', 12, 'Color', cmap.axis);
xlim([0 nyquistHz]);
set(gca,'XTick', [0:0.1:nyquistHz]);
ylabel('Power', 'FontSize', 12, 'Color', cmap.axis);
ylim([0 1.2]);
set(gca,'YTick', [0:0.1:1.2]);
title('Frequency Response: Chebyshev II', 'FontSize', 14, 'Color', cmap.axis);
    
% plot the impulse response function
subplot(3,2,4);
impz(cheb2B, cheb2A, 25);
set(gcf,'color','w');
set(gca,'FontSize', 12);
set(gca,'xColor', cmap.axis);
set(gca,'yColor', cmap.axis);
box off
%xlabel('Frequency (Hz)', 'FontSize', 12, 'Color', cmap.axis);
%xlim([0 nyquistHz]);
%set(gca,'XTick', [0:0.1:nyquistHz]);
%ylabel('Power', 'FontSize', 12, 'Color', cmap.axis);
%ylim([0 1.2]);
%set(gca,'YTick', [0:0.1:1.2]);
title('Impulse Response: Chebyshev II', 'FontSize', 14, 'Color', cmap.axis);

% % apply Savitzky-Golay filter
% output_data_SG = repmat(input_data, [3 1]);
% output_data_SG = sgolayfilt(output_data_SG, opt.filterOrder, 5);
% output_data_SG = output_data_SG(1+input_dims(1):2*input_dims(1), 1);
% 
% % apply moving average filter (span = 3)
% output_data_moving_3 = repmat(input_data, [3 1]);
% output_data_moving_3 = smooth(output_data_moving_3, 3);
% output_data_moving_3 = output_data_moving_3(1+input_dims(1):2*input_dims(1), 1);
% 
% % apply moving average filter (span = 5)
% output_data_moving_5 = repmat(input_data, [3 1]);
% output_data_moving_5 = smooth(output_data_moving_5, 5);
% output_data_moving_5 = output_data_moving_5(1+input_dims(1):2*input_dims(1), 1);

subplot(3,2,[5 6]);
plot(1:input_dims(1), input_data, 'LineWidth', 2, 'Color', 'black'); hold all;
plot(1:input_dims(1), output_data_butter, 'LineWidth', 1, 'Color', 'Red');
plot(1:input_dims(1), output_data_chebyshev, 'LineWidth', 1, 'Color', 'Green');
legend('raw', 'color', 'black', 'butterworth', 'color', 'red', 'chebyshev', 'color', 'green');

%plot(1:input_dims(1), output_data_SG, 'LineWidth', 2, 'Color', cmap.frnt(3,:));
%plot(1:input_dims(1), output_data_moving_3, 'LineWidth', 2, 'Color', cmap.frnt(4,:));
%plot(1:input_dims(1), output_data_moving_5, 'LineWidth', 2, 'Color', cmap.frnt(5,:));

title([sprintf('Lowpass at %0.5g Hz', opt.loPassHz) sprintf(', Filter Order %d', opt.filterOrder)], 'FontSize', 12, 'Color', cmap.axis);
set(gcf,'color','w');
set(gca,'FontSize', 12);
set(gca,'xColor', cmap.axis);
set(gca,'yColor', cmap.axis);
xlabel('TR', 'FontSize', 12, 'Color', cmap.axis);
xlim([1 input_dims(1)]);
xTicks =  [0:6:input_dims(1)];
xTickLabels = xTicks;
xTicks(1) = 1;
set(gca,'XTick', xTicks);
set(gca, 'XTickLabel', xTickLabels);
box off
        
