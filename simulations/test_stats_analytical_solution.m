clear all
N = 168;    % number of TRs in time series
freq = 10;  % frequency of interest
nVox = 10000;
nBin = 1:N/2;

for i = 1:nVox
    ts = normrnd(zeros(1,N), ones(1,N));
    %ts = unifrnd(-1,1, 1, N);
    ps = abs(fft(ts));
    r(i) = ps(freq+1)/sqrt(sum(ps(setdiff(nBin, freq+1)).^2));
end

figure;
subplot(2,2,1);
title('Engle Method')
hist(r,100)
xlabel('r')
rs = 0:.01:.5;

for i=1:length(rs)
    pd(i) = sum(r>=rs(i))/nVox;              % p-value based on noise simulation
    pa(i) = 1-fcdf((N/2-1)*rs(i)^2, 2, N-2); % analytical p-value
end

subplot(2,2,3);
l(1) = semilogy(rs,pd); hold on;
l(2) = semilogy(rs,pa,'r');
xlabel('r')
ylabel('p')
legend(l,{'Numerical histogram method','Analytical F-test method'})

%% Wei Method (protect 6 freqs + adjacent bins)
N = 168;    % number of TRs in time series
freq = [7, 14, 21, 10, 20, 30];  % frequency of interest
nVox = 10000;
nBin = 1:N/2;
IDXprotect = union(union(freq, freq+1), freq+2);
IDXnoise = setdiff(nBin, IDXprotect);

for i = 1:nVox
    ts = normrnd(zeros(1,N), ones(1,N));
    %ts = unifrnd(-1,1, 1, N);
    ps = abs(fft(ts));
    F(i) = ps(:, freq(1)+1).^2 ./ mean(ps(:, IDXnoise).^2, 2);
end

subplot(2,2,2);
title('Wei Method')
hist(F,100)
xlabel('F')
Fs = 0:.1:5;

for i=1:length(rs)
    pd(i) = sum(F>=Fs(i))/nVox;              % p-value based on noise simulation
    pa(i) = 1-fcdf((N/2-1)*Fs(i)^2, 2, N-2); % analytical p-value
end

subplot(2,2,4);
l(1) = semilogy(Fs,pd); hold on;
l(2) = semilogy(Fs,pa,'r');
xlabel('F')
ylabel('p')
legend(l,{'Numerical histogram method','Analytical F-test method'})