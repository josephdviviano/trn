%% multitaper test

nw     = 1.25;
hr     = 10;
TR     = 1.25;
lp     = 0;
cycles = [7, 10];

ts = dlmread('LGN_mean.1D');
t = 1:length(ts);

pxx_ft    = abs(fft(ts));
pxx_ft    = pxx_ft(2:85);
pxx_mst   = pmtm(ts, nw, length(ts));
pxx_mstHR = pmtm(ts, nw, length(ts)*hr);

pxx_mst   = pxx_mst(2:length(pxx_mst));
pxx_mstHR = pxx_mstHR(2:length(pxx_mstHR));

iter = 1;
bin  = 1;
pxx_avg = [];
while bin <= length(pxx_ft);
    pxx_avg(bin) = mean(pxx_mstHR(iter:iter+hr-1));
    bin = bin+1;
    iter = iter+hr;
end

pxx_ft  = pxx_ft  / max(pxx_ft);
pxx_mst = pxx_mst / max(pxx_mst);
pxx_avg = pxx_avg / max(pxx_avg);

%% Fourier Analysis of Data (Coherence / Voxel with o.cycles Frequencies)
ny = 1/TR/2;
if lp == 0;
    lp = ny;
end
lowRatio = lp / ny;
nBin = ceil(length(ts)/2 * lowRatio);
bins = 1:nBin;

idxS = union(union(cycles-1, cycles), cycles+1);
idxH = union(union(union((cycles*2)-1, (cycles*2)), (cycles*2)+1), ...
             union(union((cycles*3)-1, (cycles*3)), (cycles*3)+1));
idxN = setdiff(bins, union(idxS, idxH));

i = 1;
F=[];
p=[];
for c = cycles;
    F(1, i) = pxx_ft(c).^2  / (sum(pxx_ft(idxN).^2)  / (nBin-length(union(idxS, idxH))));
    F(2, i) = pxx_mst(c).^2 / (sum(pxx_mst(idxN).^2) / (nBin-length(union(idxS, idxH))));
    F(3, i) = pxx_avg(c).^2 / (sum(pxx_avg(idxN).^2) / (nBin-length(union(idxS, idxH))));
    
    p(1, i) = 1-fcdf(F(1, i), 2, 2*(nBin-length(union(idxS, idxH))));
    p(2, i) = 1-fcdf(F(2, i), 2, 2*(nBin-length(union(idxS, idxH))));
    p(3, i) = 1-fcdf(F(3, i), 2, 2*(nBin-length(union(idxS, idxH))));
    
    i = i+1;
end
    
subplot(4,1,1);
plot(t, ts);
title('timeseries');

subplot(4,1,2);
stem(pxx_ft); 
title('fft')

subplot(4,1,3);
stem(pxx_mst);
title('mpsde');

subplot(4,1,4);
stem(pxx_avg);
title('mpsde --> hr --> avg');