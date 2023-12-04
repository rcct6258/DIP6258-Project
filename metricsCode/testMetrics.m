% Script to test metrics


% Initialize object
d = DIPMetrics;

% Reset method to clear any values in the DIPMetrics properties
d.resetMetrics;

% Some fake data
% p = abs(ceil(50 + 5.*randn(128,128,3)));
% q = abs(ceil(50 + 5.*randn(128,128,3)));
p = imread('peppers.png');
q = imnoise(p,'gaussian');
% Test each metric individually
psnrMetric = d.PSNR(q,p);
ssimMetric = d.SSIM(q,p);
cwssimMetric = d.CWSSIM(reshape(p,[],size(p,2)),reshape(q,[],size(q,2)),1,16,0,0); % for some reason this one doesn't work with three channels. Has something to do with the third argument which I'm still investigating
uniqueMetric = d.UNIQUE(q,p);
msuniqueMetric = d.MSUNIQUE(q,p);
csvMetric = d.CSV(q,p);
summerMetric = d.SUMMER(q,p);

% Calculates all metrics and populates the DIPMetrics object properties
d.calculateAllMetrics(q,p);




