%% Script to enhance SIDD Dataset

% Clover Thebolt, clover.thebolt@gtri.gatech.edu
% Rohit Chakravarthy, rohit.chakravarthy@gtri.gatech.edu

%% Process SIDD
tic

% Adding metrics code path for evaluation
addpath("metricsCode\matlabPyrTools\");
addpath("metricsCode\InputWeights\");
addpath("metricsCode\csvFuncs\");
addpath("metricsCode\csvFuncs\FastEMD\");
addpath("enhancement_methodsCode\bm3d_matlab_package\bm3d\");
addpath("enhancement_methodsCode\BLS-GSM-WaveletDenoising_toolboxes\");

% NOTE: there is a "psnr.m" function in the
% BLS-GSM-WaveletDenoising_toolboxes\toolbox_signal folder. I renamed it to
% "psnr__.m" and also renamed the snr.m to "snr__.m" so that we use the
% MATLAB built in functions rather than the ones in that folder


% Save mat files 
saveMetrics = 1;

% Load SIDD Dataset
[rawImagesGaussian, rawImagesUnderexposure, gaussianImages, underexposureImages] = loadSIDD("..\data\SIDD\SIDD_Small_Raw_Only\");
disp('All SIDD data successfully loaded');
% Number of noise distortions
numDistortions = 2;
distortedImages = [underexposureImages, gaussianImages];              
default_imnoise_var = 0.001*ones(1,numDistortions);               
               
% Initialize metrics object
metrics = DIPMetrics;
metrics.resetMetrics;

% Create metrics table to populate enhancment results in
% [ Distortion Type X Image X Metric ]
SIDDMetricsNLM = cell(numDistortions,length(rawImagesGaussian),7);
SIDDMetricsBM3D = cell(numDistortions,length(rawImagesGaussian),7);
SIDDMetricsBLSGSM = cell(numDistortions,length(rawImagesGaussian),7);

% Loop through each distortion and enhance image
% Compute metrics for each denoised image
for ii = 1:numDistortions
    
    % Index first Set12 set of distorted images
    currentDistortedImages = distortedImages(:,ii);
    
    % Index noise variance required for BM3D
    noise_var = default_imnoise_var(ii);
    
    % Set rawImages variable based on which noise type is being processed first
    % This is to make sure the metrics are calculated with the correct
    % original image
    if ii == 1
        rawImages = rawImagesUnderexposure;
    elseif ii == 2
        rawImages = rawImagesGaussian;
    end
    
    for jj = 1:length(currentDistortedImages)
        
        % Apply each Enhancement method on each distorted image in Set12
        imgDistorted = currentDistortedImages{jj,1};
        
        % Non-Local-Means Enhancement
        dims = size(imgDistorted);
        %roi = [210,24,52,41];
        %roi = [128, 56, 80,41];
        roi = [floor(dims(1)/2),floor(dims(2)/2),80,41];
        patch = imcrop(imgDistorted,roi);  % Extract a homogenous LAB patch to compute the noise std dev
        patchSq = patch.^2;
        edist = sqrt(sum(patchSq,3));  % Compute the Euclidean distance from the origin
        patchSigma = sqrt(var(edist(:)));  % Calculate the standard deviation of edist to estimate the noise
        DoS = 1.5*patchSigma;  % Set the 'DegreeOfSmoothing' value to be higher than the standard deviation of the patch.
        imgDenoisedNLM = imnlmfilt(imgDistorted,'DegreeOfSmoothing',DoS);  % Filter the noisy L*a*b* image using non-local means filtering.
        
        % BM3D Enhancement
        imgDenoisedBM3D = BM3D(imgDistorted, sqrt(noise_var));
        
        % BLS-GSM Enhancement
        imgDenoisedBLSGSM = perform_blsgsm_denoising(imgDistorted);
        
        % Compute metrics for each enhancement method
     
        % NLM Metrics
        metrics.resetMetrics;
        metrics.calculateAllMetrics(single(imgDenoisedNLM),single(rawImages{jj,1}));%rawImages{jj,1});
        SIDDMetricsNLM{ii,jj,1} = metrics.psnr_;
        SIDDMetricsNLM{ii,jj,2} = metrics.ssim_;
        SIDDMetricsNLM{ii,jj,3} = metrics.cwssim_;
        SIDDMetricsNLM{ii,jj,4} = metrics.unique_;
        SIDDMetricsNLM{ii,jj,5} = metrics.msunique_;
        SIDDMetricsNLM{ii,jj,6} = metrics.csv_;
        SIDDMetricsNLM{ii,jj,7} = metrics.summer_;
        
        % BM3D Metrics
        metrics.resetMetrics;
        metrics.calculateAllMetrics(single(imgDenoisedBM3D),single(rawImages{jj,1}));
        SIDDMetricsBM3D{ii,jj,1} = metrics.psnr_;
        SIDDMetricsBM3D{ii,jj,2} = metrics.ssim_;
        SIDDMetricsBM3D{ii,jj,3} = metrics.cwssim_;
        SIDDMetricsBM3D{ii,jj,4} = metrics.unique_;
        SIDDMetricsBM3D{ii,jj,5} = metrics.msunique_;
        SIDDMetricsBM3D{ii,jj,6} = metrics.csv_;
        SIDDMetricsBM3D{ii,jj,7} = metrics.summer_;
        
        % BLS-GSM Metrics
        metrics.resetMetrics;
        metrics.calculateAllMetrics(single(imgDenoisedBLSGSM),single(rawImages{jj,1}));
        SIDDMetricsBLSGSM{ii,jj,1} = metrics.psnr_;
        SIDDMetricsBLSGSM{ii,jj,2} = metrics.ssim_;
        SIDDMetricsBLSGSM{ii,jj,3} = metrics.cwssim_;
        SIDDMetricsBLSGSM{ii,jj,4} = metrics.unique_;
        SIDDMetricsBLSGSM{ii,jj,5} = metrics.msunique_;
        SIDDMetricsBLSGSM{ii,jj,6} = metrics.csv_;
        SIDDMetricsBLSGSM{ii,jj,7} = metrics.summer_;
        
        
    end
end

% Save metrics
if saveMetrics==1
    save('SIDDMetricsNLM_var0.001_single.mat','SIDDMetricsNLM');
    save('SIDDMetricsBM3D_var0.001_single.mat','SIDDMetricsBM3D');
    save('SIDDMetricsBLSGSM_var0.001_single.mat','SIDDMetricsBLSGSM');
end

    
toc

