%% Script to enhance CURE-TSR Dataset

% Clover Thebolt, clover.thebolt@gtri.gatech.edu
% Rohit Chakravarthy, rohit.chakravarthy@gtri.gatech.edu

%% Process CURE-TSR
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

% Load CURE-TSR Dataset
[rawImages, gaussianImages, overexposureImages, blurImages,...
 dirtylensImages, rainImages, shadowImages, snowImages] = loadCURETSR("..\data\CURE-TSR\");

% Number of noise distortions
numDistortions = 7;
distortedImages = [gaussianImages, overexposureImages, blurImages,...
                   dirtylensImages, rainImages, shadowImages, snowImages];              
default_imnoise_var = 0.001*ones(1,numDistortions);               
               
% Initialize metrics object
metrics = DIPMetrics;
metrics.resetMetrics;
               
               
% Create metrics table to populate enhancment results in
% [ Distortion Type X Image X Metric ]
CURETSRMetricsNLM = cell(numDistortions,length(rawImages),7);
CURETSRMetricsBM3D = cell(numDistortions,length(rawImages),7);
CURETSRMetricsBLSGSM = cell(numDistortions,length(rawImages),7);

% Loop through each distortion and enhance image
% Compute metrics for each denoised image
for ii = 1:numDistortions
    
    % Index first Set12 set of distorted images
    currentDistortedImages = distortedImages(:,ii);
    
    % Index noise variance required for BM3D
    noise_var = default_imnoise_var(ii);
    
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
        metrics.calculateAllMetrics(uint8(imgDenoisedNLM),uint8(rawImages{jj,1}));%rawImages{jj,1});
        CURETSRMetricsNLM{ii,jj,1} = metrics.psnr_;
        CURETSRMetricsNLM{ii,jj,2} = metrics.ssim_;
        CURETSRMetricsNLM{ii,jj,3} = metrics.cwssim_;
        CURETSRMetricsNLM{ii,jj,4} = metrics.unique_;
        CURETSRMetricsNLM{ii,jj,5} = metrics.msunique_;
        CURETSRMetricsNLM{ii,jj,6} = metrics.csv_;
        CURETSRMetricsNLM{ii,jj,7} = metrics.summer_;
        
        % BM3D Metrics
        metrics.resetMetrics;
        metrics.calculateAllMetrics(uint8(imgDenoisedBM3D),uint8(rawImages{jj,1}));
        CURETSRMetricsBM3D{ii,jj,1} = metrics.psnr_;
        CURETSRMetricsBM3D{ii,jj,2} = metrics.ssim_;
        CURETSRMetricsBM3D{ii,jj,3} = metrics.cwssim_;
        CURETSRMetricsBM3D{ii,jj,4} = metrics.unique_;
        CURETSRMetricsBM3D{ii,jj,5} = metrics.msunique_;
        CURETSRMetricsBM3D{ii,jj,6} = metrics.csv_;
        CURETSRMetricsBM3D{ii,jj,7} = metrics.summer_;
        
        % BLS-GSM Metrics
        metrics.resetMetrics;
        metrics.calculateAllMetrics(uint8(imgDenoisedBLSGSM),uint8(rawImages{jj,1}));
        CURETSRMetricsBLSGSM{ii,jj,1} = metrics.psnr_;
        CURETSRMetricsBLSGSM{ii,jj,2} = metrics.ssim_;
        CURETSRMetricsBLSGSM{ii,jj,3} = metrics.cwssim_;
        CURETSRMetricsBLSGSM{ii,jj,4} = metrics.unique_;
        CURETSRMetricsBLSGSM{ii,jj,5} = metrics.msunique_;
        CURETSRMetricsBLSGSM{ii,jj,6} = metrics.csv_;
        CURETSRMetricsBLSGSM{ii,jj,7} = metrics.summer_;
        
        
    end
end

% Save metrics
if saveMetrics==1
    save('CURETSRMetricsNLM_var0.001.mat','CURETSRMetricsNLM');
    save('CURETSRMetricsBM3D_var0.001.mat','CURETSRMetricsBM3D');
    save('CURETSRMetricsBLSGSM_var0.001.mat','CURETSRMetricsBLSGSM');
end

    
toc

