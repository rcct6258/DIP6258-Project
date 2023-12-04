function [sample_raw, sample_gaussian, sample_poisson, sample_saltpepper, sample_speckle, sample_blur, sample_highintensity, sample_lowintensity, sample_underexposure, sample_overexposure, sample_contrast] = loadSAMPLE(dataFolder)

% Example Usage: 
%[raw,gauss,poiss,sp,speckle] = loadSet12("data\Set12\Set12");
%

% Noise is added to these images since they are originally noiseless

% Load image files
imageFiles = dir(fullfile(dataFolder, '*.png'));  

sample_raw = cell(length(imageFiles),1);
sample_gaussian = cell(length(imageFiles),1);
sample_poisson = cell(length(imageFiles),1);
sample_saltpepper = cell(length(imageFiles),1);
sample_speckle = cell(length(imageFiles),1);
sample_blur = cell(length(imageFiles),1);
sample_highintensity = cell(length(imageFiles),1);
sample_lowintensity =  cell(length(imageFiles),1);
sample_underexposure = cell(length(imageFiles),1);
sample_overexposure = cell(length(imageFiles),1);
sample_contrast = cell(length(imageFiles),1);

for ii = 1:length(imageFiles)
    imagePath = fullfile(dataFolder, imageFiles(ii).name);
    img = imread(imagePath);
    
    % Compute size of image to check for RGB
    s = size(img);
    
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    
    % Add noise to image
    
    % Add zero mean gaussian,var = 0.01 white noise
    img_g = imnoise(img,'gaussian',0,0.001);
    
    % Add poisson noise to data
    img_p = imnoise(img,'poisson');
    
    % Add salt and pepper noise to data
    img_sp = imnoise(img,'salt & pepper');
    
    % Add speckle noise to data
    img_speckle = imnoise(img,'speckle');
    
    % Add blurring to data
    filt_size = 7;
    light_blur = fspecial('average',filt_size);  % Light blurring effect
    img_lightblur = imfilter(img,light_blur,"replicate");
    
    % Add high intensity noise to data
    img_highintensity = imadjust(img,[0; 0.5]);
    
    % Add low intensity noise to data
    img_lowintensity = imadjust(img,[0.5; 1.0]);
    
    % Add underexposure noise to data
    %img_underexposure = img - 2*std(double(img),0,'all');
    img_underexposure = max(max(img,0)-50,0);
    
    % Add overexposure noise to data
    %img_overexposure = img + 2*std(double(img),0,'all');
    img_overexposure = min(max(img,0)+50,255);
    
    % Add contrast noise to data
    img_contrast = imadjust(img);

    
    % Append images to respective cell

    sample_raw{ii,1} = double(img);
    sample_gaussian{ii,1} = double(img_g);
    sample_poisson{ii,1} = double(img_p);
    sample_saltpepper{ii,1} = double(img_sp);
    sample_speckle{ii,1} = double(img_speckle);
    sample_blur{ii,1} = double(img_lightblur);
    sample_highintensity{ii,1} = double(img_highintensity);
    sample_lowintensity{ii,1} = double(img_lowintensity);
    sample_underexposure{ii,1} = double(img_underexposure);
    sample_overexposure{ii,1} = double(img_overexposure);
    sample_contrast{ii,1} = double(img_contrast);
    
end

end