function [set12_raw, set12_gaussian, set12_poisson, set12_saltpepper, set12_speckle, set12_blur, set12_highintensity, set12_lowintensity, set12_underexposure, set12_overexposure, set12_contrast] = loadSet12(dataFolder)

% Example Usage: 
%[raw,gauss,poiss,sp,speckle] = loadSet12("data\Set12\Set12");
%

% Noise is added to these images since they are originally noiseless

% Load image files
imageFiles = dir(fullfile(dataFolder, '*.png'));  

set12_raw = cell(length(imageFiles),1);
set12_gaussian = cell(length(imageFiles),1);
set12_poisson = cell(length(imageFiles),1);
set12_saltpepper = cell(length(imageFiles),1);
set12_speckle = cell(length(imageFiles),1);
set12_blur = cell(length(imageFiles),1);
set12_highintensity = cell(length(imageFiles),1);
set12_lowintensity =  cell(length(imageFiles),1);
set12_underexposure = cell(length(imageFiles),1);
set12_overexposure = cell(length(imageFiles),1);
set12_contrast = cell(length(imageFiles),1);

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

    set12_raw{ii,1} = double(img);
    set12_gaussian{ii,1} = double(img_g);
    set12_poisson{ii,1} = double(img_p);
    set12_saltpepper{ii,1} = double(img_sp);
    set12_speckle{ii,1} = double(img_speckle);
    set12_blur{ii,1} = double(img_lightblur);
    set12_highintensity{ii,1} = double(img_highintensity);
    set12_lowintensity{ii,1} = double(img_lowintensity);
    set12_underexposure{ii,1} = double(img_underexposure);
    set12_overexposure{ii,1} = double(img_overexposure);
    set12_contrast{ii,1} = double(img_contrast);
    
end

end