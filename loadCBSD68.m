function [cbsd_raw, cbsd_gaussian, cbsd_poisson, cbsd_saltpepper, cbsd_speckle] = loadCBSD68(dataFolder)

% Example usage:
% [cbsd_raw, cbsd_gaussian, cbsd_poisson, cbsd_saltpepper, cbsd_speckle] = loadCBSD68("data\CBSD68-dataset");

imageFolder = "CBSD68\original_png";
imageFiles = dir(fullfile(dataFolder,imageFolder,"*.png"));

cbsd_raw = cell(length(imageFiles),1);
cbsd_gaussian = cell(length(imageFiles),1);
cbsd_poisson = cell(length(imageFiles),1);
cbsd_saltpepper = cell(length(imageFiles),1);
cbsd_speckle = cell(length(imageFiles),1);


for ii = 1:length(imageFiles)
    imagePath = fullfile(dataFolder,imageFolder, imageFiles(ii).name);
    img = imread(imagePath);
    img = rgb2gray(img);
    
    % Add noise to image
    
    % Add zero mean gaussian,var = 0.01 white noise
    img_g = imnoise(img,'gaussian');
    
    % Add poisson noise to data
    img_p = imnoise(img,'poisson');
    
    % Add salt and pepper noise to data
    img_sp = imnoise(img,'salt & pepper');
    
    % Add speckle noise to data
    img_speckle = imnoise(img,'speckle');
    
    % Append images to respective cell
    cbsd_raw{ii,1} = img;
    cbsd_gaussian{ii,1} = img_g;
    cbsd_poisson{ii,1} = img_p;
    cbsd_saltpepper{ii,1} = img_sp;
    cbsd_speckle{ii,1} = img_speckle;

end
end

