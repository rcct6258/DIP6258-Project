%function [cureor_raw, cureor_noisy] = loadCUREOR(dataFolder,nImages,challenge,level,noiseType,cameraType)
%function [cureor_raw, cureor_gaussian, cureor_poisson,cureor_saltpepper, cureor_speckle] = loadCUREOR(dataFolder,nImages,challenge,level,noiseType,cameraType)
function [cureor_raw, cureor_underexposure, cureor_overexposure, cureor_blur, cureor_contrast, cureor_dirtylens, cureor_saltpepper] = loadCUREOR(dataFolder)
% Example Usage:
% curor = loadCUREOR("data\CURE-OR");

% No additional noise added to these images since noise is already present

% Note the rawFolder data is always loaded since it is the noiseless case
% and will be needed for comparison 

% Hand-picked images
imageNames = ["1_006","1_007","1_009","1_010","1_046","1_067","1_071","1_072","1_075","1_082"];

% Folder Names
rawFolder = "10_grayscale_no_challenge\white\iPhone\";
underexposureFolder = "12_grayscale_underexposure\Level_1\white\iPhone\";
overexposureFolder = "13_grayscale_overexposure\Level_1\white\iPhone\";
blurFolder = "14_grayscale_blur\Level_1\white\iPhone\";
contrastFolder = "15_grayscale_contrast\Level_1\white\iPhone\";
dirtylensFolder = "16_grayscale_dirtylens1\Level_1\white\iPhone\";
saltpepperFolder = "18_grayscale_saltpepper\Level_1\white\iPhone\";

% Cells to store loaded images
numImages = 10;
cureor_raw = cell(numImages,1);
cureor_underexposure = cell(numImages,1);
cureor_overexposure = cell(numImages,1);
cureor_blur = cell(numImages,1);
cureor_contrast = cell(numImages,1);
cureor_dirtylens = cell(numImages,1);
cureor_saltpepper = cell(numImages,1);


% Load raw data paths
rawImagePath = fullfile(dataFolder,rawFolder);
rawImageFiles = dir(fullfile(rawImagePath,"*.jpg"));
for ii = 1:length(imageNames)
    % Find and load image
    imgFileMatch = string(rawImageFiles(contains({rawImageFiles.name},imageNames(ii))).name);
    imgFile = fullfile(rawImagePath,imgFileMatch);
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    cureor_raw{ii,1} = double(img);
end

% Load underexposure data paths
underexposureImagePath = fullfile(dataFolder,underexposureFolder);
underexposureImageFiles = dir(fullfile(underexposureImagePath,"*.jpg"));
for ii = 1:length(imageNames)
    % Find and load image
    imgFileMatch = string(underexposureImageFiles(contains({underexposureImageFiles.name},imageNames(ii))).name);
    imgFile = fullfile(underexposureImagePath,imgFileMatch);
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    cureor_underexposure{ii,1} = double(img);
end

% Load overexposure data paths
overexposureImagePath = fullfile(dataFolder,overexposureFolder);
overexposureImageFiles = dir(fullfile(overexposureImagePath,"*.jpg"));
for ii = 1:length(imageNames)
    % Find and load image
    imgFileMatch = string(overexposureImageFiles(contains({overexposureImageFiles.name},imageNames(ii))).name);
    imgFile = fullfile(overexposureImagePath,imgFileMatch);
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    cureor_overexposure{ii,1} = double(img);
end

% Load blur data paths
blurImagePath = fullfile(dataFolder,blurFolder);
blurImageFiles = dir(fullfile(blurImagePath,"*.jpg"));
for ii = 1:length(imageNames)
    % Find and load image
    imgFileMatch = string(blurImageFiles(contains({blurImageFiles.name},imageNames(ii))).name);
    imgFile = fullfile(blurImagePath,imgFileMatch);
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    cureor_blur{ii,1} = double(img);
end

% Load contrast data paths
contrastImagePath = fullfile(dataFolder,contrastFolder);
contrastImageFiles = dir(fullfile(contrastImagePath,"*.jpg"));
for ii = 1:length(imageNames)
    % Find and load image
    imgFileMatch = string(contrastImageFiles(contains({contrastImageFiles.name},imageNames(ii))).name);
    imgFile = fullfile(contrastImagePath,imgFileMatch);
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    cureor_contrast{ii,1} = double(img);
end

% Load dirty lens 1 data paths
dirtylensImagePath = fullfile(dataFolder,dirtylensFolder);
dirtylensImageFiles = dir(fullfile(dirtylensImagePath,"*.jpg"));
for ii = 1:length(imageNames)
    % Find and load image
    imgFileMatch = string(dirtylensImageFiles(contains({dirtylensImageFiles.name},imageNames(ii))).name);
    imgFile = fullfile(dirtylensImagePath,imgFileMatch);
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    cureor_dirtylens{ii,1} = double(img);
end

% Load salt & pepper data paths
saltpepperImagePath = fullfile(dataFolder,saltpepperFolder);
saltpepperImageFiles = dir(fullfile(saltpepperImagePath,"*.jpg"));
for ii = 1:length(imageNames)
    % Find and load image
    imgFileMatch = string(saltpepperImageFiles(contains({saltpepperImageFiles.name},imageNames(ii))).name);
    imgFile = fullfile(saltpepperImagePath,imgFileMatch);
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    cureor_saltpepper{ii,1} = double(img);
end






