%function [curetsr_raw, curetsr_noisy] = loadCURETSR(dataFolder,nImages,imageCategory,noiseType)
function [curetsr_raw, curetsr_gaussian, curetsr_overexposure, curetsr_blur, curetsr_dirtylens, curetsr_rain, curetsr_shadow, curetsr_snow] = loadCURETSR(dataFolder)




% Example Usage:
% [raw,noisy] = loadCURETSR("data\CURE-TSR",100,"Real_Test","DirtyLens-1");


% Hand-picked images

rawFileNames = ["01_01_00_00_0147","01_01_00_00_0538","01_05_00_00_0050",...
                 "01_05_00_00_0226","01_05_00_00_0360","01_06_00_00_0041",...
                 "01_07_00_00_0049","01_08_00_00_0091","01_10_00_00_0143",...
                 "01_14_00_00_0372"];
             
gaussianFileNames = ["01_01_08_01_0147","01_01_08_01_0538","01_05_08_01_0050",...
                      "01_05_08_01_0226","01_05_08_01_0360","01_06_08_01_0041",...
                      "01_07_08_01_0049","01_08_08_01_0091","01_10_08_01_0143",...
                      "01_14_08_01_0372"];

overexposureFileNames = ["01_01_06_01_0147","01_01_06_01_0538","01_05_06_01_0050",...
                         "01_05_06_01_0226","01_05_06_01_0360","01_06_06_01_0041",...
                         "01_07_06_01_0049","01_08_06_01_0091","01_10_06_01_0143",...
                         "01_14_06_01_0372"];
                     
blurFileNames = ["01_01_07_01_0147","01_01_07_01_0538","01_05_07_01_0050",...
                 "01_05_07_01_0226","01_05_07_01_0360","01_06_07_01_0041",...
                 "01_07_07_01_0049","01_08_07_01_0091","01_10_07_01_0143",...
                 "01_14_07_01_0372"];
             
dirtylensFileNames = ["01_01_05_01_0147","01_01_05_01_0538","01_05_05_01_0050",...
                      "01_05_05_01_0226","01_05_05_01_0360","01_06_05_01_0041",...
                      "01_07_05_01_0049","01_08_05_01_0091","01_10_05_01_0143",...
                      "01_14_05_01_0372"];

rainFileNames = ["01_01_09_01_0147","01_01_09_01_0538","01_05_09_01_0050",...
                 "01_05_09_01_0226","01_05_09_01_0360","01_06_09_01_0041",...
                 "01_07_09_01_0049","01_08_09_01_0091","01_10_09_01_0143",...
                 "01_14_09_01_0372"];
             
shadowFileNames = ["01_01_10_01_0147","01_01_10_01_0538","01_05_10_01_0050",...
                   "01_05_10_01_0226","01_05_10_01_0360","01_06_10_01_0041",...
                   "01_07_10_01_0049","01_08_10_01_0091","01_10_10_01_0143",...
                   "01_14_10_01_0372"];
               
snowFileNames = ["01_01_11_01_0147","01_01_11_01_0538","01_05_11_01_0050",...
                 "01_05_11_01_0226","01_05_11_01_0360","01_06_11_01_0041",...
                 "01_07_11_01_0049","01_08_11_01_0091","01_10_11_01_0143",...
                 "01_14_11_01_0372"];


% Cells to store loaded images
numImages = 10;
curetsr_raw = cell(numImages,1);
curetsr_gaussian = cell(numImages,1);
curetsr_overexposure = cell(numImages,1);
curetsr_blur = cell(numImages,1);
curetsr_dirtylens = cell(numImages,1);
curetsr_rain = cell(numImages,1);
curetsr_shadow = cell(numImages,1);
curetsr_snow = cell(numImages,1);



% Load raw data paths
rawImagePath = fullfile(dataFolder,"Real_Test\ChallengeFree\");
for ii = 1:length(rawFileNames)
    imgFile = fullfile(rawImagePath,rawFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_raw{ii,1} = double(img);
end


% Load gaussian data paths
gaussianImagePath = fullfile(dataFolder,"Real_Test\Noise-1\");
for ii = 1:length(gaussianFileNames)
    imgFile = fullfile(gaussianImagePath,gaussianFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_gaussian{ii,1} = double(img);
end
             

% Load overexpsoure data paths
overexposureImagePath = fullfile(dataFolder,"Real_Test\Exposure-1\");
for ii = 1:length(overexposureFileNames)
    imgFile = fullfile(overexposureImagePath,overexposureFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_overexposure{ii,1} = double(img);
end
             

% Load blur data paths
blurImagePath = fullfile(dataFolder,"Real_Test\GaussianBlur-1\");
for ii = 1:length(blurFileNames)
    imgFile = fullfile(blurImagePath,blurFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_blur{ii,1} = double(img);
end
           

% Load dirty lens data paths
dirtylensImagePath = fullfile(dataFolder,"Real_Test\DirtyLens-1\");
for ii = 1:length(dirtylensFileNames)
    imgFile = fullfile(dirtylensImagePath,dirtylensFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_dirtylens{ii,1} = double(img);
end


% Load rain data paths
rainImagePath = fullfile(dataFolder,"Real_Test\Rain-1\");
for ii = 1:length(rainFileNames)
    imgFile = fullfile(rainImagePath,rainFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_rain{ii,1} = double(img);
end

% Load shadow data paths
shadowImagePath = fullfile(dataFolder,"Real_Test\Shadow-1\");
for ii = 1:length(shadowFileNames)
    imgFile = fullfile(shadowImagePath,shadowFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_shadow{ii,1} = double(img);
end


% Load snow data paths
snowImagePath = fullfile(dataFolder,"Real_Test\Snow-1\");
for ii = 1:length(snowFileNames)
    imgFile = fullfile(snowImagePath,snowFileNames(ii)+".bmp");
    img = imread(imgFile);
    % Compute size of image to check for RGB
    s = size(img);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img = rgb2gray(img);
    end
    % Append image to cell
    curetsr_snow{ii,1} = double(img);
end


end