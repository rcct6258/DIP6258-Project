function [sidd_raw_gaussian, sidd_raw_underexposure,sidd_gaussian, sidd_underexposure] = loadSIDD(dataFolder)

% 1 - Noisy Raw-RGB image (.MAT).
% Black Level subtracted, normalized to [0, 1].
% 
% 2 - Ground truth Raw-RGB image (.MAT).
% Black Level subtracted, normalized to [0, 1].

% L - low light
% N - normal brightness
% H - high brightness
% 
% iphoneFolders_L = dir("..\data\SIDD\SIDD_Small_Raw_Only\Data\*IP*L");
% iphoneFolders_N = dir("..\data\SIDD\SIDD_Small_Raw_Only\Data\*IP*N");

iphoneFolders_L = flipud(dir(fullfile(dataFolder+"\Data","*IP*L")));
iphoneFolders_N = flipud(dir(fullfile(dataFolder+"\Data","*IP*N")));

numImagesL = 5;%length(iphoneFolders_L);
numImagesN = 5;%length(iphoneFolders_N);


sidd_raw_gaussian = cell(numImagesN,1);
sidd_gaussian = cell(numImagesN,1);
sidd_raw_underexposure = cell(numImagesL,1);
sidd_underexposure = cell(numImagesL,1);


% Load noisy gaussian data
for ii = 1:numImagesN
    imageFolder = fullfile(iphoneFolders_N(ii).folder,iphoneFolders_N(ii).name);
    img_raw_path = fullfile(imageFolder,"GT_RAW_010.mat");
    img_gaussian_path = fullfile(imageFolder,"NOISY_RAW_010.mat");
    
    img_raw = load(img_raw_path);
    s = size(img_raw);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img_raw = rgb2gray(img_raw);
    end
    
    img_gaussian = load(img_gaussian_path);
    s = size(img_gaussian);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img_gaussian = rgb2gray(img_gaussian);
    end
    
    % Append images to cells
    sidd_raw_gaussian{ii,1} = double(img_raw.x);
    sidd_gaussian{ii,1} = double(img_gaussian.x);
end






% Load noisy underexposure data
for ii = 1:numImagesL
    imageFolder = fullfile(iphoneFolders_L(ii).folder,iphoneFolders_L(ii).name);
    img_raw_path = fullfile(imageFolder,"GT_RAW_010.mat");
    img_underexposure_path = fullfile(imageFolder,"NOISY_RAW_010.mat");
    
    img_raw = load(img_raw_path);
    s = size(img_raw);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img_raw = rgb2gray(img_raw);
    end
    
    img_underexposure = load(img_underexposure_path);
    s = size(img_underexposure);
    % Convert three channel RGB image to grayscale if condition applies
    if ismember(3,s)
        img_underexposure = rgb2gray(img_underexposure);
    end
    
    % Append images to cells
    sidd_raw_underexposure{ii,1} = double(img_raw.x);
    sidd_underexposure{ii,1} = double(img_underexposure.x);
end


