classdef DIPMetrics < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        psnr_;
        ssim_;
        cwssim_
        unique_
        msunique_;
        csv_;
        summer_;
        
        
    end
    %%
    methods
        %% Compute PSNR
        function psnr_ = PSNR(obj,inputImage,refImage)
            [peaksnr, snr_] = psnr(inputImage,refImage);
            psnr_ = peaksnr;
        end
        
        %% Compute SSIM
        
        function ssim_ = SSIM(obj,inputImage,refImage)
            [ssimval, ssimmap] = ssim(inputImage,refImage);
            ssim_ = ssimval;
        end
        
        %% Compute CW-SSIM
        
        function cwssim_ = CWSSIM(obj,img1, img2, level, or, guardb, K)
            addpath("matlabPyrTools");
            %========================================================================
            %CW-SSIM Index, Version 1.0
            %Copyright(c) 2010  Zhou Wang, Mehul Sampat and Alan Bovik
            %All Rights Reserved.
            %----------------------------------------------------------------------
            %Permission to use, copy, or modify this software and its documentation
            %for educational and research purposes only and without fee is hereby
            %granted, provided that this copyright notice and the original authors'
            %names appear on all copies and supporting documentation. This program
            %shall not be used, rewritten, or adapted as the basis of a commercial
            %software or hardware product without first obtaining permission of the
            %authors. The authors make no representations about the suitability of
            %this software for any purpose. It is provided "as is" without express
            %or implied warranty.
            %----------------------------------------------------------------------
            %
            %This is an implementation of the algorithm for calculating the Complex-Wavelet 
            %Structural SIMilarity (CW-SSIM) index between two images. Please refer
            %to the following paper:
            %
            % M. P. Sampat, Z. Wang, S. Gupta, A. C. Bovik, M. K. Markey. 
            % "Complex Wavelet Structural Similarity: A New Image Similarity Index", 
            % IEEE Transactions on Image Processing, 18(11), 2385-401, 2009.
            %
            % *** Important code dependencies: ***
            % This code requires the "matlabPyrTools" package developed by Prof. Eero Simoncelli. 
            % This package can be downloaded from: http://www.cns.nyu.edu/~lcv/software.php 
            %
            % Kindly report any suggestions or corrections to mehul.sampat@ieee.org
            %
            %----------------------------------------------------------------------
            %
            % Input :::
            %            (1) img1......the first image being compared
            %            (2) img2......the second image being compared
            %            (3) level......the number of levels to used in the complex steerable pyramid decomposition   
            %            (4) or..........the number of orientations to be used in the complex steerable pyramid decomposition     
            %            (5) guardb...this parameter is used to control how much is discarded from the four image boundaries. 
            %            (6) K...........the constant in the CWSSIM index formula (see the above reference) default value: K=0
            %
            % Output ::: 
            %            (1) cwssim...the CWSSIM index value between 2 images. If img1 = img2, then cwssim = 1.
            %
            % Example Usage: Given 2 test images img1 and img2
            %
            % cwssim = cwssim_index(img1, img2,1,16,0,0);
            %
            % See the results: "cwssim" gives the CW-SSIM index value
            %========================================================================

            [pyr1, pind] = buildSCFpyr(img1, level, or-1);%........% decompose img1 using a complex steerable pyramid decomposition
            [pyr2, pind] = buildSCFpyr(img2, level, or-1);%........% decompose img2 using a complex steerable pyramid decomposition

            winsize = 7;
            window = ones(7);%..............................................% The CW-SSIM indices are computed locally using a sliding 
            %                                                                         % 7-by-7 window that moves across each wavelet subband.

            window = window./sum(sum(window));%................% normalize the window

            gb = guardb/(2^(level-1));%...................................% the gb parameter is used to control how much is discarded from the four image boundaries. 

            s = pind((level-1)*or+2, :);
            w = fspecial('gaussian', s-winsize+1, s(1)/4);%........% The CW-SSIM index map is combined into a scalar similarity measure using a
                                                                                       % weighted summation.The weighting function is obtained using a Gaussian 
                                                                                       % profile with a standard deviation equaling a quarter of the image size at 
                                                                                       % finest level of pyramid decomposition.

            for i=1:or
               bandind = i+(level-1)*or+1;
               band1 = pyrBand(pyr1, pind, bandind);%............% Access a subband from a pyramid (see help pyrBand)
               band2 = pyrBand(pyr2, pind, bandind);%............% Access a subband from a pyramid (see help pyrBand)
               band1 = band1(gb+1:end-gb, gb+1:end-gb);
               band2 = band2(gb+1:end-gb, gb+1:end-gb);

               corr = band1.*conj(band2);
               varr = abs(band1).^2 + abs(band2).^2;
               corr_band = filter2(window, corr, 'valid');%.........% The CW-SSIM indices are computed locally using a sliding 7-by-7 
               varr_band = filter2(window, varr, 'valid');%.........% window that moves across each wavelet subband.

               cssim_map = ...
               (2*abs(corr_band) + K)./(varr_band + K);%........% The purpose of the small constant K is mainly to improve the robustness of 
                                                                                       % the CW-SSIM measure when the local signal to noise ratios are low.

               band_cssim(i) = sum(sum(cssim_map.*w));%......% The CW-SSIM index map is combined into a scalar similarity measure using a
                                                                                       % weighted summation.
            end

            cwssim_ = mean(band_cssim);%...............................% value is 1 for two identical images.
        end
        
        %% Compute UNIQUE
        
        function unique_ = UNIQUE(obj,inputImage,refImage)
            unique_ = mslUNIQUE(obj,inputImage,refImage);
        end
        
        function result = mslUNIQUE(obj,img1,img2)
            
            %  Author:              Mohit Prabhushankar
            %  PI:                  Ghassan AlRegib
            %  Version:             1.0
            %  Published in:        Signal Processing Letter October 2016
            %  Publication details: 

            
            %Resize image to 0.5 times the original size - Faster processing
            %and HVS is more adapted to low frequency components
            image1 = imresize(img1,0.5);
            image2 = imresize(img2,0.5);
            
            if length(size(image1)) == 2
                image1 = cat(3,image1,image1,image1);
            end
            if length(size(image2)) == 2
                image2 = cat(3,image2,image2,image2);
            end
                
            %Loading Precalculated weights and bias
            workspace = load('InputWeights/ImageNet_Weights_YGCr.mat');        
            weight = workspace.W;  
            bias = workspace.b;

            %ColorSpace transformation
            img1 = rgb2ycbcr(image1);
            img2 = rgb2ycbcr(image2);
            img1(:,:,2) = image1(:,:,2);
            img2(:,:,2) = image2(:,:,2);     
            

            %Preparing images (Zero centering and ZCA whitening) and
            %multiplying by weights and adding bias
            img1_s = mslProcessUNIQUE(obj,img1,weight,bias);
            img2_s = mslProcessUNIQUE(obj,img2,weight,bias);

            %Discount those features that are much lesser than average
            %activation(0.035) - Suppression
            res = find(img1_s < 0.025);
            if (sum(res) > 0)              
                img1_s(res) = 0;              
            end

            res = find(img2_s < 0.025);
            if (sum(res) > 0)              
                img2_s(res) = 0;              
            end  

            %Pooling using 10th power of Spearman Correlation coefficient
            result = abs(corr(img1_s,img2_s,'type','Spearman'))^10;
        end
        
        function feature = mslProcessUNIQUE(obj,img,W,b)
            
            %  Author:              Mohit Prabhushankar
            %  PI:                  Ghassan AlRegib
            %  Version:             1.0
            %  Published in:        Signal Processing Letter October 2016
            %  Publication details: 
            
            I = im2double(img);

            %Parameter Initialisation
            [m,n,~] = size(I);
            epsilon = 0.1; 
            count = 1; 
            scale = 8;

            %Convert m x n x 3 image into [(8x8x3) x count] patches
            i = 1;
            while (i < m - (scale - 2))
                j = 1;
                while (j< n-(scale-2)) %(j < 512)
                    patch_temp = I(i:i+(scale-1),j:j+(scale-1),:);
                    patches(:,count) = reshape(patch_temp,[],1);
                    count = count+1;
                    j = j+scale;
                end    
                i = i+scale;
            end

            % Subtract mean patch (hence zeroing the mean of the patches)
            meanPatch = mean(patches,2);  
            patches = bsxfun(@minus, patches, meanPatch);

            % Apply ZCA whitening
            sigma = patches * patches' / (count-1);
            [u, s, ~] = svd(sigma);
            ZCAWhite = u * diag(1 ./ sqrt(diag(s) + epsilon)) * u';
            patches = ZCAWhite * patches;

            %Process the patches using the different models and multiply 
            %resultant by the sharpness indices
            feature = mslComputeUNIQUE(obj,patches,W,b);

            %Reshaping back to a single vector
            feature = reshape(feature,[],1);
        end
        
        function feature = mslComputeUNIQUE(obj,patches,weight,bias)
            feature = weight * patches + repmat(bias,1,size(patches,2));
            feature = 1./(1 + exp(-(feature)));    
        end
        %% Compute MS-UNIQUE
        
        function msunique_ = MSUNIQUE(obj,inputImage,refImage)
            msunique_ = mslMSUNIQUE(obj,inputImage,refImage);
        end
        
        function result = mslMSUNIQUE(obj,img1,img2)
            
            %  Author:              Mohit Prabhushankar
            %  PI:                  Ghassan AlRegib
            %  Version:             1.0
            %  Published in:        Electronic Imaging 2017
            %  Publication details: 

            
            %Resize image to 0.5 times the original size - Faster processing
            %and HVS is more adapted to low frequency components
            image1 = imresize(img1,0.5);
            image2 = imresize(img2,0.5);
            
            if length(size(image1)) == 2
                image1 = cat(3,image1,image1,image1);
            end
            if length(size(image2)) == 2
                image2 = cat(3,image2,image2,image2);
            end

            %Loading Precalculated weights and bias
            [W1,b1,Ind1] = mslLoadVariables(obj,'InputWeights/ImageNet_Weights_YGCr');
            [W2,b2,Ind2] = mslLoadVariables(obj,'InputWeights/ImageNet_Weights_YGCr_625x192');
            [W3,b3,Ind3] = mslLoadVariables(obj,'InputWeights/ImageNet_Weights_YGCr_81x192');
            [W4,b4,Ind4] = mslLoadVariables(obj,'InputWeights/ImageNet_Weights_YGCr_121x192');
            [W5,b5,Ind5] = mslLoadVariables(obj,'InputWeights/ImageNet_Weights_YGCr_169x192');

            %ColorSpace transformation
            img1 = rgb2ycbcr(image1);
            img2 = rgb2ycbcr(image2);
            img1(:,:,2) = image1(:,:,2);
            img2(:,:,2) = image2(:,:,2);     

            %Preparing images (Zero centering and ZCA whitening) and
            %multiplying by weights and adding bias
            img1_s = mslProcessMSUNIQUE(obj,img1,W1,b1,Ind1,W2,b2,Ind2,W3,b3,Ind3,W4,b4,Ind4,W5,b5,Ind5);
            img2_s = mslProcessMSUNIQUE(obj,img2,W1,b1,Ind1,W2,b2,Ind2,W3,b3,Ind3,W4,b4,Ind4,W5,b5,Ind5);

            %Discount those features that are much lesser than average
            %activation(0.035) - Suppression
            res = find(img1_s < 0.025);
            if (sum(res) > 0)              
                img1_s(res) = 0;              
            end

            res = find(img2_s < 0.025);
            if (sum(res) > 0)              
                img2_s(res) = 0;              
            end  

            %Pooling using 10th power of Spearman Correlation coefficient
            result = abs(corr(img1_s,img2_s,'type','Spearman'))^10;
        end  
        
        function feature = mslProcessMSUNIQUE(obj,img,W1,b1,Ind1,W2,b2,Ind2,W3,b3,Ind3,W4,b4,Ind4,W5,b5,Ind5)
            
            %  Author:              Mohit Prabhushankar
            %  PI:                  Ghassan AlRegib
            %  Version:             1.0
            %  Published in:        
            %  Publication details: 
            
            I = im2double(img);

            %Parameter Initialisation
            [m,n,~] = size(I);
            epsilon = 0.1; 
            count = 1; 
            scale = 8;

            %Convert m x n x 3 image into [(8x8x3) x count] patches
            i = 1;
            while (i < m - (scale - 2))
                j = 1;
                while (j< n-(scale-2)) %(j < 512)
                    patch_temp = I(i:i+(scale-1),j:j+(scale-1),:);
                    patches(:,count) = reshape(patch_temp,[],1);
                    count = count+1;
                    j = j+scale;
                end    
                i = i+scale;
            end

            % Subtract mean patch (hence zeroing the mean of the patches)
            meanPatch = mean(patches,2);  
            patches = bsxfun(@minus, patches, meanPatch);

            % Apply ZCA whitening
            sigma = patches * patches' / (count-1);
            [u, s, ~] = svd(sigma);
            ZCAWhite = u * diag(1 ./ sqrt(diag(s) + epsilon)) * u';
            patches = ZCAWhite * patches;

            %Process the patches using the different models and multiply 
            %resultant by the sharpness indices
            feature1 = mslComputeMSUNIQUE(obj,patches,W1,b1,Ind1);
            feature2 = mslComputeMSUNIQUE(obj,patches,W2,b2,Ind2);
            feature3 = mslComputeMSUNIQUE(obj,patches,W3,b3,Ind3);
            feature4 = mslComputeMSUNIQUE(obj,patches,W4,b4,Ind4);
            feature5 = mslComputeMSUNIQUE(obj,patches,W5,b5,Ind5);

            feature_full = [feature1;feature2;feature3;feature4;feature5];

            %Reshaping back to a single vector
            feature = reshape(feature_full,[],1);
        end
        
        function feature = mslComputeMSUNIQUE(obj,patches,weight,bias,Indices)
            feature = weight * patches + repmat(bias,1,size(patches,2));
            feature = 1./(1 + exp(-(feature)));    
            feature = feature.*(repmat(Indices,1,size(patches,2)));
        end
        
        function [weight,bias,Indices] = mslLoadVariables(obj,var)

        % Load workspace
            workspace = load([var,'.mat']);        
            weight = workspace.W;
            bias = workspace.b;
            numFilters = length(bias);

        % Separate filters and weigh them according to sharpness

            weight_process = weight(:,1:64);
            weight_process = weight_process - (sum(weight_process(:))/length(weight_process(:)));
            weight_process = weight_process./((ones(size(weight_process,1),1))*max(abs(weight_process)));

            kurt = kurtosis(weight_process,0,2);

            Ind = find((kurt > 5));
            Indices = ones(numFilters,1);
            Indices(Ind(1:length(Ind)),1) = 2;

            Ind2 = find((kurt < 2));
            Indices(Ind2(1:length(Ind2)),1) = 0.5;
        end
        %% Compute CSV
        
        function csv_ = CSV(obj,img1, img2)
            
            addpath(genpath("csvFuncs"));
            
            if length(size(img1)) == 2
                img1 = cat(3,img1,img1,img1);
            end
            if length(size(img2)) == 2
                img2 = cat(3,img2,img2,img2);
            end
            
            %Perceptual color weights (d in Eq. 4)Section 3.2, Figure 3
            D= [0 1 0.9497 0.76544 1 1 1 1 1 1 1;...
            1 0 1 1 1 1 1 1 1 1 1;...
            0.9497 1 0 0.93538 1 1 1 1 1 1 1;...
            0.76544 1 0.93538 0 1 1 1 1 1 0.68809 1;...
            1 1 1 1 0 1 1 1 1 1 1;...
            1 1 1 1 1 0 1 1 0.92114 1 1;...
            1 1 1 1 1 1 0 1 1 1 1;...
            1 1 1 1 1 1 1 0 1 1 1;...
            1 1 1 1 1 0.92114 1 1 0 1 1;...
            1 1 1 0.68809 1 1 1 1 1 0 1;...
            1 1 1 1 1 1 1 1 1 1 0];
            T=20;

            % Section 6.3 - Parameters 
            % blockSize is the W in the paper
            blockSize=[20 20];
            % weightCND is the A in the paper 
            weightCND=0.9;
            % P is the power of the mean value
            P=4;
            % stdVal is the \sigma
            stdVal=50;
            % K_L, K_C and K_H are 
            K_L=1;
            K_C=1;
            K_H=1;
            K=[K_L,K_C,K_H];
            % n is not explicitly defined since 11 word dictionary is directly used in
            % RGB_to_color_terms function.




            % Section 5. - Retinal Ganglion cell-based difference (RGCD)
            %Define filter type as Laplacian of Gaussian (LoG) 
            filterType='log';
            %Initialize standard deviation (mentioned in Section 6.3.)

            %Generate filter
            h = fspecial(filterType, blockSize,stdVal);
            %Filter the R channels of compared images 
            imgInd=1;
            img1LoG_R=imfilter(double(img1(:,:,imgInd)),h,'replicate');
            img2LoG_R=imfilter(double(img2(:,:,imgInd)),h,'replicate');
            %Filter the G channels of compared images 
            imgInd=2;
            img1LoG_G=imfilter(double(img1(:,:,imgInd)),h,'replicate');
            img2LoG_G=imfilter(double(img2(:,:,imgInd)),h,'replicate');
            %Filter the B channels of compared images 
            imgInd=3;
            img1LoG_B=imfilter(double(img1(:,:,imgInd)),h,'replicate');
            img2LoG_B=imfilter(double(img2(:,:,imgInd)),h,'replicate');
            % Calculate the difference between filter maps
            imgLoG_R=abs(img1LoG_R-img2LoG_R);
            imgLoG_G=abs(img1LoG_G-img2LoG_G);
            imgLoG_B=abs(img1LoG_B-img2LoG_B);
            % Pool the individual maps using geometric mean to obtain the RGCD map
            csv_rgcd=nthroot(imgLoG_R.*imgLoG_G.*imgLoG_B,3);

            % Section 4 - Blockwise Structural Difference  (BSD)
            %Define local normalization operator 
            funNorm=@(block_struct) (block_struct.data-mean2(block_struct.data))/...
            (std(reshape(block_struct.data,[numel(block_struct.data) 1]))+0.001);
            %Calculate BSD in each color channel separately
            imgBSD_R=abs(blockproc(double(img1(:,:,1)),blockSize,funNorm)-blockproc(double(img2(:,:,1)),blockSize,funNorm));
            imgBSD_G=abs(blockproc(double(img1(:,:,2)),blockSize,funNorm)-blockproc(double(img2(:,:,2)),blockSize,funNorm));
            imgBSD_B=abs(blockproc(double(img1(:,:,3)),blockSize,funNorm)-blockproc(double(img2(:,:,3)),blockSize,funNorm));
            % Pool the individual maps using geometric mean to obtain the BSD map
            csv_bsd=nthroot(imgBSD_R.*imgBSD_G.*imgBSD_B,3);

            % Section 3 Perceptual Color Difference
            %Define mean pooling operation
            colorConv=makecform('srgb2lab');




            funMean=@(block_struct) (mean2(block_struct.data));




            % Mean pool the color channels of compared images
            for ii=1:3    
            img1Mean(:,:,ii)=blockproc(double(img1(:,:,ii)),blockSize,funMean);
            img2Mean(:,:,ii)=blockproc(double(img2(:,:,ii)),blockSize,funMean);
            end

            lab_1= applycform(double(squeeze(img1Mean))./255,colorConv );
            lab_2= applycform(double(squeeze(img2Mean))./255,colorConv);


            % Initilize CIEDE and CND maps
            [s1,s2,~]=size(img1Mean);
            ciedeMap=zeros(s1,s2);
            cndMap=ciedeMap;
            % Preprocess data to avoid out of range issues and calculate CIEDE and CND
            % values for each pixel
            for i=1:s1
            for j=1:s2          
            ciedeMap(i,j)= min([deltaE2000(squeeze(lab_1(i,j,:))',squeeze(lab_2(i,j,:))',K), T])/T;
            cn_1= RGB_to_color_terms(squeeze(img1Mean(i,j,:))');
            cn_2= RGB_to_color_terms(squeeze(img2Mean(i,j,:))');
            cndMap(i,j)= emd_hat_gd_metric_mex(cn_1',cn_2',D);   
            end    
            end
            % Interpolate CIEDE and CND maps to input image resolution
            % Mean is subtracted from the resized maps since interpolation can lead to
            % negative values.
            ciedeMap=imresize(ciedeMap,size(imgBSD_R));
            ciedeMap=ciedeMap-min(min(ciedeMap));
            cndMap=imresize(cndMap,size(imgBSD_R));
            cndMap=cndMap-min(min(cndMap));

            % Section 6.1. -  PCD map is calculated by combining ciede and cnd maps
            csv_pcd=weightCND*cndMap+(1-weightCND)*ciedeMap;

            % PCD,RGCD and BSD maps are multiplicatively combined
            csvMap=csv_bsd.*csv_rgcd.*csv_pcd;
            % Section 6.1. - CSV is obtained from the pooled map as in Eq 14.
            csv_=(1-nthroot(mean2(csvMap),P));
        end
        
          
        
        %% Compute SUMMER
        function summer_ = SUMMER(obj,im1,im2)
            %
            %  Author:              Dogancan Temel
            %  PI:                  Ghassan AlRegib
            %  Version:             1.0
            %  Date:                August 8, 2018
            %  Submitted to:        Elsevier Signal Processing: Image Communication 
            %
            %  Input:               Reference and distorted images (color)
            %                       im1: reference image
            %                       im2: distorted image
            %  Output:              SUMMER score
            
            
            if length(size(im1)) == 2
                im1 = cat(3,im1,im1,im1);
            end
            if length(size(im2)) == 2
                im2 = cat(3,im2,im2,im2);
            end
            
            
            outCT=0;
            for c=1:3      
                inp = double(im1(:,:,c))/255- double(im2(:,:,c))/255;
                for j = 1:4
                    inpR = imresize(inp, 1/2^j,'box');
                    out=mean2(abs(log(1+abs(fft2(abs(inpR)))))) ;
                    outCT = outCT + out;
                end
            end          
            for c=1:3
                for j=3:4
                    x1 = abs(fft2(imresize(double(im1(:,:,c))/255, 1/2^(j),'bilinear')));
                    x2 = abs(fft2(imresize(double(im2(:,:,c))/255, 1/2^(j),'bilinear')));
                    outCT = outCT*log(1+mean2(x1./x2)) ;

                end
            end
            summer_=5*nthroot(1./(outCT+1),3);
        end
        
        
                   
        
        %% Compute All Metrics
        
        % Populate properties with calculated metrics
        function calculateAllMetrics(obj,inputImage,refImage)
            
            obj.psnr_ = obj.PSNR(inputImage,refImage);
            obj.ssim_ = obj.SSIM(inputImage,refImage);
            obj.cwssim_ = obj.CWSSIM(reshape(inputImage,[],size(inputImage,2)),reshape(refImage,[],size(refImage,2)),1,16,0,0);
            obj.unique_ = obj.UNIQUE(inputImage,refImage);
            obj.msunique_ = obj.MSUNIQUE(inputImage,refImage);
            obj.csv_ = obj.CSV(inputImage,refImage);
            obj.summer_ = obj.SUMMER(inputImage,refImage);
            
            
        end
        
        function resetMetrics(obj)
            obj.psnr_ = [];
            obj.ssim_ = [];
            obj.cwssim_ = [];
            obj.unique_ = [];
            obj.msunique_ = [];
            obj.csv_ = [];
            obj.summer_ = [];
            
        end
        
    end
end

