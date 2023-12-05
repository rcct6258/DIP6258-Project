%% Script to Analyze Image Quality Metrics Results

% Clover Thebolt, clover.thebolt@gtri.gatech.edu
% Rohit Chakravarthy, rohit.chakravarthy@gtri.gatech.edu

clc
clear all
close all

%% General Parameters
MetricTypes = {'PSNR (dB)','SSIM Score','CWSSIM Score','UNIQUE Score',...
    'MSUNIQUE Score','CSV Score','SUMMER Score'};
EnhMethTypes = {'BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)',...
    'NLM (Spatial Method)'};
EnhMethColors = {'g','b','r'};
EnhMethMarkers = {'^','square','o'};

%% Load Set 12 Results and Initialize Parameters
% Load Results for Set12. Key:
% 1st dim is # of distortion types (10 total)
% 2nd dim is # of images (12 total)
% 3rd dim is # of metrics (7 total)
load('set12MetricsNLM_var0.001.mat','set12MetricsNLM');
load('set12MetricsBM3D_var0.001.mat','set12MetricsBM3D');
load('set12MetricsBLSGSM_var0.001.mat','set12MetricsBLSGSM');
% Combine cells into one
set12MetricsResults(:,:,:,1) = set12MetricsBLSGSM;
set12MetricsResults(:,:,:,2) = set12MetricsBM3D;
set12MetricsResults(:,:,:,3) = set12MetricsNLM;
set12MetricsResults_rev_order = set12MetricsResults;
set12MetricsResults_rev_order(6,:,:,:) = set12MetricsResults(7,:,:,:);
set12MetricsResults_rev_order(7,:,:,:) = set12MetricsResults(6,:,:,:);
set12MetricsResults = set12MetricsResults_rev_order;
% Dataset-specific Parameters
num_images_set12 = 12;
% Set12DistortionTypes = {'Gaussian (var=0.001)','Poisson','Salt & Pepper',...
%     'Speckle','Blurring','High Intensity','Low Intensity',...
%     'Under Exposure','Over Exposure','Contrast'};  % Orig order in public drive
set12DistortionTypes = {'Gaussian','Poisson','Salt & Pepper',...
    'Speckle','Blurring','Low Intensity','High Intensity',...
    'Under-Exposure','Over-Exposure','Contrast'};  % Revised order

%% (Not Used) Plot Individual Image Quality Enhancement Results
% for i_metric_type = 1:1  %length(MetricTypes)
%     for i_distortion_type = 1:1  %length(Set12DistortionTypes)  
%         figure
%         for i_enh_meth = 1:3
%             plot(1:num_images_set12,[set12MetricsResults{i_distortion_type,:,i_metric_type,i_enh_meth}],...
%                 'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
%                 'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
%             hold on
%         end
%         title('Comparison of Image Quality Enhancement Methods',...
%             strcat('Dataset=Set12, Distortion Type=',set12DistortionTypes(i_distortion_type)))
%         xlabel('Image Number'), ylabel(MetricTypes(i_metric_type))
%         legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)','NLM (Spatial Method)',...
%             'Location','best')
%     end
% end

%% (Plots Not Used) Get Enhancement Method Rankings
indiv_im_rankings = zeros(num_images_set12,3,length(set12DistortionTypes),length(MetricTypes));
all_im_avg_rankings = zeros(3,length(set12DistortionTypes),length(MetricTypes));
for i_metric_type = 1:length(MetricTypes)
    for i_distortion_type = 1:length(set12DistortionTypes)  
        for i_image = 1:num_images_set12
            [~,m_ind] = sort([set12MetricsResults{i_distortion_type,i_image,i_metric_type,:}],'descend');
            indiv_im_rankings(i_image,1:3,i_distortion_type,i_metric_type) = m_ind;
        end
        all_im_avg_rankings(:,i_distortion_type,i_metric_type) = mean(indiv_im_rankings(:,:,i_distortion_type,i_metric_type),1);
    end
end
% % Plot Rankings
% for i_metric_type = 1:1  %length(MetricTypes)
%     figure
%     for i_enh_meth = 1:3   
%         plot(1:length(set12DistortionTypes),all_im_avg_rankings(i_enh_meth,:,i_metric_type),...
%             'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
%             'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
%         hold on
%     end
%     set(gca,'Ydir','reverse')
%     title('Comparison of Image Quality Enhancement Method Rankings',...
%         strcat('Dataset=Set12, All Distortion Types, Metric Type=',MetricTypes(i_metric_type)))
%     xlabel('Distortion Type'), ylabel('Ranking (1-3)')
%     xticklabels(set12DistortionTypes)
%     legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)',...
%         'NLM (Spatial Method)','Location','best')
% end

%% (Not Used) Get Average Scores
% all_im_avg_scores = zeros(3,length(set12DistortionTypes),length(MetricTypes));
% for i_metric_type = 1:length(MetricTypes)
%     for i_distortion_type = 1:length(set12DistortionTypes)
%         for i_enh_meth = 1:3
%             all_im_avg_scores(i_enh_meth,i_distortion_type,i_metric_type) = mean([set12MetricsResults{i_distortion_type,:,i_metric_type,i_enh_meth}],2);
%         end
%     end
% end
% % Plot Average Scores Across Images
% for i_metric_type = 1:1  %length(MetricTypes)
%     figure
%     for i_enh_meth = 1:3   
%         plot(1:length(set12DistortionTypes),all_im_avg_scores(i_enh_meth,:,i_metric_type),...
%             'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
%             'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
%         hold on
%     end
%     title('Comparison of Image Quality Enhancement Method Scores',...
%         strcat('Dataset=Set12, All Distortion Types, Metric Type=',MetricTypes(i_metric_type)))
%     xlabel('Distortion Type'), ylabel('Metric Type Score')
%     xticklabels(set12DistortionTypes)
%     legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)',...
%         'NLM (Spatial Method)','Location','best')
% end

%% Normalize Scores for each image to be between 0 and 1
indiv_im_scores_norm = zeros(length(set12DistortionTypes),num_images_set12,length(MetricTypes),3);
all_im_avg_scores_norm = zeros(3,length(set12DistortionTypes),length(MetricTypes));
for i_metric_type = 1:length(MetricTypes)
    for i_distortion_type = 1:length(set12DistortionTypes)  
        for i_image = 1:12
            indiv_im_scores_norm(i_distortion_type,i_image,i_metric_type,:) = normalize([set12MetricsResults{i_distortion_type,i_image,i_metric_type,:}],"norm",Inf);
        end
        all_im_avg_scores_norm(:,i_distortion_type,i_metric_type) = mean(indiv_im_scores_norm(i_distortion_type,:,i_metric_type,:),2);
    end
end
% Plot Normalized Scores
figure
MetricType_Labels = {'PSNR','SSIM','CWSSIM','UNIQUE','MSUNIQUE','CSV','SUMMER'};
set12DistortionTypes_xaxis = repmat(set12DistortionTypes,1,length(MetricTypes));  % Revised order
for i_metric_type = 1:length(MetricTypes)
    for i_enh_meth = 1:3   
        plot((i_metric_type-1)*length(set12DistortionTypes)+1:(i_metric_type-1)*length(set12DistortionTypes)+length(set12DistortionTypes),all_im_avg_scores_norm(i_enh_meth,:,i_metric_type),...
            'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
            'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
        set(gca,'XTick',(1:length(set12DistortionTypes)*length(MetricTypes)),'XTickLabel',repmat(set12DistortionTypes,1,length(MetricTypes)))
        hold on
    end
end
title('Comparison of Image Quality Enhancement Method Scores (Normalized)',...
    strcat('Dataset=Set12, All Distortion Types, All Metric Types'))
xlabel('Distortion Type'), ylabel('Normalized Scores (0-1)')
legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)',...
    'NLM (Spatial Method)','Location','east')

%% Calc Totals of Best/#1 Rankings by Distortion Type
all_im_avg_rankings_rounded = round(all_im_avg_rankings);
enh_meths_best_totals = zeros(3,length(set12DistortionTypes));
for i_distortion_type = 1:length(set12DistortionTypes)
    for i_enh_meth = 1:3
        enh_meths_best_totals(i_enh_meth,i_distortion_type) = ...
            sum(all_im_avg_rankings_rounded(i_enh_meth,i_distortion_type,:)==1);
    end
end
% Bar graph of the totals
figure
h_bar = bar(enh_meths_best_totals');
set(gca,'XTick',(1:length(set12DistortionTypes)),'XTickLabel',set12DistortionTypes)
h_ax = gca;
h_ax.ColorOrder = [0 1 0; 0 0 1; 1 0 0];
h_text = [];
for i_bar = 1:length(h_bar)
    h_text = [h_text text(h_bar(i_bar).XData+h_bar(i_bar).XOffset,h_bar(i_bar).YData,num2str(h_bar(i_bar).YData.'), ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
end
ylim([0 7.5])
title('Image Quality Enhancement Method Rankings By Distortion Type',...
    strcat('Dataset=Set12, Summing #1 Rankings over All 7 Metric Types'))
xlabel('Distortion Type'), ylabel('#1 Rankings (max 7 possible)')
legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)','NLM (Spatial Method)',...
    'Location','north')


%% Cross-Dataset Comparison
CompData = [14 6 1; 5 12 7; 6 13 8; 5 9 10];
CompDistortionTypes = {'Gaussian','Under-Exposure','Over-Exposure','Blurring'};
figure
h_bar = bar(CompData);
set(gca,'XTick',(1:length(CompDistortionTypes)),'XTickLabel',CompDistortionTypes)
h_ax = gca;
h_ax.ColorOrder = [0 1 0; 0 0 1; 1 0 0];
h_text = [];
for i_bar = 1:length(h_bar)
    h_text = [h_text text(h_bar(i_bar).XData+h_bar(i_bar).XOffset,h_bar(i_bar).YData,num2str(h_bar(i_bar).YData.'), ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
end
ylim([0 17])
title('Image Quality Enhancement Method Rankings Across Datasets',...
    'Summing #1 Rankings over All 7 Metric Types for Selected Distortion Types')
xlabel('Distortion Type'), ylabel('#1 Rankings')
legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)','NLM (Spatial Method)',...
    'Location','north')
