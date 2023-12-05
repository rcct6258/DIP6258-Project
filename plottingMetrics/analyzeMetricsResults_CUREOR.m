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

%% Load Results and Initialize Parameters for CURE-OR 
% Key:
% 1st dim is # of distortion types (6 total)
% 2nd dim is # of images (10 total)
% 3rd dim is # of metrics (7 total)
load('CUREORMetricsNLM_var0.001.mat','CUREORMetricsNLM');
load('CUREORMetricsBM3D_var0.001.mat','CUREORMetricsBM3D');
load('CUREORMetricsBLSGSM_var0.001.mat','CUREORMetricsBLSGSM');
% Combine cells into one
CUREORMetricsResults(:,:,:,1) = CUREORMetricsBLSGSM;
CUREORMetricsResults(:,:,:,2) = CUREORMetricsBM3D;
CUREORMetricsResults(:,:,:,3) = CUREORMetricsNLM;
% Dataset-specific Parameters
num_images = 10;
CUREORDistortionTypes = {'Under-Exposure','Over-Exposure','Blurring',...
    'Contrast','Dirty Lens','Salt & Pepper'};

%% (Not Used) Plot Individual Image Quality Enhancement Results
% for i_metric_type = 1:1  %length(MetricTypes)
%     for i_distortion_type = 1:1  %length(CUREORDistortionTypes)  
%         figure
%         for i_enh_meth = 1:3
%             plot(1:num_images,[CUREORMetricsResults{i_distortion_type,:,i_metric_type,i_enh_meth}],...
%                 'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
%                 'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
%             hold on
%         end
%         title('Comparison of Image Quality Enhancement Methods',...
%             strcat('Dataset=CURE-OR, Distortion Type=',CUREORDistortionTypes(i_distortion_type)))
%         xlabel('Image Number'), ylabel(MetricTypes(i_metric_type))
%         legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)','NLM (Spatial Method)',...
%             'Location','best')
%     end
% end

%% (Plots Not Used) Get CURE-OR Enhancement Method Rankings
indiv_im_rankings = zeros(num_images,3,length(CUREORDistortionTypes),length(MetricTypes));
all_im_avg_rankings = zeros(3,length(CUREORDistortionTypes),length(MetricTypes));
for i_metric_type = 1:length(MetricTypes)
    for i_distortion_type = 1:length(CUREORDistortionTypes)  
        for i_image = 1:num_images
            [~,m_ind] = sort([CUREORMetricsResults{i_distortion_type,i_image,i_metric_type,:}],'descend');
            indiv_im_rankings(i_image,1:3,i_distortion_type,i_metric_type) = m_ind;
        end
        all_im_avg_rankings(:,i_distortion_type,i_metric_type) = mean(indiv_im_rankings(:,:,i_distortion_type,i_metric_type),1);
    end
end
% % Plot Rankings
% for i_metric_type = 1:length(MetricTypes)
%     figure
%     for i_enh_meth = 1:3   
%         plot(1:length(CUREORDistortionTypes),all_im_avg_rankings(i_enh_meth,:,i_metric_type),...
%             'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
%             'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
%         hold on
%     end
%     set(gca,'Ydir','reverse')
%     title('Comparison of Image Quality Enhancement Method Rankings',...
%         strcat('Dataset=CURE-OR, All Distortion Types, Metric Type=',MetricTypes(i_metric_type)))
%     xlabel('Distortion Type'), ylabel('Ranking (1-3)')
%     xticklabels(CUREORDistortionTypes)
%     legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)',...
%         'NLM (Spatial Method)','Location','best')
% end

%% (Not Used) Get Average Scores
% all_im_avg_scores = zeros(3,length(CUREORDistortionTypes),length(MetricTypes));
% for i_metric_type = 1:length(MetricTypes)
%     for i_distortion_type = 1:length(CUREORDistortionTypes)
%         for i_enh_meth = 1:3
%             all_im_avg_scores(i_enh_meth,i_distortion_type,i_metric_type) = mean([CUREORMetricsResults{i_distortion_type,:,i_metric_type,i_enh_meth}],2);
%         end
%     end
% end
% % Plot Average Scores
% for i_metric_type = 1:length(MetricTypes)
%     figure
%     for i_enh_meth = 1:3   
%         plot(1:length(CUREORDistortionTypes),all_im_avg_scores(i_enh_meth,:,i_metric_type),...
%             'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
%             'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
%         hold on
%     end
%     title('Comparison of Image Quality Enhancement Method Scores',...
%         strcat('Dataset=CURE-OR, All Distortion Types, Metric Type=',MetricTypes(i_metric_type)))
%     xlabel('Distortion Type'), ylabel('Metric Type Score')
%     xticklabels(CUREORDistortionTypes)
%     legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)',...
%         'NLM (Spatial Method)','Location','best')
% end

%% Normalize Scores for each image to be between 0 and 1
indiv_im_scores_norm = zeros(length(CUREORDistortionTypes),num_images,length(MetricTypes),3);
all_im_avg_scores_norm = zeros(3,length(CUREORDistortionTypes),length(MetricTypes));
for i_metric_type = 1:length(MetricTypes)
    for i_distortion_type = 1:length(CUREORDistortionTypes)  
        for i_image = 1:num_images
            indiv_im_scores_norm(i_distortion_type,i_image,i_metric_type,:) = normalize([CUREORMetricsResults{i_distortion_type,i_image,i_metric_type,:}],"norm",Inf);
        end
        all_im_avg_scores_norm(:,i_distortion_type,i_metric_type) = mean(indiv_im_scores_norm(i_distortion_type,:,i_metric_type,:),2);
    end
end
% Plot Normalized Scores
figure
MetricType_Labels = {'PSNR','SSIM','CWSSIM','UNIQUE','MSUNIQUE','CSV','SUMMER'};
for i_metric_type = 1:length(MetricTypes)
    for i_enh_meth = 1:3   
        plot((i_metric_type-1)*length(CUREORDistortionTypes)+1:(i_metric_type-1)*length(CUREORDistortionTypes)+length(CUREORDistortionTypes),all_im_avg_scores_norm(i_enh_meth,:,i_metric_type),...
            'Color',EnhMethColors{i_enh_meth},'LineStyle','none',...
            'Marker',EnhMethMarkers{i_enh_meth},'MarkerSize',8,'MarkerFaceColor',EnhMethColors{i_enh_meth})
        set(gca,'XTick',(1:length(CUREORDistortionTypes)*length(MetricTypes)),'XTickLabel',repmat(CUREORDistortionTypes,1,length(MetricTypes)))
        hold on
    end
end
title('Comparison of Image Quality Enhancement Method Scores (Normalized)',...
    strcat('Dataset=CURE-OR, All Distortion Types, All Metric Types'))
xlabel('Distortion Type'), ylabel('Normalized Scores (0-1)')
legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)',...
    'NLM (Spatial Method)','Location','east')

%% Calc Totals of Best/#1 Rankings by Distortion Type
all_im_avg_rankings_rounded = round(all_im_avg_rankings);
enh_meths_best_totals = zeros(3,length(CUREORDistortionTypes));
for i_distortion_type = 1:length(CUREORDistortionTypes)
    for i_enh_meth = 1:3
        enh_meths_best_totals(i_enh_meth,i_distortion_type) = ...
            sum(all_im_avg_rankings_rounded(i_enh_meth,i_distortion_type,:)==1);
    end
end
% Bar graph of the totals
figure
h_bar = bar(enh_meths_best_totals');
set(gca,'XTick',(1:length(CUREORDistortionTypes)),'XTickLabel',CUREORDistortionTypes)
h_ax = gca;
h_ax.ColorOrder = [0 1 0; 0 0 1; 1 0 0];
h_text = [];
for i_bar = 1:length(h_bar)
    h_text = [h_text text(h_bar(i_bar).XData+h_bar(i_bar).XOffset,h_bar(i_bar).YData,num2str(h_bar(i_bar).YData.'), ...
                          'VerticalAlignment','bottom','horizontalalign','center')];
end
ylim([0 7.5])
title('Image Quality Enhancement Method Rankings By Distortion Type',...
    strcat('Dataset=CURE-OR, Summing #1 Rankings over All 7 Metric Types'))
xlabel('Distortion Type'), ylabel('#1 Rankings (max 7 possible)')
legend('BLS-GSM (Wavelet-Transform Method)','BM3D (Blended Method)','NLM (Spatial Method)',...
    'Location','north')


