
version = 19;
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Figures/';
AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
AllSetsVersionCombinedEmbryosPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'CompiledEmbryos/'];
AllSetsCombinedEmbryosPath = ['S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/'];
load(SizeDataPath,'APLengths','DVLengths','VentralDistances','DorsalDistances','NGoodEmbryos',...
    'MeanAPLength','StdAPLength','MeanDVLength','StdDVLength','AspectRatios');
MeanAspectRatio = mean(AspectRatios);
StdAspectRatio = std(AspectRatios);
[DubuisTimes, DubuisMeanProfile, ~] = getMembraneFurrowProfiles( 'dubuis');
[yw25CTimes, yw25CProfile, yw25CSE] = getMembraneFurrowProfiles( 'yw25csquished_nopv');
[hisrfp25CTimes, hisrfp25CProfile, hisrfp25CSE] = getMembraneFurrowProfiles( 'hisrfp25c_nopv');
AllSetInfo = GetFixedSetPrefixInfo;

NumSets = length(AllSetInfo.Temperatures);

NChannels = 5;



APbins = 0:0.025:1;
NumAPbins = length(APbins);

NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);
NumAPbins = length(APbins);

ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32;
    76 0 153;
    255 0 255]/255;%0.9290 0.6940 0.1250]/255;

rainbow_colors = [212 84 66;
    255 128 0;
    237 176 32;
    100 149 93;
    66 87 168;
    76 0 153;
    255 0 255]/255;
    
SetLineStyles = {'-', '--', ':'};
%
%load([AllCompiledPath, 'AllCompiledEmbryos.Mat'], 'AllCompiledEmbryos')
%%

master_index = 7;
for temp = [25, 27.5, 22.5, 20, 17.5]
    for ch_index = 5%[3 5]
        exp_indices = find(AllSetInfo.Temperatures == temp);
        exp_index = exp_indices(1);
        exp_index2 = exp_indices(2);
        exp_index3 = exp_indices(3);
        
        x_sample = {};
        UseTF = {};
        ys = {};
        BadTF = {};
        xfits = {};
        yfits = {};
        x_sample_control = {};
        UseTF_control = {};
        ys_control = {};
        BadTF_control = {};
        xfits_control = {};
        yfits_control = {};
        SetLabel = {};
        SlideIDs = {};
        SlideIDs_control = {};
        
        CompiledEmbryos = AllCompiledEmbryos{exp_index};
        x_sample{1} = CompiledEmbryos.DubuisEmbryoTimes;    
        SlideIDs{1} = CompiledEmbryos.SlideIDs;
        UseTF{1} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1});
        x_sample{1} =x_sample{1}(UseTF{1});
        SlideIDs{1}  = SlideIDs{1}(UseTF{1});
        %ys{1} =  CompiledEmbryos.UnivScaledProfiles.Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep1(UseTF{1},:,ch_index, master_index);
        ys{1} =  CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{1},:,ch_index, master_index);
        BadTF{1} =  sum(isnan(ys{1}), 2).' >min( sum(isnan(ys{1}), 2).');
        x_sample{1} = x_sample{1}(~BadTF{1});
        SlideIDs{1} = SlideIDs{1}(~BadTF{1});
        ys{1} = ys{1}(~BadTF{1},:);
        %xfits{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep1.x;
        %yfits{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep1.TestSet.mean(:,:,ch_index, master_index);
        yfits{1} =CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        xfits{1} =CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        
        x_sample_control{1} = CompiledEmbryos.DubuisEmbryoTimes;
        SlideIDs_control{1} = CompiledEmbryos.SlideIDs;
        UseTF_control{1} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1});
        x_sample_control{1} =x_sample_control{1}(UseTF_control{1});
        SlideIDs_control{1}  = SlideIDs_control{1}(UseTF_control{1});
        %ys_control{1} =  CompiledEmbryos.UnivScaledProfiles.Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep1(UseTF_control{1},:,ch_index, master_index);
        ys_control{1} = CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{1},:,ch_index, master_index);
        BadTF_control{1} =  sum(isnan(ys_control{1}), 2).' >min( sum(isnan(ys_control{1}), 2).');
        x_sample_control{1} = x_sample_control{1}(~BadTF_control{1});
        SlideIDs_control{1} = SlideIDs_control{1}(~BadTF_control{1});
        ys_control{1} = ys_control{1}(~BadTF_control{1},:);
        %xfits_control{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep1.x;
        %yfits_control{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep1.ControlSet.mean(:,:,ch_index, master_index);
        yfits_control{1} =CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        xfits_control{1} =CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        
        SetLabel{1} = AllSetInfo.SetLabels{exp_index};
        
        CompiledEmbryos2 = AllCompiledEmbryos{exp_index2};
        x_sample{2} = CompiledEmbryos2.DubuisEmbryoTimes;
        SlideIDs{2} = CompiledEmbryos2.SlideIDs;
        UseTF{2} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2});
        x_sample{2} =x_sample{2}(UseTF{2});
        SlideIDs{2} =SlideIDs{2}(UseTF{2});
        %ys{2} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF{2},:,ch_index, master_index);
        ys{2} =  CompiledEmbryos2.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{2},:,ch_index, master_index);
        BadTF{2} =  sum(isnan(ys{2}), 2).' >min( sum(isnan(ys{2}), 2).');
        x_sample{2} = x_sample{2}(~BadTF{2});
        SlideIDs{2} = SlideIDs{2}(~BadTF{2});
        ys{2} = ys{2}(~BadTF{2},:);
        %xfits{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep2.x;
        %yfits{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:,:,ch_index, master_index);
        yfits{2} =CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        xfits{2} =CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        
        
        x_sample_control{2} = CompiledEmbryos2.DubuisEmbryoTimes;
        SlideIDs_control{2} = CompiledEmbryos2.SlideIDs;
        UseTF_control{2} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2});
        x_sample_control{2} =x_sample_control{2}(UseTF_control{2});
        SlideIDs_control{2} =SlideIDs_control{2}(UseTF_control{2});
        %ys_control{2} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control{2},:,ch_index, master_index);
        ys_control{2} =  CompiledEmbryos2.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{2},:,ch_index, master_index);
        BadTF_control{2}=  sum(isnan(ys_control{2}), 2).' >min( sum(isnan(ys_control{2}), 2).');
        x_sample_control{2} = x_sample_control{2}(~BadTF_control{2});
        SlideIDs_control{2} = SlideIDs_control{2}(~BadTF_control{2});
        ys_control{2} = ys_control{2}(~BadTF_control{2},:);
        %xfits_control{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep2.x;
        %yfits_control{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean(:,:,ch_index, master_index);
        yfits_control{2} =CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        xfits_control{2} =CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        
        SetLabel{2} = AllSetInfo.SetLabels{exp_index2};
        
        CompiledEmbryos3 = AllCompiledEmbryos{exp_index3};
        x_sample{3} = CompiledEmbryos3.DubuisEmbryoTimes;
        SlideIDs{3} = CompiledEmbryos3.SlideIDs;
        UseTF{3} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3});
        x_sample{3} =x_sample{3}(UseTF{3});
        SlideIDs{3} =SlideIDs{3}(UseTF{3});
        %ys{3} =  CompiledEmbryos3.UnivScaledProfiles.Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles.Flipped(UseTF{3},:,ch_index, master_index);
        ys{3} =  CompiledEmbryos3.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{3},:,ch_index, master_index);
        BadTF{3} =  sum(isnan(ys{3}), 2).' >min( sum(isnan(ys{3}), 2).');
        x_sample{3} = x_sample{3}(~BadTF{3});
        SlideIDs{3} = SlideIDs{3}(~BadTF{3});
        ys{3} = ys{3}(~BadTF{3},:);
        %xfits{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles.Flipped.x;
        %yfits{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.mean(:,:,ch_index, master_index);
        yfits{3} =CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        xfits{3} =CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        
        
        x_sample_control{3} = CompiledEmbryos3.DubuisEmbryoTimes;
        SlideIDs_control{3} = CompiledEmbryos3.SlideIDs;
        UseTF_control{3} = CompiledEmbryos3.ControlSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample_control{3});
        x_sample_control{3} =x_sample_control{3}(UseTF_control{3});
        SlideIDs_control{3} =SlideIDs_control{3}(UseTF_control{3});
        %ys_control{3} =   CompiledEmbryos3.UnivScaledProfiles.Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles.Flipped(UseTF_control{3},:,ch_index, master_index);
        ys_control{3} =  CompiledEmbryos3.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{3},:,ch_index, master_index);
        BadTF_control{3}=  sum(isnan(ys_control{3}), 2).' >min( sum(isnan(ys_control{3}), 2).');
        x_sample_control{3} = x_sample_control{3}(~BadTF_control{3});
        SlideIDs_control{3} = SlideIDs_control{3}(~BadTF_control{3});
        ys_control{3} = ys_control{3}(~BadTF_control{3},:);
        %xfits_control{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles.Flipped.x;
        %yfits_control{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.mean(:,:,ch_index, master_index);
        yfits_control{3} =CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        xfits_control{3} =CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        
        SetLabel{3} = AllSetInfo.SetLabels{exp_index3};
        
        for APindex =5:4:37
            ymax = max([max( yfits_control{1}(:,APindex)), max( yfits_control{2}(:,APindex)), max( yfits_control{3}(:,APindex)),...
                max( yfits{1}(:,APindex)), max( yfits{2}(:,APindex)), max( yfits{3}(:,APindex)),...
                max(ys_control{1}(:,APindex)), max(ys_control{2}(:,APindex)), max(ys_control{3}(:,APindex)),...
                max(ys{1}(:,APindex)),max(ys{2}(:,APindex)),max(ys{3}(:,APindex))])*1.1;
            ymin = min([min( yfits_control{1}(:,APindex)), min( yfits_control{2}(:,APindex)), min( yfits_control{3}(:,APindex)),...
                min( yfits{1}(:,APindex)), min( yfits{2}(:,APindex)), min( yfits{3}(:,APindex)),...
                min(ys_control{1}(:,APindex)), min(ys_control{2}(:,APindex)), min(ys_control{3}(:,APindex)),...
                min(ys{1}(:,APindex)),min(ys{2}(:,APindex)),min(ys{3}(:,APindex))]);
            DeltaFCAx = cell(2,4);
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .8]);
            set(gcf,'color','w');
            for sub_index = 1:3
                DeltaFCAx{1,sub_index} = subplot(2, 4, sub_index);
                for slide_index = unique(SlideIDs{sub_index})
                scatter(x_sample{sub_index}(SlideIDs{sub_index} == slide_index), ys{sub_index}(SlideIDs{sub_index} == slide_index,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(slide_index+1,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(slide_index+1,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                     hold on
                end
           
                plot(xfits{sub_index}, yfits{sub_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                %ymax = max([max(yfits{sub_index}(:, APindex)),max(ys{sub_index}(:,APindex))]) *1.1;
                
                grid on
                hold off
                xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
                %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
                xlim([0, 60])
                ylim([0, ymax-ymin])
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
                
                %ylim([0, min([ymax, 100])])
                
                
                DeltaFCAx{1,sub_index}.YAxis.FontSize = 14;
                DeltaFCAx{1,sub_index}.XAxis.FontSize = 14;
                
                if AllSetInfo.Replicates(exp_indices(sub_index)) ~= 0
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_indices(sub_index))),' Test Set'];
                else
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Flipped Stain Test Set'];
                end
                title(DeltaFCAx{1,sub_index}, [PlotLabel], 'FontSize', 14)
                %
                
                
                DeltaFCAx{1,sub_index}.FontSize = 14;
                
                
                
                DeltaFCAx{2,sub_index} = subplot(2, 4, sub_index+4);
                for slide_index = unique(SlideIDs{sub_index})
                scatter(x_sample_control{sub_index}(SlideIDs_control{sub_index} == slide_index), ys_control{sub_index}(SlideIDs_control{sub_index} == slide_index,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(slide_index+1,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(slide_index+1,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                     hold on
                end
                hold on
                plot(xfits_control{sub_index}, yfits_control{sub_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                %ymax = max([max(yfits_control{sub_index}(:, APindex)),max(ys_control{sub_index}(:,APindex))]) *1.1;
                
                grid on
                hold off
                xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
                
                
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
                
                ylim([0,  ymax-ymin])
                xlim([0 60])
                
                DeltaFCAx{2, sub_index}.YAxis.FontSize = 14;
                DeltaFCAx{2, sub_index}.XAxis.FontSize = 14;
                
                if AllSetInfo.Replicates(exp_indices(sub_index)) ~= 0
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_indices(sub_index))),' Control Set'];
                else
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Flipped Stain Control Set'];
                end
                title(DeltaFCAx{2, sub_index}, [PlotLabel], 'FontSize', 14)
                %
                
                
                DeltaFCAx{2, sub_index}.FontSize = 14;
            end
            % all sets plotted on single plot 
            DeltaFCAx{1,4} = subplot(2, 4, 4);
            for sub_index = 1:3
                scatter(x_sample{sub_index}, ys{sub_index}(:,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(sub_index+2,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(sub_index+2,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                hold on
            end 
            for sub_index = 1:3
                 plot(xfits{sub_index}, yfits{sub_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                 hold on 
            end
            grid on
            xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
            xlim([0, 60])
            ylim([0,  ymax-ymin])
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)   
            DeltaFCAx{1,4}.YAxis.FontSize = 14;
            DeltaFCAx{1,4}.XAxis.FontSize = 14;
            title(DeltaFCAx{1,4}, 'All Test Sets', 'FontSize', 14)
            DeltaFCAx{1,4}.FontSize = 14;

            DeltaFCAx{2,4} = subplot(2, 4, 8);
            for sub_index = 1:3
                scatter(x_sample_control{sub_index}, ys_control{sub_index}(:,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(sub_index+2,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(sub_index+2,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                hold on
            end
            for sub_index = 1:3
                plot(xfits_control{sub_index}, yfits_control{sub_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                hold on
            end
            grid on
            hold off
            xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            ylim([0,  ymax-ymin])
            xlim([0 60])
            DeltaFCAx{2, 4}.YAxis.FontSize = 14;
            DeltaFCAx{2, 4}.XAxis.FontSize = 14;
            title(DeltaFCAx{2, 4}, 'All Control Sets', 'FontSize', 14)
            DeltaFCAx{2, 4}.FontSize = 14;
      
            
            sgtitle(['AP: ', num2str(round(APbins(APindex), 2))], 'FontSize', 16)
            if ~isdir([AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'SlideRescaledDorsalAvgAPProfiles', filesep, 'MasterSet', num2str(master_index)])
                mkdir([AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'SlideRescaledDorsalAvgAPProfiles', filesep, 'MasterSet', num2str(master_index)])
            end
            outpath2 = [AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'SlideRescaledDorsalAvgAPProfiles',...
                filesep, 'MasterSet', num2str(master_index), filesep,...
                ChannelNames{ch_index}, '_T', strrep(num2str(temp), '.', '_'), 'C_AP',num2str(APindex), '.png' ];
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
    end
end

%%
i_list = [1 1 2];
j_list = [2 3 3];

i_list = [1 1 2 4 4 5 7 7 8 10 10 11 13 13 14];
j_list = [2 3 3 5 6 6 8 9 9 11 12 12 14 15 15];

i_list = [1 4 7 10 13];
j_list = [2 5 8 11 14];
%
% i_list = [4 4 5 7 7 8 10 10 11 13 13 14];
% j_list = [5 6 6 8 9 9 11 12 12 14 15 15];
% %
% 
% i_list = [1 2 4 5 7 8 10 11 13 14];
% j_list = [3 3 6 6 9 9 12 12 15 15];

% i_list = [1 1 2 4 4 5 7 7 8 10 10 11 13 13 14];
% j_list = [2 3 3 5 6 6 8 9 9 11 12 12 14 15 15];


for ch_index = [3 5]
for list_index = 1:length(i_list)
    for prof_index = 1:8
        exp_index = i_list(list_index);
        CompiledEmbryos = AllCompiledEmbryos{exp_index};
        x_sample = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample);
        x_sample =x_sample(UseTF);
        ys = CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF, :, ch_index, prof_index);
        BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        xfits = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:, :, ch_index, prof_index);
        
        
        exp_index2 = j_list(list_index);
        CompiledEmbryos2 = AllCompiledEmbryos{exp_index2};
        x_sample_control = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF_control = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control);
        x_sample_control =x_sample_control(UseTF_control);
        ys_control = CompiledEmbryos2.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control, :, ch_index, prof_index);BadTF_control=  sum(isnan(ys_control), 2).' >min( sum(isnan(ys_control), 2).');
        x_sample_control = x_sample_control(~BadTF_control);
        ys_control = ys_control(~BadTF_control,:);
        xfits_control = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control =  CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:, :, ch_index, prof_index);
        
        APmat = repmat(APbins, length(xfits), 1);
        xmat = repmat(xfits.', 1, NumAPbins);
        xindex_mat = repmat((1:length(xfits)).', 1, NumAPbins);
        
        yfits = yfits(1:size(yfits, 1),:).';
        yfits = yfits(:).';
        yfits_control = yfits_control(1:size(yfits_control, 1),:).';
        yfits_control = yfits_control(:).';
        xflat = xmat(1:size(xmat, 1),:);
        xflat = xflat(:).';
        xindex_flat = xindex_mat(1:size(xindex_mat, 1),:);
        xindex_flat = xindex_flat(:).';
        
        TFkeep = ~isnan(yfits) & ~isnan(yfits_control);
        x_plot = xindex_flat(TFkeep);
        yfits = yfits(TFkeep);
        yfits2 = yfits_control(TFkeep);
        
        SS_res = sum((yfits2-yfits).^2);
        SS_tot = sum((yfits2-mean(yfits2)).^2);
        Rsquared = 1-SS_res/SS_tot;
        
        SetLabel = AllSetInfo.SetLabels{exp_index};
        SetLabel2 = AllSetInfo.SetLabels{exp_index2};
        DubuisTimeRange = xfits_control;
        colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
        FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
        BinsToFit = zeros(NChannels, NumAPbins, 'logical');
        BinsToFit(2,9:33) = true;
        BinsToFit(3,5:13) = true;
        BinsToFit(4,22:37) = true;
        BinsToFit(5,5:21) = true;
        
        
        close all
        DeltaFCFixCorrectedFig =figure(1);
        set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.05, 0.05, .6, .6]);
        set(gcf,'color','w');
        DeltaFCAx = axes(DeltaFCFixCorrectedFig);
        map = colormap(colorsDubuisTimes);
        h = colorbar;
        % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
        
        colorTitleHandle = get(h,'Title');
        titleString = '\del (\mum)';
        titleString = 'time (min)';
        set(colorTitleHandle ,'String',titleString);
        h.Ticks =  FractonalDubuisTimeRange(1:5:71)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
        h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0','65.0', '70.0'};
        hold on
        for plot_index = 1:5:length(x_plot)
            scatter(yfits(plot_index),yfits2(plot_index),...
                75, 'MarkerFaceColor', colorsDubuisTimes(x_plot(plot_index),:),...
                'MarkerEdgeColor', colorsDubuisTimes(x_plot(plot_index),:));
            hold on
        end
        
        ymax = max([max(yfits),max(yfits2)]) *1.1;
        plot([0 ymax], [0 ymax], 'k', 'LineWidth', 2.0);
        
        
        grid on
        hold off
        if AllSetInfo.Replicates(exp_index) ~= 0
            xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Test Set'];
        else
            xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Test Set'];
        end
        xlabel(xLabel, 'FontSize', 16)
        
        if AllSetInfo.Replicates(exp_index2) ~= 0
            yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index2)),' Test Set'];
        else
            yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Flipped Stain Test Set'];
        end
        ylabel(yLabel, 'FontSize', 16)
        %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
        xlim([0, min([ymax, 150])])
        
        title(DeltaFCAx, {['Master Profile ', num2str(prof_index),', R^2: ',num2str(round(Rsquared, 3))],['Control Scaling  ', ChannelNames{ch_index}, ' (AU)']}, 'FontSize', 18);
        
        ylim([0, min([ymax, 150])])
        
        DeltaFCAx.YAxis.FontSize = 16;
        DeltaFCAx.XAxis.FontSize = 16;
        
        %
        
        
        DeltaFCAx.FontSize = 16;
        
        %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
          outpath = [AllSetsProfFigPath, filesep,'TwoSlideRescaledSetComps' , filesep, SetLabel, '_Test_',SetLabel2,'_Test_ControlScaling_',ChannelNames{ch_index},'_SlideRescaledMaster',num2str(prof_index), '_DubuisColored.png'];
        saveas(DeltaFCFixCorrectedFig,outpath);
    end
end
end
% i_list = [1 1 2];
% j_list = [2 3 3];
%

%%
ch_index = 3;
for list_index = 1:length(i_list)
    for prof_index = 1:8
        exp_index = i_list(list_index);
        CompiledEmbryos = AllCompiledEmbryos{exp_index};
        x_sample = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample);
        x_sample =x_sample(UseTF);
        ys = CompiledEmbryos.UnivScaledProfiles.FitZeroedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF, :, ch_index, prof_index);
        BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
        x_sample = x_sample(~BadTF);
        ys = ys(~BadTF,:);
        xfits = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitZeroedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitZeroedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:, :, ch_index, prof_index);
        
        
        exp_index2 = j_list(list_index);
        CompiledEmbryos2 = AllCompiledEmbryos{exp_index2};
        x_sample_control = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF_control = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control);
        x_sample_control =x_sample_control(UseTF_control);
        ys_control = CompiledEmbryos2.UnivScaledProfiles.FitZeroedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control, :, ch_index, prof_index);BadTF_control=  sum(isnan(ys_control), 2).' >min( sum(isnan(ys_control), 2).');
        x_sample_control = x_sample_control(~BadTF_control);
        ys_control = ys_control(~BadTF_control,:);
        xfits_control = CompiledEmbryos2.BootstrappedUnivScaledProfiles.FitZeroedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control =  CompiledEmbryos2.BootstrappedUnivScaledProfiles.FitZeroedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:, :, ch_index, prof_index);
        
        APmat = repmat(APbins, length(xfits), 1);
        xmat = repmat(xfits.', 1, NumAPbins);
        xindex_mat = repmat((1:length(xfits)).', 1, NumAPbins);
        
        yfits = yfits(1:size(yfits, 1),:).';
        yfits = yfits(:).';
        yfits_control = yfits_control(1:size(yfits_control, 1),:).';
        yfits_control = yfits_control(:).';
        xflat = xmat(1:size(xmat, 1),:);
        xflat = xflat(:).';
        xindex_flat = xindex_mat(1:size(xindex_mat, 1),:);
        xindex_flat = xindex_flat(:).';
        
        TFkeep = ~isnan(yfits) & ~isnan(yfits_control);
        x_plot = xindex_flat(TFkeep);
        yfits = yfits(TFkeep);
        yfits2 = yfits_control(TFkeep);
        
        SS_res = sum((yfits2-yfits).^2);
        SS_tot = sum((yfits2-mean(yfits2)).^2);
        Rsquared = 1-SS_res/SS_tot;
        
        SetLabel = AllSetInfo.SetLabels{exp_index};
        SetLabel2 = AllSetInfo.SetLabels{exp_index2};
        DubuisTimeRange = xfits_control;
        colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
        FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
        BinsToFit = zeros(NChannels, NumAPbins, 'logical');
        BinsToFit(2,9:33) = true;
        BinsToFit(3,5:13) = true;
        BinsToFit(4,22:37) = true;
        BinsToFit(5,5:21) = true;
        
        
        close all
        DeltaFCFixCorrectedFig =figure(1);
        set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.05, 0.05, .6, .6]);
        set(gcf,'color','w');
        DeltaFCAx = axes(DeltaFCFixCorrectedFig);
        map = colormap(colorsDubuisTimes);
        h = colorbar;
        % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
        
        colorTitleHandle = get(h,'Title');
        titleString = '\del (\mum)';
        titleString = 'time (min)';
        set(colorTitleHandle ,'String',titleString);
        h.Ticks =  FractonalDubuisTimeRange(1:5:71)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
        h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0','65.0', '70.0'};
        hold on
        for plot_index = 1:5:length(x_plot)
            scatter(yfits(plot_index),yfits2(plot_index),...
                75, 'MarkerFaceColor', colorsDubuisTimes(x_plot(plot_index),:),...
                'MarkerEdgeColor', colorsDubuisTimes(x_plot(plot_index),:));
            hold on
        end
        
        ymax = max([max(yfits),max(yfits2)]) *1.1;
        plot([0 ymax], [0 ymax], 'k', 'LineWidth', 2.0);
        
        
        grid on
        hold off
        if AllSetInfo.Replicates(exp_index) ~= 0
            xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Test Set'];
        else
            xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Test Set'];
        end
        xlabel(xLabel, 'FontSize', 16)
        
        if AllSetInfo.Replicates(exp_index2) ~= 0
            yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index2)),' Test Set'];
        else
            yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Flipped Stain Test Set'];
        end
        ylabel(yLabel, 'FontSize', 16)
        %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
        xlim([0, min([ymax, 150])])
        
        title(DeltaFCAx, {['Master Profile ', num2str(prof_index),', R^2: ',num2str(round(Rsquared, 3))],['Control Scaling  ', ChannelNames{ch_index}, ' (AU)']}, 'FontSize', 18);
        
        ylim([0, min([ymax, 150])])
        
        DeltaFCAx.YAxis.FontSize = 16;
        DeltaFCAx.XAxis.FontSize = 16;
        
        %
        
        
        DeltaFCAx.FontSize = 16;
        
        %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
        outpath = [AllSetsProfFigPath, filesep,'TwoFitZeroedSlideRescaledSetComps' , filesep, SetLabel, '_Test_',SetLabel2,'_Test_ControlScaling_',ChannelNames{ch_index},'_FitZeroedSlideRescaledMaster',num2str(prof_index), '_DubuisColored.png'];
        saveas(DeltaFCFixCorrectedFig,outpath);
    end
end

%%


i_list = [1 1 2];
j_list = [2 3 3];

i_list = [1 1 2 4 4 5 7 7 8 10 10 11 13 13 14];
j_list = [2 3 3 5 6 6 8 9 9 11 12 12 14 15 15];
%
% i_list = [4 4 5 7 7 8 10 10 11 13 13 14];
% j_list = [5 6 6 8 9 9 11 12 12 14 15 15];
%
%
% i_list = [1 1 2];
% j_list = [2 3 3];
%



ch_index = 3;
for list_index = 1:length(i_list)
    
    exp_index = i_list(list_index);
    CompiledEmbryos = AllCompiledEmbryos{exp_index};
    x_sample = CompiledEmbryos.DubuisEmbryoTimes;
    UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample);
    x_sample =x_sample(UseTF);
    ys = CompiledEmbryos.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF, :, ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    xfits = CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.x;
    yfits = CompiledEmbryos.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:, :, ch_index);
    
    
    exp_index2 = j_list(list_index);
    CompiledEmbryos2 = AllCompiledEmbryos{exp_index2};
    x_sample_control = CompiledEmbryos2.DubuisEmbryoTimes;
    UseTF_control = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control);
    x_sample_control =x_sample_control(UseTF_control);
    ys_control = CompiledEmbryos2.UnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control, :, ch_index);
    BadTF_control=  sum(isnan(ys_control), 2).' >min( sum(isnan(ys_control), 2).');
    x_sample_control = x_sample_control(~BadTF_control);
    ys_control = ys_control(~BadTF_control,:);
    xfits_control = CompiledEmbryos2.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep1.x;
    yfits_control =  CompiledEmbryos2.BootstrappedUnivScaledProfiles.TestSetFitSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:, :, ch_index);
    
    APmat = repmat(APbins, length(xfits), 1);
    xmat = repmat(xfits.', 1, NumAPbins);
    xindex_mat = repmat((1:length(xfits)).', 1, NumAPbins);
    
    yfits = yfits(1:size(yfits, 1),:).';
    yfits = yfits(:).';
    yfits_control = yfits_control(1:size(yfits_control, 1),:).';
    yfits_control = yfits_control(:).';
    xflat = xmat(1:size(xmat, 1),:);
    xflat = xflat(:).';
    xindex_flat = xindex_mat(1:size(xindex_mat, 1),:);
    xindex_flat = xindex_flat(:).';
    
    TFkeep = ~isnan(yfits) & ~isnan(yfits_control);
    x_plot = xindex_flat(TFkeep);
    yfits = yfits(TFkeep);
    yfits2 = yfits_control(TFkeep);
    
    SS_res = sum((yfits2-yfits).^2);
    SS_tot = sum((yfits2-mean(yfits2)).^2);
    Rsquared = 1-SS_res/SS_tot;
    
    SetLabel = AllSetInfo.SetLabels{exp_index};
    SetLabel2 = AllSetInfo.SetLabels{exp_index2};
    DubuisTimeRange = xfits_control;
    colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
    FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
    BinsToFit = zeros(NChannels, NumAPbins, 'logical');
    BinsToFit(2,9:33) = true;
    BinsToFit(3,5:13) = true;
    BinsToFit(4,22:37) = true;
    BinsToFit(5,5:21) = true;
    
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.05, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    map = colormap(colorsDubuisTimes);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    
    colorTitleHandle = get(h,'Title');
    titleString = '\del (\mum)';
    titleString = 'time (min)';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks =  FractonalDubuisTimeRange(1:50:701)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
    h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0','65.0', '70.0'};
    hold on
    for plot_index = 1:20:length(x_plot)
        scatter(yfits(plot_index),yfits2(plot_index),...
            75, 'MarkerFaceColor', colorsDubuisTimes(x_plot(plot_index),:),...
            'MarkerEdgeColor', colorsDubuisTimes(x_plot(plot_index),:));
        hold on
    end
    
    ymax = max([max(yfits),max(yfits2)]) *1.1;
    plot([0 ymax], [0 ymax], 'k', 'LineWidth', 2.0);
    
    
    grid on
    hold off
    if AllSetInfo.Replicates(exp_index) ~= 0
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Test Set'];
    else
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Test Set'];
    end
    xlabel(xLabel, 'FontSize', 16)
    
    if AllSetInfo.Replicates(exp_index2) ~= 0
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index2)),' Test Set'];
    else
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Flipped Stain Test Set'];
    end
    ylabel(yLabel, 'FontSize', 16)
    %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
    xlim([0, min([ymax, 150])])
    
    title(DeltaFCAx, {['Rep2 Master Profile, R^2: ',num2str(round(Rsquared, 3))],['Test Scaling  ', ChannelNames{ch_index}, ' (AU)']}, 'FontSize', 18);
    
    ylim([0, min([ymax, 150])])
    
    DeltaFCAx.YAxis.FontSize = 16;
    DeltaFCAx.XAxis.FontSize = 16;
    
    %
    
    
    DeltaFCAx.FontSize = 16;
    if ~isdir([AllSetsProfFigPath, filesep,'Rep2Comps'])
        mkdir([AllSetsProfFigPath, filesep,'Rep2Comps' ])
    end
    %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
    outpath = [AllSetsProfFigPath, filesep,'Rep2Comps' , filesep, SetLabel, '_Test_',SetLabel2,'_Test_SlideRescaled_Rep2Scaling_',ChannelNames{ch_index}, '_DubuisColored.png'];
    saveas(DeltaFCFixCorrectedFig,outpath);
end



%% Inter-temperature plots
temps = [17.5 20 22.5 25 27.5];
for master_index = 1:8
for ch_index = [3 5]
    exp_indices = cell(1, length(temps));
    exp_index = zeros(1, length(temps));
    exp_index2 = zeros(1, length(temps));
    exp_index3 = zeros(1, length(temps));
    x_sample = cell(3, length(temps));
    UseTF = cell(3, length(temps));
    ys = cell(3, length(temps));
    BadTF = cell(3, length(temps));
    xfits = cell(3, length(temps));
    yfits = cell(3, length(temps));
    x_sample_control = cell(3, length(temps));
    UseTF_control = cell(3, length(temps));
    ys_control = cell(3, length(temps));
    BadTF_control = cell(3, length(temps));
    xfits_control = cell(3, length(temps));
    yfits_control = cell(3, length(temps));
    SetLabel = cell(3, length(temps));
    for temp_index = 1:length(temps)
        exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
        exp_index(temp_index) =  exp_indices{temp_index}(1);
        exp_index2(temp_index) =  exp_indices{temp_index}(2);
        exp_index3(temp_index) =  exp_indices{temp_index}(3);
 
        
        CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
        x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
        x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
        ys{1, temp_index} =  CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{1, temp_index},:,ch_index, master_index);
        BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).' >min( sum(isnan(ys{1, temp_index}), 2).');
        x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
        ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
        xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        
        
        x_sample_control{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF_control{1, temp_index} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1, temp_index});
        x_sample_control{1, temp_index} =x_sample_control{1, temp_index}(UseTF_control{1, temp_index});
        ys_control{1, temp_index} =  CompiledEmbryos.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{1, temp_index},:,ch_index, master_index);
        BadTF_control{1, temp_index} =  sum(isnan(ys_control{1, temp_index}), 2).' >min( sum(isnan(ys_control{1, temp_index}), 2).');
        x_sample_control{1, temp_index} = x_sample_control{1, temp_index}(~BadTF_control{1, temp_index});
        ys_control{1, temp_index} = ys_control{1, temp_index}(~BadTF_control{1, temp_index},:);
        xfits_control{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
        
        CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
        x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
        x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
        %ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF{2, temp_index},:,ch_index, master_index);
        ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{2, temp_index},:,ch_index, master_index);
        BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' >min( sum(isnan(ys{2, temp_index}), 2).');
        x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
        ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
        %xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
        %yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:,:,ch_index, master_index);
        xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        
        
        x_sample_control{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF_control{2, temp_index} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2, temp_index});
        x_sample_control{2, temp_index} =x_sample_control{2, temp_index}(UseTF_control{2, temp_index});
        %ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control{2, temp_index},:,ch_index, master_index);
        ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{2, temp_index},:,ch_index, master_index);
        BadTF_control{2, temp_index}=  sum(isnan(ys_control{2, temp_index}), 2).' >min( sum(isnan(ys_control{2, temp_index}), 2).');
        x_sample_control{2, temp_index} = x_sample_control{2, temp_index}(~BadTF_control{2, temp_index});
        ys_control{2, temp_index} = ys_control{2, temp_index}(~BadTF_control{2, temp_index},:);
        %xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
        %yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean(:,:,ch_index, master_index);
        xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
        
        CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
        x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
        x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
        %ys{3, temp_index} =  CompiledEmbryos3.UnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped(UseTF{3, temp_index},:,ch_index, master_index);
        ys{3, temp_index} =  CompiledEmbryos3.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{3, temp_index},:,ch_index, master_index);
        BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).' >min( sum(isnan(ys{3, temp_index}), 2).');
        x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
        ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
        %xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.x;
        %yfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.mean(:,:,ch_index, master_index);
        xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        
        
        x_sample_control{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF_control{3, temp_index} = CompiledEmbryos3.ControlSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample_control{3, temp_index});
        x_sample_control{3, temp_index} =x_sample_control{3, temp_index}(UseTF_control{3, temp_index});
        %ys_control{3, temp_index} =   CompiledEmbryos3.UnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped(UseTF_control{3, temp_index},:,ch_index, master_index);
        ys_control{3, temp_index} =  CompiledEmbryos3.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{3, temp_index},:,ch_index, master_index);
        BadTF_control{3, temp_index}=  sum(isnan(ys_control{3, temp_index}), 2).' >min( sum(isnan(ys_control{3, temp_index}), 2).');
        x_sample_control{3, temp_index} = x_sample_control{3, temp_index}(~BadTF_control{3, temp_index});
        ys_control{3, temp_index} = ys_control{3, temp_index}(~BadTF_control{3, temp_index},:);
        %xfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.x;
        %yfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.mean(:,:,ch_index, master_index);
        xfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
    end
    for APindex =5:4:37
        
        DeltaFCAx = cell(2,5);
        close all
        DeltaFCFixCorrectedFig =figure(1);
        set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .8]);
        set(gcf,'color','w');
        
        for temp_index = 1:length(temps)
           ymaxes(temp_index) = max([max( yfits_control{1, temp_index}(:)), max( yfits_control{2, temp_index}(:)), max( yfits_control{3, temp_index}(:)),...
                max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:)),...
                max(ys_control{1, temp_index}(:)), max(ys_control{2, temp_index}(:)), max(ys_control{3, temp_index}(:)),...
                max(ys{1, temp_index}(:)),max(ys{2, temp_index}(:)),max(ys{3, temp_index}(:))])*1.1;
            ymins(temp_index) = min([min( yfits_control{1, temp_index}(:)), min( yfits_control{2, temp_index}(:)), min( yfits_control{3, temp_index}(:)),...
                min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{3, temp_index}(:)),...
                min(ys_control{1, temp_index}(:)), min(ys_control{2, temp_index}(:)), min(ys_control{3, temp_index}(:)),...
                min(ys{1, temp_index}(:)),min(ys{2, temp_index}(:)),min(ys{3, temp_index}(:))]);
        end
        
        for temp_index = 1:length(temps)
            ymax = max(ymaxes);
            ymin = min(ymins);
            % all sets plotted on single plot
            DeltaFCAx{1,5} = subplot(2, 5, temp_index);
            for sub_index = 1:3
                scatter(x_sample{sub_index, temp_index}, ys{sub_index, temp_index}(:,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(sub_index+2,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(sub_index+2,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                hold on
            end
            for sub_index = 1:3
                plot(xfits{sub_index, temp_index}, yfits{sub_index, temp_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                hold on
            end
            grid on
            xlabel({'Time (Dubuis) into', 'cycle 14 (min)'}, 'FontSize', 14)
            xlim([0, 60])
            ylim([0,  ymax-ymin])
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            DeltaFCAx{1,5}.YAxis.FontSize = 14;
            DeltaFCAx{1,5}.XAxis.FontSize = 14;
            title(DeltaFCAx{1,5}, ['T = ',num2str(temps(temp_index)), '°C Test Sets'], 'FontSize', 14)
            DeltaFCAx{1,5}.FontSize = 14;
            
            DeltaFCAx{2,5} = subplot(2, 5, temp_index+5);
            for sub_index = 1:3
                scatter(x_sample_control{sub_index, temp_index}, ys_control{sub_index, temp_index}(:,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(sub_index+2,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(sub_index+2,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                hold on
            end
            for sub_index = 1:3
                plot(xfits_control{sub_index, temp_index}, yfits_control{sub_index, temp_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                hold on
            end
            grid on
            hold off
            xlabel({'Time (Dubuis) into', 'cycle 14 (min)'}, 'FontSize', 14)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            ylim([0,  ymax-ymin])
            xlim([0 60])
            DeltaFCAx{2, 5}.YAxis.FontSize = 14;
            DeltaFCAx{2, 5}.XAxis.FontSize = 14;
            title(DeltaFCAx{2,5}, ['T = ',num2str(temps(temp_index)), '°C Control Sets'], 'FontSize', 14)
            DeltaFCAx{2, 5}.FontSize = 14;
            
        end
        sgtitle(['AP: ', num2str(round(APbins(APindex), 2)), ', Master Set: ', num2str(master_index)], 'FontSize', 16)
        if ~isdir([AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'SlideRescaledDorsalAvgAPProfiles', filesep, 'MasterSet', num2str(master_index)])
            mkdir([AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'SlideRescaledDorsalAvgAPProfiles', filesep, 'MasterSet', num2str(master_index)])
        end
        outpath2 = [AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'SlideRescaledDorsalAvgAPProfiles',...
            filesep, 'MasterSet', num2str(master_index), filesep,...
            ChannelNames{ch_index}, '_AllTemps_AP',num2str(APindex), '.png' ];
        saveas(DeltaFCFixCorrectedFig,outpath2);
        
    end
end
end



for master_index = 1:12
for ch_index = 3
    exp_indices = cell(1, length(temps));
    exp_index = zeros(1, length(temps));
    exp_index2 = zeros(1, length(temps));
    exp_index3 = zeros(1, length(temps));
    x_sample = cell(3, length(temps));
    UseTF = cell(3, length(temps));
    ys = cell(3, length(temps));
    BadTF = cell(3, length(temps));
    xfits = cell(3, length(temps));
    yfits = cell(3, length(temps));
    x_sample_control = cell(3, length(temps));
    UseTF_control = cell(3, length(temps));
    ys_control = cell(3, length(temps));
    BadTF_control = cell(3, length(temps));
    xfits_control = cell(3, length(temps));
    yfits_control = cell(3, length(temps));
    SetLabel = cell(3, length(temps));
    for temp_index = 1:length(temps)
        exp_indices{temp_index} = find(AllSetInfo.Temperatures == temps(temp_index));
        exp_index(temp_index) =  exp_indices{temp_index}(1);
        exp_index2(temp_index) =  exp_indices{temp_index}(2);
        exp_index3(temp_index) =  exp_indices{temp_index}(3);
 
        
        CompiledEmbryos = AllCompiledEmbryos{exp_index(temp_index)};
        x_sample{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF{1, temp_index} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1, temp_index});
        x_sample{1, temp_index} =x_sample{1, temp_index}(UseTF{1, temp_index});
        ys{1, temp_index} =  CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{1, temp_index},:,ch_index, master_index);
        BadTF{1, temp_index} =  sum(isnan(ys{1, temp_index}), 2).' >min( sum(isnan(ys{1, temp_index}), 2).');
        x_sample{1, temp_index} = x_sample{1, temp_index}(~BadTF{1, temp_index});
        ys{1, temp_index} = ys{1, temp_index}(~BadTF{1, temp_index},:);
        xfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        
        
        x_sample_control{1, temp_index} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF_control{1, temp_index} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1, temp_index});
        x_sample_control{1, temp_index} =x_sample_control{1, temp_index}(UseTF_control{1, temp_index});
        ys_control{1, temp_index} =  CompiledEmbryos.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{1, temp_index},:,ch_index, master_index);
        BadTF_control{1, temp_index} =  sum(isnan(ys_control{1, temp_index}), 2).' >min( sum(isnan(ys_control{1, temp_index}), 2).');
        x_sample_control{1, temp_index} = x_sample_control{1, temp_index}(~BadTF_control{1, temp_index});
        ys_control{1, temp_index} = ys_control{1, temp_index}(~BadTF_control{1, temp_index},:);
        xfits_control{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control{1, temp_index} = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        SetLabel{1, temp_index} = AllSetInfo.SetLabels{exp_index(temp_index)};
        
        CompiledEmbryos2 = AllCompiledEmbryos{exp_index2(temp_index)};
        x_sample{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF{2, temp_index} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2, temp_index});
        x_sample{2, temp_index} =x_sample{2, temp_index}(UseTF{2, temp_index});
        %ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF{2, temp_index},:,ch_index, master_index);
        ys{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{2, temp_index},:,ch_index, master_index);
        BadTF{2, temp_index} =  sum(isnan(ys{2, temp_index}), 2).' >min( sum(isnan(ys{2, temp_index}), 2).');
        x_sample{2, temp_index} = x_sample{2, temp_index}(~BadTF{2, temp_index});
        ys{2, temp_index} = ys{2, temp_index}(~BadTF{2, temp_index},:);
        %xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
        %yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.TestSet.mean(:,:,ch_index, master_index);
        xfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        
        
        x_sample_control{2, temp_index} = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF_control{2, temp_index} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2, temp_index});
        x_sample_control{2, temp_index} =x_sample_control{2, temp_index}(UseTF_control{2, temp_index});
        %ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2(UseTF_control{2, temp_index},:,ch_index, master_index);
        ys_control{2, temp_index} =  CompiledEmbryos2.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{2, temp_index},:,ch_index, master_index);
        BadTF_control{2, temp_index}=  sum(isnan(ys_control{2, temp_index}), 2).' >min( sum(isnan(ys_control{2, temp_index}), 2).');
        x_sample_control{2, temp_index} = x_sample_control{2, temp_index}(~BadTF_control{2, temp_index});
        ys_control{2, temp_index} = ys_control{2, temp_index}(~BadTF_control{2, temp_index},:);
        %xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.x;
        %yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Rep2.ControlSet.mean(:,:,ch_index, master_index);
        xfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control{2, temp_index} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        SetLabel{2, temp_index} = AllSetInfo.SetLabels{exp_index2(temp_index)};
        
        CompiledEmbryos3 = AllCompiledEmbryos{exp_index3(temp_index)};
        x_sample{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF{3, temp_index} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3, temp_index});
        x_sample{3, temp_index} =x_sample{3, temp_index}(UseTF{3, temp_index});
        %ys{3, temp_index} =  CompiledEmbryos3.UnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped(UseTF{3, temp_index},:,ch_index, master_index);
        ys{3, temp_index} =  CompiledEmbryos3.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{3, temp_index},:,ch_index, master_index);
        BadTF{3, temp_index} =  sum(isnan(ys{3, temp_index}), 2).' >min( sum(isnan(ys{3, temp_index}), 2).');
        x_sample{3, temp_index} = x_sample{3, temp_index}(~BadTF{3, temp_index});
        ys{3, temp_index} = ys{3, temp_index}(~BadTF{3, temp_index},:);
        %xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.x;
        %yfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.TestSet.mean(:,:,ch_index, master_index);
        xfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index, master_index);
        
        
        x_sample_control{3, temp_index} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF_control{3, temp_index} = CompiledEmbryos3.ControlSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample_control{3, temp_index});
        x_sample_control{3, temp_index} =x_sample_control{3, temp_index}(UseTF_control{3, temp_index});
        %ys_control{3, temp_index} =   CompiledEmbryos3.UnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped(UseTF_control{3, temp_index},:,ch_index, master_index);
        ys_control{3, temp_index} =  CompiledEmbryos3.UnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{3, temp_index},:,ch_index, master_index);
        BadTF_control{3, temp_index}=  sum(isnan(ys_control{3, temp_index}), 2).' >min( sum(isnan(ys_control{3, temp_index}), 2).');
        x_sample_control{3, temp_index} = x_sample_control{3, temp_index}(~BadTF_control{3, temp_index});
        ys_control{3, temp_index} = ys_control{3, temp_index}(~BadTF_control{3, temp_index},:);
        %xfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.x;
        %yfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles.Flipped.ControlSet.mean(:,:,ch_index, master_index);
        xfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
        yfits_control{3, temp_index} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index, master_index);
        SetLabel{3, temp_index} = AllSetInfo.SetLabels{exp_index3(temp_index)};
    end
    for APindex =5:4:37
        
        DeltaFCAx = cell(2,5);
        close all
        DeltaFCFixCorrectedFig =figure(1);
        set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .9, .8]);
        set(gcf,'color','w');
        
        for temp_index = 1:length(temps)
           ymaxes(temp_index) = max([max( yfits_control{1, temp_index}(:)), max( yfits_control{2, temp_index}(:)), max( yfits_control{3, temp_index}(:)),...
                max( yfits{1, temp_index}(:)), max( yfits{2, temp_index}(:)), max( yfits{3, temp_index}(:)),...
                max(ys_control{1, temp_index}(:)), max(ys_control{2, temp_index}(:)), max(ys_control{3, temp_index}(:)),...
                max(ys{1, temp_index}(:)),max(ys{2, temp_index}(:)),max(ys{3, temp_index}(:))])*1.1;
            ymins(temp_index) = min([min( yfits_control{1, temp_index}(:)), min( yfits_control{2, temp_index}(:)), min( yfits_control{3, temp_index}(:)),...
                min( yfits{1, temp_index}(:)), min( yfits{2, temp_index}(:)), min( yfits{3, temp_index}(:)),...
                min(ys_control{1, temp_index}(:)), min(ys_control{2, temp_index}(:)), min(ys_control{3, temp_index}(:)),...
                min(ys{1, temp_index}(:)),min(ys{2, temp_index}(:)),min(ys{3, temp_index}(:))]);
        end
        
        for temp_index = 1:length(temps)
            ymax = max(ymaxes);
            ymin = min(ymins);
            % all sets plotted on single plot
            DeltaFCAx{1,5} = subplot(2, 5, temp_index);
            for sub_index = 1:3
                scatter(x_sample{sub_index, temp_index}, ys{sub_index, temp_index}(:,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(sub_index+2,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(sub_index+2,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                hold on
            end
            for sub_index = 1:3
                plot(xfits{sub_index, temp_index}, yfits{sub_index, temp_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                hold on
            end
            grid on
            xlabel({'Time (Dubuis) into', 'cycle 14 (min)'}, 'FontSize', 14)
            xlim([0, 60])
            ylim([0,  ymax-ymin])
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            DeltaFCAx{1,5}.YAxis.FontSize = 14;
            DeltaFCAx{1,5}.XAxis.FontSize = 14;
            title(DeltaFCAx{1,5}, ['T = ',num2str(temps(temp_index)), '°C Test Sets'], 'FontSize', 14)
            DeltaFCAx{1,5}.FontSize = 14;
            
            DeltaFCAx{2,5} = subplot(2, 5, temp_index+5);
            for sub_index = 1:3
                scatter(x_sample_control{sub_index, temp_index}, ys_control{sub_index, temp_index}(:,APindex)-ymin,...
                    100, 'MarkerFaceColor',gap_colors(sub_index+2,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(sub_index+2,:),...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);% [.7 .7 .7]);
                hold on
            end
            for sub_index = 1:3
                plot(xfits_control{sub_index, temp_index}, yfits_control{sub_index, temp_index}(:,APindex)-ymin,['k', SetLineStyles{sub_index}], 'LineWidth', 2.0);
                hold on
            end
            grid on
            hold off
            xlabel({'Time (Dubuis) into', 'cycle 14 (min)'}, 'FontSize', 14)
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
            ylim([0,  ymax-ymin])
            xlim([0 60])
            DeltaFCAx{2, 5}.YAxis.FontSize = 14;
            DeltaFCAx{2, 5}.XAxis.FontSize = 14;
            title(DeltaFCAx{2,5}, ['T = ',num2str(temps(temp_index)), '°C Control Sets'], 'FontSize', 14)
            DeltaFCAx{2, 5}.FontSize = 14;
            
        end
        sgtitle(['AP: ', num2str(round(APbins(APindex), 2)), ', Master Set: ', num2str(master_index)], 'FontSize', 16)
        if ~isdir([AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'FitSlideRescaledDorsalAvgAPProfiles', filesep, 'MasterSet', num2str(master_index)])
            mkdir([AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'FitSlideRescaledDorsalAvgAPProfiles', filesep, 'MasterSet', num2str(master_index)])
        end
        outpath2 = [AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'FitSlideRescaledDorsalAvgAPProfiles',...
            filesep, 'MasterSet', num2str(master_index), filesep,...
            ChannelNames{ch_index}, '_AllTemps_AP',num2str(APindex), '.png' ];
        saveas(DeltaFCFixCorrectedFig,outpath2);
        
    end
end
end

%%
ProfMaster2s = NaN(71, 41, 15);
ProfMaster2s(:,:,1) = AllCompiledEmbryos{1}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean(:,:,3,1);
ProfMaster2s(:,:,2) = AllCompiledEmbryos{2}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean(:,:,3,1);
ProfMaster2s(:,:,3) = AllCompiledEmbryos{3}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean(:,:,3,1);
ProfMaster2s(:,:,10) = AllCompiledEmbryos{10}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,3,1);
ProfMaster2s(:,:,11) = AllCompiledEmbryos{11}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,3,1);
ProfMaster2s(:,:,12) = AllCompiledEmbryos{12}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,3,1);
ProfMaster2s(:,:,14) = AllCompiledEmbryos{14}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,3,1);
ProfMaster2s(:,:,15) = AllCompiledEmbryos{15}.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,3,1);
temp_indices = [4 4 4 5 5 5 3 3 3 2 2 2 1 1 1];
DubuisTimeRange = 0:70;
colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
TempRange = 1:5;
colorsTemperature = hsv(length(TempRange)); % Colormap "jet" is another option
FractonalTempRange = TempRange/max(TempRange);
close all
for time_index = 20:5:60
    
    
    DeltaFCFixCorrectedFig =figure(time_index);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.05, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    map = colormap(colorsTemperature);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    
    colorTitleHandle = get(h,'Title');
    titleString = 'T (ºC)';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks =  FractonalTempRange-0.5*(FractonalTempRange(2)-FractonalTempRange(1)); %Create 8 ticks from zero to 1
    h.TickLabels = {'17.5', '20.0','22.5','25.0','27.5'};%,'25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0','65.0', '70.0'};
    hold on
    for i = 1:15
        plot(APbins, ProfMaster2s(time_index+1, :, i), 'Color', colorsTemperature(temp_indices(i),:),...
            'LineWidth', 2.0)
        
        
        
    end
    
    ymax = max(ProfMaster2s(:)) *1.1;
    %plot([0 ymax], [0 ymax], 'k', 'LineWidth', 2.0);
    
    
    grid on
    hold off
    if AllSetInfo.Replicates(exp_index) ~= 0
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Test Set'];
    else
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Test Set'];
    end
    xlabel('x/L', 'FontSize', 16)
    
    if AllSetInfo.Replicates(exp_index2) ~= 0
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index2)),' Test Set'];
    else
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Flipped Stain Test Set'];
    end
    ylabel('Bicoid (AU)', 'FontSize', 16)
    %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
    xlim([.1, .9])
    
    title(DeltaFCAx,[num2str(DubuisTimeRange(time_index+1)), ' min into cycle 14'] , 'FontSize', 18);
    
    ylim([0, min([ymax, 150])])
    
    DeltaFCAx.YAxis.FontSize = 16;
    DeltaFCAx.XAxis.FontSize = 16;
    
    %
    
    
    DeltaFCAx.FontSize = 16;
    outpath = [AllSetsProfFigPath, filesep,'BicoidProfComps' , filesep,'TimeBin', num2str(DubuisTimeRange(time_index+1)),  '.png'];
    saveas(DeltaFCFixCorrectedFig,outpath);
    
    
    
end

%%
HunchbackProfs = NaN(71, 41, 15);
HunchbackProfs(:,:,1) = AllCompiledEmbryos{1}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);

HunchbackProfs(:,:,3) = AllCompiledEmbryos{3}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,4) = AllCompiledEmbryos{4}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,5) = AllCompiledEmbryos{5}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,6) = AllCompiledEmbryos{6}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,7) = AllCompiledEmbryos{7}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,8) = AllCompiledEmbryos{8}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,9) = AllCompiledEmbryos{9}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,11) = AllCompiledEmbryos{11}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,12) = AllCompiledEmbryos{12}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,14) = AllCompiledEmbryos{14}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
HunchbackProfs(:,:,15) = AllCompiledEmbryos{15}.BootstrappedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,5,4);
temp_indices = [4 4 4 5 5 5 3 3 3 2 2 2 1 1 1];
DubuisTimeRange = 0:70;
colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);
TempRange = 1:5;
colorsTemperature = hsv(length(TempRange)); % Colormap "jet" is another option
FractonalTempRange = TempRange/max(TempRange);
close all
for time_index = 5:5:60
    
    
    DeltaFCFixCorrectedFig =figure(time_index);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.05, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    map = colormap(colorsTemperature);
    h = colorbar;
    % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
    
    colorTitleHandle = get(h,'Title');
    titleString = 'T (ºC)';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks =  FractonalTempRange-0.5*(FractonalTempRange(2)-FractonalTempRange(1)); %Create 8 ticks from zero to 1
    h.TickLabels = {'17.5', '20.0','22.5','25.0','27.5'};%,'25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0','65.0', '70.0'};
    hold on
    for i = 1:15
        plot(APbins, HunchbackProfs(time_index+1, :, i), 'Color', colorsTemperature(temp_indices(i),:),...
            'LineWidth', 2.0)
        
        
        
    end
    
    ymax = max(HunchbackProfs(:)) *1.1;
    %plot([0 ymax], [0 ymax], 'k', 'LineWidth', 2.0);
    
    
    grid on
    hold off
    if AllSetInfo.Replicates(exp_index) ~= 0
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Test Set'];
    else
        xLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Test Set'];
    end
    xlabel('x/L', 'FontSize', 16)
    
    if AllSetInfo.Replicates(exp_index2) ~= 0
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index2)),' Test Set'];
    else
        yLabel = [num2str(AllSetInfo.Temperatures(exp_index2)), '°C Flipped Stain Test Set'];
    end
    ylabel('Hunchback (AU)', 'FontSize', 16)
    %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
    xlim([.1, .9])
    
    title(DeltaFCAx,[num2str(DubuisTimeRange(time_index+1)), ' min into cycle 14'] , 'FontSize', 18);
    
    ylim([0, min([ymax, 150])])
    
    DeltaFCAx.YAxis.FontSize = 16;
    DeltaFCAx.XAxis.FontSize = 16;
    
    %
    
    
    DeltaFCAx.FontSize = 16;
    
    outpath = [AllSetsProfFigPath, filesep,'HunchbackProfComps' , filesep,'TimeBin', num2str(DubuisTimeRange(time_index+1)),  '.png'];
    saveas(DeltaFCFixCorrectedFig,outpath);
    
end

%%
master_index = 9;
for temp = [25, 27.5, 22.5, 20, 17.5]
    for ch_index = [5]
        exp_indices = find(AllSetInfo.Temperatures == temp);
        exp_index = exp_indices(1);
        exp_index2 = exp_indices(2);
        exp_index3 = exp_indices(3);
        
        x_sample = {};
        UseTF = {};
        ys = {};
        BadTF = {};
        xfits = {};
        yfits = {};
        x_sample_control = {};
        UseTF_control = {};
        ys_control = {};
        BadTF_control = {};
        xfits_control = {};
        yfits_control = {};
        SetLabel = {};
        
        CompiledEmbryos = AllCompiledEmbryos{exp_index};
        x_sample{1} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF{1} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1});
        x_sample{1} =x_sample{1}(UseTF{1});
        if temp ~= 25
            ys{1} =  CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{1},:,ch_index,master_index);
        else
            ys{1} =  CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling(UseTF{1},:,ch_index,master_index);
        end
        BadTF{1} =  sum(isnan(ys{1}), 2).' >min( sum(isnan(ys{1}), 2).');
        x_sample{1} = x_sample{1}(~BadTF{1});
        ys{1} = ys{1}(~BadTF{1},:);
        if temp ~= 25
            xfits{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index,master_index);
        else
            xfits{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.x;
            yfits{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean(:,:,ch_index,master_index);
        end
        
        x_sample_control{1} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF_control{1} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1});
        x_sample_control{1} =x_sample_control{1}(UseTF_control{1});
        if temp ~= 25
            ys_control{1} =  CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{1},:,ch_index,master_index);
        else
            ys_control{1} =  CompiledEmbryos.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling(UseTF_control{1},:,ch_index,master_index);
        end
        BadTF_control{1} =  sum(isnan(ys_control{1}), 2).' >min( sum(isnan(ys_control{1}), 2).');
        x_sample_control{1} = x_sample_control{1}(~BadTF_control{1});
        ys_control{1} = ys_control{1}(~BadTF_control{1},:);
        if temp ~= 25
            xfits_control{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits_control{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index,master_index);
        else
            xfits_control{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.x;
            yfits_control{1} = CompiledEmbryos.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.mean(:,:,ch_index,master_index);
        end
        SetLabel{1} = AllSetInfo.SetLabels{exp_index};
        
        CompiledEmbryos2 = AllCompiledEmbryos{exp_index2};
        x_sample{2} = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF{2} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2});
        x_sample{2} =x_sample{2}(UseTF{2});
        if temp ~= 25
            ys{2} =  CompiledEmbryos2.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{2},:,ch_index,master_index);
        else
            ys{2} =  CompiledEmbryos2.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling(UseTF{2},:,ch_index,master_index);
        end
        BadTF{2} =  sum(isnan(ys{2}), 2).' >min( sum(isnan(ys{2}), 2).');
        x_sample{2} = x_sample{2}(~BadTF{2});
        ys{2} = ys{2}(~BadTF{2},:);
        if temp ~= 25
            xfits{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index,master_index);
        else
            xfits{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.x;
            yfits{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean(:,:,ch_index,master_index);
        end
        
        x_sample_control{2} = CompiledEmbryos2.DubuisEmbryoTimes;
        
        UseTF_control{2} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2});
        x_sample_control{2} =x_sample_control{2}(UseTF_control{2});
        if temp ~= 25
            ys_control{2} =  CompiledEmbryos2.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{2},:,ch_index,master_index);
        else
            ys_control{2} =  CompiledEmbryos2.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling(UseTF_control{2},:,ch_index,master_index);
        end
        BadTF_control{2}=  sum(isnan(ys_control{2}), 2).' >min( sum(isnan(ys_control{2}), 2).');
        x_sample_control{2} = x_sample_control{2}(~BadTF_control{2});
        ys_control{2} = ys_control{2}(~BadTF_control{2},:);
        if temp ~= 25
            xfits_control{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits_control{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index,master_index);
        else
            xfits_control{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.x;
            yfits_control{2} = CompiledEmbryos2.BootstrappedUnivScaledProfiles.FitSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.mean(:,:,ch_index,master_index);
        end
        
        SetLabel{2} = AllSetInfo.SetLabels{exp_index2};
        
        CompiledEmbryos3 = AllCompiledEmbryos{exp_index3};
        x_sample{3} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF{3} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3});
        x_sample{3} =x_sample{3}(UseTF{3});
        if temp ~= 25
            ys{3} =  CompiledEmbryos3.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF{3},:,ch_index,master_index);
        else
            ys{3} =  CompiledEmbryos3.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling(UseTF{3},:,ch_index,master_index);
        end
        BadTF{3} =  sum(isnan(ys{3}), 2).' >min( sum(isnan(ys{3}), 2).');
        x_sample{3} = x_sample{3}(~BadTF{3});
        ys{3} = ys{3}(~BadTF{3},:);
        if temp ~= 25
            xfits{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.TestSet.mean(:,:,ch_index,master_index);
        else
            xfits{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.x;
            yfits{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.TestSet.mean(:,:,ch_index,master_index);
            
        end
        
        x_sample_control{3} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF_control{3} = CompiledEmbryos3.ControlSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample_control{3});
        x_sample_control{3} =x_sample_control{3}(UseTF_control{3});
        if temp ~= 25
            ys_control{3} =  CompiledEmbryos3.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling(UseTF_control{3},:,ch_index,master_index);
        else
            ys_control{3} =  CompiledEmbryos3.UnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling(UseTF_control{3},:,ch_index,master_index);
        end
        BadTF_control{3}=  sum(isnan(ys_control{3}), 2).' >min( sum(isnan(ys_control{3}), 2).');
        x_sample_control{3} = x_sample_control{3}(~BadTF_control{3});
        ys_control{3} = ys_control{3}(~BadTF_control{3},:);
        if temp ~= 25
            xfits_control{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.x;
            yfits_control{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.ControlScaling.ControlSet.mean(:,:,ch_index,master_index);
        else
            xfits_control{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.x;
            yfits_control{3} = CompiledEmbryos3.BootstrappedUnivScaledProfiles.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles.TestScaling.ControlSet.mean(:,:,ch_index,master_index);
            
        end
        
        SetLabel{3} = AllSetInfo.SetLabels{exp_index3};
        
        for APindex =5:4:37
            ymax = max([max( yfits_control{1}(:,APindex)), max( yfits_control{2}(:,APindex)), max( yfits_control{3}(:,APindex)),...
                max( yfits{1}(:,APindex)), max( yfits{2}(:,APindex)), max( yfits{3}(:,APindex)),...
                max(ys_control{1}(:,APindex)), max(ys_control{2}(:,APindex)), max(ys_control{3}(:,APindex)),...
                max(ys{1}(:,APindex)),max(ys{2}(:,APindex)),max(ys{3}(:,APindex))])*1.1;
            DeltaFCAx = cell(2,3);
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .8]);
            set(gcf,'color','w');
            for sub_index = 1:3
                DeltaFCAx{1,sub_index} = subplot(2, 3, sub_index);
                scatter(x_sample{sub_index}, ys{sub_index}(:,APindex),...
                    100, 'MarkerFaceColor',gap_colors(ch_index,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(ch_index,:));% [.7 .7 .7]);
                hold on
                plot(xfits{sub_index}, yfits{sub_index}(:,APindex), 'k', 'LineWidth', 2.0);
                %ymax = max([max(yfits{sub_index}(:, APindex)),max(ys{sub_index}(:,APindex))]) *1.1;
                
                grid on
                hold off
                xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
                %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
                xlim([0, 60])
                ylim([0, ymax])
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
                
                %ylim([0, min([ymax, 100])])
                
                
                DeltaFCAx{1,sub_index}.YAxis.FontSize = 14;
                DeltaFCAx{1,sub_index}.XAxis.FontSize = 14;
                
                if AllSetInfo.Replicates(exp_indices(sub_index)) ~= 0
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_indices(sub_index))),' Test Set'];
                else
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Flipped Stain Test Set'];
                end
                title(DeltaFCAx{1,sub_index}, [PlotLabel], 'FontSize', 14)
                %
                
                
                DeltaFCAx{1,sub_index}.FontSize = 14;
                
                
                
                DeltaFCAx{2,sub_index} = subplot(2, 3, sub_index+3);
                scatter(x_sample_control{sub_index}, ys_control{sub_index}(:,APindex),...
                    100, 'MarkerFaceColor',gap_colors(ch_index,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(ch_index,:));% [.7 .7 .7]);
                hold on
                plot(xfits_control{sub_index}, yfits_control{sub_index}(:,APindex), 'k', 'LineWidth', 2.0);
                %ymax = max([max(yfits_control{sub_index}(:, APindex)),max(ys_control{sub_index}(:,APindex))]) *1.1;
                
                grid on
                hold off
                xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
                
                
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
                
                ylim([0, ymax])
                xlim([0 60])
                
                DeltaFCAx{2, sub_index}.YAxis.FontSize = 14;
                DeltaFCAx{2, sub_index}.XAxis.FontSize = 14;
                
                if AllSetInfo.Replicates(exp_index) ~= 0
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_index)),' Control Set'];
                else
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_index)), '°C Flipped Stain Control Set'];
                end
                title(DeltaFCAx{2, sub_index}, [PlotLabel], 'FontSize', 14)
                %
                
                
                DeltaFCAx{2, sub_index}.FontSize = 14;
            end
            sgtitle(['AP: ', num2str(round(APbins(APindex), 2))], 'FontSize', 16)
            outpath2 = [AllSetsProfFigPath, filesep, 'BootstrappedProfiles', filesep, 'FitSlideRescaledAvgAP', filesep,...
                ChannelNames{ch_index}, '_T', strrep(num2str(temp), '.', '_'), 'C_AP',num2str(APindex), '.png' ]
            saveas(DeltaFCFixCorrectedFig,outpath2);
        end
    end
end
%% Plot Bootstrapped unscaled profiles


for temp = [25, 27.5, 22.5, 20, 17.5]
    for ch_index = [3]
        exp_indices = find(AllSetInfo.Temperatures == temp);
        exp_index = exp_indices(1);
        exp_index2 = exp_indices(2);
        exp_index3 = exp_indices(3);
        
        x_sample = {};
        UseTF = {};
        ys = {};
        BadTF = {};
        xfits = {};
        yfits = {};
        x_sample_control = {};
        UseTF_control = {};
        ys_control = {};
        BadTF_control = {};
        xfits_control = {};
        yfits_control = {};
        SetLabel = {};
        
        CompiledEmbryos = AllCompiledEmbryos{exp_index};
        x_sample{1} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF{1} = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample{1});
        x_sample{1} =x_sample{1}(UseTF{1});
        ys{1} =  CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(UseTF{1},:,ch_index);
        BadTF{1} =  sum(isnan(ys{1}), 2).' >min( sum(isnan(ys{1}), 2).');
        x_sample{1} = x_sample{1}(~BadTF{1});
        ys{1} = ys{1}(~BadTF{1},:);
        xfits{1} = CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
        yfits{1} = CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index);
        
        
        x_sample_control{1} = CompiledEmbryos.DubuisEmbryoTimes;
        UseTF_control{1} = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(x_sample_control{1});
        x_sample_control{1} =x_sample_control{1}(UseTF_control{1});
        ys_control{1} =  CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(UseTF_control{1},:,ch_index);
        BadTF_control{1} =  sum(isnan(ys_control{1}), 2).' >min( sum(isnan(ys_control{1}), 2).');
        x_sample_control{1} = x_sample_control{1}(~BadTF_control{1});
        ys_control{1} = ys_control{1}(~BadTF_control{1},:);
        xfits_control{1} =  CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
        yfits_control{1} = CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.mean(:,:,ch_index);
        SetLabel{1} = AllSetInfo.SetLabels{exp_index};
        
        CompiledEmbryos2 = AllCompiledEmbryos{exp_index2};
        x_sample{2} = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF{2} = CompiledEmbryos2.TestSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample{2});
        x_sample{2} =x_sample{2}(UseTF{2});
        ys{2} =  CompiledEmbryos2.FitSlideRescaledDorsalAvgAPProfiles(UseTF{2},:,ch_index);
        BadTF{2} =  sum(isnan(ys{2}), 2).' >min( sum(isnan(ys{2}), 2).');
        x_sample{2} = x_sample{2}(~BadTF{2});
        ys{2} = ys{2}(~BadTF{2},:);
        xfits{2} = CompiledEmbryos2.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
        yfits{2} = CompiledEmbryos2.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index);
        
        x_sample_control{2} = CompiledEmbryos2.DubuisEmbryoTimes;
        UseTF_control{2} = CompiledEmbryos2.ControlSetEmbryos & CompiledEmbryos2.IsNC14 & ~isnan(x_sample_control{2});
        x_sample_control{2} =x_sample_control{2}(UseTF_control{2});
        ys_control{2} =  CompiledEmbryos2.FitSlideRescaledDorsalAvgAPProfiles(UseTF_control{2},:,ch_index);
        BadTF_control{2}=  sum(isnan(ys_control{2}), 2).' >min( sum(isnan(ys_control{2}), 2).');
        x_sample_control{2} = x_sample_control{2}(~BadTF_control{2});
        ys_control{2} = ys_control{2}(~BadTF_control{2},:);
        xfits_control{2} =  CompiledEmbryos2.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
        yfits_control{2} = CompiledEmbryos2.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.mean(:,:,ch_index);
        SetLabel{2} = AllSetInfo.SetLabels{exp_index2};
        
        CompiledEmbryos3 = AllCompiledEmbryos{exp_index3};
        x_sample{3} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF{3} = CompiledEmbryos3.TestSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample{3});
        x_sample{3} =x_sample{3}(UseTF{3});
        ys{3} =  CompiledEmbryos3.FitSlideRescaledDorsalAvgAPProfiles(UseTF{3},:,ch_index);
        BadTF{3} =  sum(isnan(ys{3}), 2).' >min( sum(isnan(ys{3}), 2).');
        x_sample{3} = x_sample{3}(~BadTF{3});
        ys{3} = ys{3}(~BadTF{3},:);
        xfits{3} = CompiledEmbryos3.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
        yfits{3} = CompiledEmbryos3.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index);
        
        x_sample_control{3} = CompiledEmbryos3.DubuisEmbryoTimes;
        UseTF_control{3} = CompiledEmbryos3.ControlSetEmbryos & CompiledEmbryos3.IsNC14 & ~isnan(x_sample_control{3});
        x_sample_control{3} =x_sample_control{3}(UseTF_control{3});
        ys_control{3} =  CompiledEmbryos3.FitSlideRescaledDorsalAvgAPProfiles(UseTF_control{3},:,ch_index);
        BadTF_control{3}=  sum(isnan(ys_control{3}), 2).' >min( sum(isnan(ys_control{3}), 2).');
        x_sample_control{3} = x_sample_control{3}(~BadTF_control{3});
        ys_control{3} = ys_control{3}(~BadTF_control{3},:);
        xfits_control{3} =  CompiledEmbryos3.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x;
        yfits_control{3} = CompiledEmbryos3.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.mean(:,:,ch_index);
        
        SetLabel{3} = AllSetInfo.SetLabels{exp_index3};
        
        for APindex =5:4:37
            ymax = max([max( yfits_control{1}(:,APindex)), max( yfits_control{2}(:,APindex)), max( yfits_control{3}(:,APindex)),...
                max( yfits{1}(:,APindex)), max( yfits{2}(:,APindex)), max( yfits{3}(:,APindex)),...
                max(ys_control{1}(:,APindex)), max(ys_control{2}(:,APindex)), max(ys_control{3}(:,APindex)),...
                max(ys{1}(:,APindex)),max(ys{2}(:,APindex)),max(ys{3}(:,APindex))])*1.1;
            
            DeltaFCAx = cell(2,3);
            close all
            DeltaFCFixCorrectedFig =figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .8]);
            set(gcf,'color','w');
            for sub_index = 1:3
                ymax = max([max(ys_control{sub_index}(:, APindex)), max(yfits_control{sub_index}(:, APindex)),...
                    max(yfits{sub_index}(:, APindex)),  max(ys{sub_index}(:, APindex))])*1.1;
                
                DeltaFCAx{1,sub_index} = subplot(2, 3, sub_index);
                scatter(x_sample{sub_index}, ys{sub_index}(:,APindex),...
                    100, 'MarkerFaceColor',gap_colors(ch_index,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(ch_index,:));% [.7 .7 .7]);
                hold on
                plot(xfits{sub_index}, yfits{sub_index}(:,APindex), 'k', 'LineWidth', 2.0);
                %ymax = max([max(yfits{sub_index}(:, APindex)),max(ys{sub_index}(:,APindex))]) *1.1;
                
                grid on
                hold off
                xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
                %xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
                xlim([0, 60])
                ylim([0, ymax])
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
                
                %ylim([0, min([ymax, 100])])
                
                
                DeltaFCAx{1,sub_index}.YAxis.FontSize = 14;
                DeltaFCAx{1,sub_index}.XAxis.FontSize = 14;
                
                if AllSetInfo.Replicates(exp_indices(sub_index)) ~= 0
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_indices(sub_index))),' Test Set'];
                else
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Flipped Stain Test Set'];
                end
                title(DeltaFCAx{1,sub_index}, [PlotLabel], 'FontSize', 14)
                %
                
                
                DeltaFCAx{1,sub_index}.FontSize = 14;
                
                
                
                DeltaFCAx{2,sub_index} = subplot(2, 3, sub_index+3);
                scatter(x_sample_control{sub_index}, ys_control{sub_index}(:,APindex),...
                    100, 'MarkerFaceColor',gap_colors(ch_index,:),...% [.7 .7 .7],...
                    'MarkerEdgeColor',gap_colors(ch_index,:));% [.7 .7 .7]);
                hold on
                plot(xfits_control{sub_index}, yfits_control{sub_index}(:,APindex), 'k', 'LineWidth', 2.0);
                %ymax = max([max(yfits_control{sub_index}(:, APindex)),max(ys_control{sub_index}(:,APindex))]) *1.1;
                
                grid on
                hold off
                xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 14)
                
                ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 14)
                
                ylim([0, ymax])
                xlim([0 60])
                
                DeltaFCAx{2, sub_index}.YAxis.FontSize = 14;
                DeltaFCAx{2, sub_index}.XAxis.FontSize = 14;
                
                if AllSetInfo.Replicates(exp_indices(sub_index)) ~= 0
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Replicate ', num2str(AllSetInfo.Replicates(exp_indices(sub_index))),' Control Set'];
                else
                    PlotLabel = [num2str(AllSetInfo.Temperatures(exp_indices(sub_index))), '°C Flipped Stain Control Set'];
                end
                title(DeltaFCAx{2, sub_index}, [PlotLabel], 'FontSize', 14)
                %
                
                
                DeltaFCAx{2, sub_index}.FontSize = 14;
            end
            sgtitle(['AP: ', num2str(round(APbins(APindex), 2))], 'FontSize', 16)
            outpath2 = [AllSetsProfFigPath, filesep, 'BootstrappedProfiles2', filesep, 'FitSlideRescaledAvgAP', filesep,...
                ChannelNames{ch_index}, '_T', strrep(num2str(temp), '.', '_'), 'C_AP',num2str(APindex), '.png' ]
            saveas(DeltaFCFixCorrectedFig,outpath2);
            
        end
    end
end