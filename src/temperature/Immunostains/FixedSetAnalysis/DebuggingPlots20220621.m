for i = 1:15
    SetLabel = AllSetInfo.SetLabels{i};
    for ch_index = [3, 5] % 2:5
        if ch_index == 3
            RefBins =  7:4:21;
        else
            RefBins =  [7:4:21 31:4:37];
        end
        for RefBin =RefBins
       
            close all
            DeltaFCFixCorrectedFig = figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            scatter(AllCompiledEmbryos{i}.DubuisEmbryoTimes(AllCompiledEmbryos{i}.AllDubuisValidProfilesTestTF),...
                AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(AllCompiledEmbryos{i}.AllDubuisValidProfilesTestTF,RefBin,ch_index).',...
                100, 'MarkerFaceColor', [.7 .7 .7],...
                'MarkerEdgeColor', [.7 .7 .7]); 
            
            
            x_unsorted = AllCompiledEmbryos{i}.DubuisEmbryoTimes(AllCompiledEmbryos{i}.AllDubuisValidProfilesTestTF);
            p_unsorted = AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(AllCompiledEmbryos{i}.AllDubuisValidProfilesTestTF,RefBin,ch_index).';
            [x, sort_order] = sort(x_unsorted);
            p = p_unsorted(sort_order);
            %pdsg =  smoothdata(p,'sgolay','degree',2,'samplepoints',x);
            hold on 
            MeanProfile = AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.Test.mean(:,RefBin,ch_index).';
            plot(AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.x, MeanProfile, 'k', 'LineWidth', 2.0);
            hold on 
            counts_1sigma = AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.Test.count(:,RefBin,ch_index).';
            %counts_1sigma = AllCompiledEmbryos{i}.DubuisTimesSmoothedProfiles.TestCounts.counts_within_1sigma;
            %plot(x, p, 'r-', 'LineWidth', 2.0);
            scatter(AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.x(counts_1sigma >= 3), MeanProfile(counts_1sigma >= 3),...
                75, 'MarkerFaceColor', gap_colors(ch_index,:), 'MarkerEdgeColor', gap_colors(ch_index,:));
            
            

          
            ymax = max(max(AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(:,RefBin,ch_index)))*1.1;
            %     plot([MinTKeyRegion, MinTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
            %      plot([MaxTKeyRegion, MaxTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
            %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
            %     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
            grid on
            hold off
            xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
            xlim([0, 65])
            
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)
            
            ylim([0, ymax])
            
            
            DeltaFCAx.YAxis.FontSize = 16;
            DeltaFCAx.XAxis.FontSize = 16;
            
            if AllSetInfo.Replicates(i) ~= 0
            PlotLabel = [num2str(AllSetInfo.Temperatures(i)), '°C Replicate ', num2str(AllSetInfo.Replicates(i)),' Test Set'];
            else
                 PlotLabel = [num2str(AllSetInfo.Temperatures(i)), '°C Flipped Stain Test Set'];
            end
            title(DeltaFCAx, [PlotLabel,' AP: ',num2str(APbins(RefBin)), ' (Test)'], 'FontSize', 18)
            %
            
            
            DeltaFCAx.FontSize = 16;
            %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
            outpath = [AllSetsProfFigPath, filesep,'SetComps' , filesep, SetLabel, '_Test_',ChannelNames{ch_index}, '_TimeWindowed_AP',num2str(RefBin),'.png'];
            saveas(DeltaFCFixCorrectedFig,outpath);
        end
    end
end
%%
for i = 11:15
    SetLabel = AllSetInfo.SetLabels{i};
    for ch_index = [3, 5] % 2:5
        if ch_index == 3
            RefBins =  9:21;
        else
            RefBins =  [7:21 29:37];
        end
        for RefBin =RefBins
       
            close all
            DeltaFCFixCorrectedFig = figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            scatter(AllCompiledEmbryos{i}.DubuisEmbryoTimes(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF),...
                AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF,RefBin,ch_index).',...
                100, 'MarkerFaceColor', [.7 .7 .7],...
                'MarkerEdgeColor', [.7 .7 .7]); 
            
            
            x_unsorted = AllCompiledEmbryos{i}.DubuisEmbryoTimes(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF);
            p_unsorted = AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF,RefBin,ch_index).';
            [x, sort_order] = sort(x_unsorted);
            p = p_unsorted(sort_order);
            %pdsg =  smoothdata(p,'sgolay','degree',2,'samplepoints',x);
            hold on 
            MeanProfile = AllCompiledEmbryos{i}.DubuisTimeSmoothedAvgAPProfiles.Control(:,RefBin,ch_index).';
            plot(AllCompiledEmbryos{i}.DubuisSmoothedTimes, MeanProfile, 'k', 'LineWidth', 2.0);
            hold on 
            counts_1sigma = AllCompiledEmbryos{i}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_1sigma;
            %counts_1sigma = AllCompiledEmbryos{i}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_halfsigma;
            %plot(x, p, 'r-', 'LineWidth', 2.0);
            scatter(AllCompiledEmbryos{i}.DubuisSmoothedTimes(counts_1sigma >= 5), MeanProfile(counts_1sigma >= 5),...
                75, 'MarkerFaceColor', gap_colors(ch_index,:), 'MarkerEdgeColor', gap_colors(ch_index,:));
            
            ymax = max(max(AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(:,:,ch_index)))*1.1;
      
            grid on
            hold off
            xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
            xlim([0, 65])
            
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)
            
            ylim([0, ymax])
            
            
            DeltaFCAx.YAxis.FontSize = 16;
            DeltaFCAx.XAxis.FontSize = 16;
            
            if AllSetInfo.Replicates(i) ~= 0
            PlotLabel = [num2str(AllSetInfo.Temperatures(i)), '°C Replicate ', num2str(AllSetInfo.Replicates(i)),' Control Set'];
            else
                 PlotLabel = [num2str(AllSetInfo.Temperatures(i)), '°C Flipped Stain Control Set'];
            end
            title(DeltaFCAx, [PlotLabel,' AP: ',num2str(APbins(RefBin)), ' (Control)'], 'FontSize', 18)
            %
            
            
            DeltaFCAx.FontSize = 16;
            %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
            outpath = [AllSetsProfFigPath, filesep,'SetComps' , filesep, SetLabel, '_Control_',ChannelNames{ch_index}, '_TimeSmoothed_AP',num2str(RefBin),'.png'];
            saveas(DeltaFCFixCorrectedFig,outpath);
        end
    end
end