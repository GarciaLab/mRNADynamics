for i = 2:15
    SetLabel = AllSetInfo.SetLabels{i};
    for ch_index = [3] % 2:5
        if ch_index == 3
            RefBins =  13;
        else
            RefBins =  13;
        end
        for RefBin =RefBins
       
            close all
            DeltaFCFixCorrectedFig = figure(1);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
%             scatter(AllCompiledEmbryos{i}.DubuisEmbryoTimes(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF),...
%                 AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF,RefBin,ch_index).',...
%                 100, 'MarkerFaceColor', [.7 .7 .7],...
%                 'MarkerEdgeColor', [.7 .7 .7]); 
            
            scatter(AllCompiledEmbryos{i}.FixCorrectedDeltaFC_um.mean(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF),...
                AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF,RefBin,ch_index).',...
                100, 'MarkerFaceColor', [.7 .7 .7],...
                'MarkerEdgeColor', [.7 .7 .7]); 
            
            
            x_unsorted = AllCompiledEmbryos{i}.DubuisEmbryoTimes(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF);
            p_unsorted = AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(AllCompiledEmbryos{i}.AllDubuisValidProfilesControlTF,RefBin,ch_index).';
            [x, sort_order] = sort(x_unsorted);
            p = p_unsorted(sort_order);
            %pdsg =  smoothdata(p,'sgolay','degree',2,'samplepoints',x);
            hold on 
%             MeanProfile = AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.Control.mean(:,RefBin,ch_index).';
%             plot(AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.x, MeanProfile, 'k', 'LineWidth', 2.0);
%             hold on 
%             counts_1sigma = AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.Control.count(:,RefBin,ch_index).';
%             counts_1sigma = AllCompiledEmbryos{i}.DubuisTimesSmoothedProfiles.TestCounts.counts_within_1sigma;
%             plot(x, p, 'r-', 'LineWidth', 2.0);
%             scatter(AllCompiledEmbryos{i}.WindowedProfiles.DubuisTime.x(counts_1sigma >= 3), MeanProfile(counts_1sigma >= 3),...
%                 75, 'MarkerFaceColor', gap_colors(ch_index,:), 'MarkerEdgeColor', gap_colors(ch_index,:));
%             
%             
            MeanProfile = AllCompiledEmbryos{i}.FixCorrectedSmoothedAvgAPProfiles.Control(:,RefBin,ch_index).';
            plot(AllCompiledEmbryos{i}.FixCorrectedSmoothedDeltaFCs, MeanProfile, 'k', 'LineWidth', 2.0);
            hold on 
%             counts_1sigma = AllCompiledEmbryos{i}.FixCorrectedDeltaFC_um(:,RefBin,ch_index).';
%             %counts_1sigma = AllCompiledEmbryos{i}.DubuisTimesSmoothedProfiles.TestCounts.counts_within_1sigma;
%             %plot(x, p, 'r-', 'LineWidth', 2.0);
%             scatter(AllCompiledEmbryos{i}.WindowedProfiles.DeltaFC.x, MeanProfile(counts_1sigma >= 3),...
%                 75, 'MarkerFaceColor', gap_colors(ch_index,:), 'MarkerEdgeColor', gap_colors(ch_index,:));
%             
%             
            
            

          
            ymax = max(max(AllCompiledEmbryos{i}.SlideRescaledDorsalAvgAPProfiles(:,RefBin,ch_index)))*1.1;
            %     plot([MinTKeyRegion, MinTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
            %      plot([MaxTKeyRegion, MaxTKeyRegion], [0, ymax], 'k--','LineWidth', 2.0);
            %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
            %     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
            grid on
            hold off
            xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
            xlabel('\delta_{FC} (\mum)', 'FontSize', 16)
            xlim([0, 45])
            
            ylabel([ChannelNames{ch_index}, ' (AU)'], 'FontSize', 16)
            
            ylim([0, ymax])
            
            
            DeltaFCAx.YAxis.FontSize = 16;
            DeltaFCAx.XAxis.FontSize = 16;
            
            if AllSetInfo.Replicates(i) ~= 0
            PlotLabel = [num2str(AllSetInfo.Temperatures(i)), '°C Replicate ', num2str(AllSetInfo.Replicates(i)),' Control Set'];
            else
                 PlotLabel = [num2str(AllSetInfo.Temperatures(i)), '°C Flipped Stain Control Set'];
            end
            title(DeltaFCAx, [PlotLabel,' AP: ',num2str(APbins(RefBin)), ' (Test)'], 'FontSize', 18)
            %
            
            
            DeltaFCAx.FontSize = 16;
            %outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test','_'];
            outpath = [AllSetsProfFigPath, filesep,'SetComps' , filesep, SetLabel, '_Control_',ChannelNames{ch_index}, '_DeltaFCWindowed_AP',num2str(RefBin),'.png'];
            saveas(DeltaFCFixCorrectedFig,outpath);
        end
    end
end
%%
for i = 11:15
    SetLabel = AllSetInfo.SetLabels{i};
    for ch_index = [3, 5] % 2:5
        if ch_index == 3
            RefBins =  13'%9:21;
        else
            RefBins =  13;%[7:21 29:37];
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


%%
%%
if ~isdir([AllSetsProfFigPath, filesep, 'FlippedSetComps'])
    mkdir([AllSetsProfFigPath, filesep, 'FlippedSetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate 1'];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_FlippedRep1_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(2); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate 2'];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_FlippedRep2_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate 1'];
ylab = ['Control for ',num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_FlippedRep1_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(2); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate 2'];
ylab = ['Control for ',num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_FlippedRep2_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%% Use DeltaFC Profiles
if ~isdir([AllSetsProfFigPath, 'DeltaFCSetComps'])
    mkdir([AllSetsProfFigPath, 'DeltaFCSetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    exp2 = ReplicateIndices(2);
for ch_index = [3:5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
    end

    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, ' Test Sets with windowed Delta FC profiles'], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'DeltaFCSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_StandardReplicates_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    exp2 = ReplicateIndices(2);
for ch_index = 3:5
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, ' Control Sets with windowed Delta FC profiles'], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'DeltaFCSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_StandardReplicates_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
if ~isdir([AllSetsProfFigPath, filesep, 'FlippedSetComps'])
    mkdir([AllSetsProfFigPath, filesep, 'FlippedSetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate 1'];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_FlippedRep1_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(2); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate 2'];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_FlippedRep2_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate 1'];
ylab = ['Control for ',num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_FlippedRep1_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(2); 
    FlippedIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 1);
    exp2 = FlippedIndices(1);
for ch_index = [3, 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 3 & ProfCounts2 >= 3);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 3 & ProfCounts2 >= 3);
if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate 2'];
ylab = ['Control for ',num2str(unique_temperatures(temp_index)), 'ºC Flipped Condition '];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'FlippedSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_FlippedRep2_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
%%
%%
if ~isdir([AllSetsProfFigPath, filesep, 'DeltaFCSetCompsTo17_5CRep2'])
    mkdir([AllSetsProfFigPath, filesep, 'DeltaFCSetCompsTo17_5CRep2'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for exp2 = [1:13 15]%length(unique_temperatures)
   
    exp1 = 14;
for ch_index = 3:5
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(AllSetInfo.Temperatures(exp1)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = ['Control for ', num2str(AllSetInfo.Temperatures(exp2)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, ' Control Sets with windowed Delta FC profiles', AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'DeltaFCSetCompsTo17_5CRep2',filesep,'T', strrep(num2str(num2str(AllSetInfo.Temperatures(exp2))),'.','_'),'C_Replicate',num2str(AllSetInfo.Replicates(exp2)),'control_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
%%
%%
if ~isdir([AllSetsProfFigPath, filesep, 'DubuisTimeSetCompsTo17_5CRep2'])
    mkdir([AllSetsProfFigPath, filesep, 'DubuisTimeSetCompsTo17_5CRep2'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
for exp2 = [1:13 15]%length(unique_temperatures)
   
    exp1 = 14;
for ch_index = 3:5
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.mean(:,:, ch_index);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DubuisTimesWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        AddTitleStr = ' (Low Counts)';
    else
        AddTitleStr = '';
    end
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);

    scatter(PlotProfSet1,PlotProfSet2,...
        100, 'MarkerFaceColor', gap_colors(ch_index,:),...
        'MarkerEdgeColor', gap_colors(ch_index,:));
    hold on 


    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(AllSetInfo.Temperatures(exp1)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = ['Control for ', num2str(AllSetInfo.Temperatures(exp2)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];

xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, [ChannelNames{ch_index}, ' Control Sets with windowed Dubuis Time profiles', AddTitleStr], 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'DubuisTimeSetCompsTo17_5CRep2',filesep,'T', strrep(num2str(num2str(AllSetInfo.Temperatures(exp2))),'.','_'),'C_Replicate',num2str(AllSetInfo.Replicates(exp2)),'control_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
%% Comps Control Set Replicates for each temperature with spots colored by AP position 
if ~isdir([AllSetsProfFigPath, 'DeltaFCSetComps'])
    mkdir([AllSetsProfFigPath, 'DeltaFCSetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
APRange = 0.1:0.025:0.9;
colors = hsv(length(APRange)); % Colormap "jet" is another option
MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};
FractionalAPRange = (1:length(APRange))/length(APRange);

y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);

for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    exp2 = ReplicateIndices(2);
for ch_index = [3 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    APmat1 = repmat(1:41, size(ProfSet1,1), 1);
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    APmat2 = repmat(1:41, size(ProfSet2,1), 1);
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    ProfAPs1 = APmat1(:)-4;
    ProfAPs2 = APmat2(:)-4;
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotAPs1 = ProfAPs1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotAPs2 = ProfAPs2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotAPs1 = ProfAPs1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotAPs2 = ProfAPs2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
    end
        
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = 'x/L';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalAPRange(1:4:33); %Create 8 ticks from zero to 1
h.TickLabels = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'};%'20.0','25.0', '30.0', '35.0', '40.0', '45.0'} ;
hold on

    for plot_index = 1:length(PlotProfSet1)
    scatter(PlotProfSet1(plot_index),PlotProfSet2(plot_index),...
        75, 'MarkerFaceColor', colors(PlotAPs1(plot_index),:),...
        'MarkerEdgeColor', colors(PlotAPs1(plot_index),:));
    hold on 
    end

    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k')
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];
xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, ChannelNames{ch_index}, 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath, filesep, 'DeltaFCSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_StandardReplicates_APcolor_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%% Comps Control Set Replicates for each temperature with spots colored by DeltaFC position 
if ~isdir([AllSetsProfFigPath, 'DeltaFCSetComps'])
    mkdir([AllSetsProfFigPath, 'DeltaFCSetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
DeltaFCRange = 0:45;
colorsDeltaFC = hsv(length(DeltaFCRange)); % Colormap "jet" is another option
FractionalDeltaFCRange = DeltaFCRange/max(DeltaFCRange);

y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);

for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    exp2 = ReplicateIndices(2);
for ch_index = [3 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    DeltaFCmat1 = repmat(AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.x.',  1, size(ProfSet1,2));
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.mean(:,:, ch_index);
    DeltaFCmat2 = repmat(AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.x.',  1, size(ProfSet2,2));
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Control.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    DeltaFCs1 = DeltaFCmat1(:);
    DeltaFCs2 = DeltaFCmat2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotDeltaFCs1 = DeltaFCs1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotDeltaFCs2 = DeltaFCs2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotDeltaFCs1 = DeltaFCs1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotDeltaFCs2 = DeltaFCs2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
    end
        
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = '\delta_{FC} (\mum)';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalDeltaFCRange(6:5:41)+0.5*(FractionalDeltaFCRange(2)-FractionalDeltaFCRange(1)); %Create 8 ticks from zero to 1
h.TickLabels = {'5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0'} ;
hold on

    for plot_index = 1:length(PlotProfSet1)
    scatter(PlotProfSet1(plot_index),PlotProfSet2(plot_index),...
        75, 'MarkerFaceColor', colorsDeltaFC(PlotDeltaFCs1(plot_index)+1,:),...
        'MarkerEdgeColor', colorsDeltaFC(PlotDeltaFCs1(plot_index)+1,:));
    hold on 
    end

    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k');
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = ['Control for ', num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];
xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, ChannelNames{ch_index}, 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath,filesep, 'HernanChat' filesep, 'DeltaFCSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ccontrol_StandardReplicates_DeltaFCcolor_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end

%% Comps Test Set Replicates for each temperature with spots colored by DeltaFC position 
if ~isdir([AllSetsProfFigPath, filesep 'HernanChat', filesep, 'DeltaFCSetComps'])
    mkdir([AllSetsProfFigPath, filesep 'HernanChat', filesep, 'DeltaFCSetComps'])
end
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
DeltaFCRange = 0:45;
colorsDeltaFC = hsv(length(DeltaFCRange)); % Colormap "jet" is another option
FractionalDeltaFCRange = DeltaFCRange/max(DeltaFCRange);

y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);

for temp_index = 1:length(unique_temperatures)
    ReplicateIndices = find(AllSetInfo.Temperatures == unique_temperatures(temp_index) & AllSetInfo.Flipped == 0);
    exp1 = ReplicateIndices(1); 
    exp2 = ReplicateIndices(2);
for ch_index = [3 5]
    ProfSet1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.mean(:,:, ch_index);
    DeltaFCmat1 = repmat(AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.x.',  1, size(ProfSet1,2));
    ProfSet2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.mean(:,:, ch_index);
    DeltaFCmat2 = repmat(AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.x.',  1, size(ProfSet2,2));
    ProfCounts1 = AllCompiledEmbryos{exp1}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.count(:,:, ch_index);
    ProfCounts2 = AllCompiledEmbryos{exp2}.NormalizedProfiles.DeltaFCWindowedAvgAP.Test.count(:,:, ch_index);
    ProfSet1 = ProfSet1(:);
    ProfSet2 = ProfSet2(:);
    ProfCounts1 = ProfCounts1(:);
    ProfCounts2 = ProfCounts2(:);
    DeltaFCs1 = DeltaFCmat1(:);
    DeltaFCs2 = DeltaFCmat2(:);
    PlotProfSet1 = ProfSet1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotProfSet2 = ProfSet2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotDeltaFCs1 = DeltaFCs1(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    PlotDeltaFCs2 = DeltaFCs2(ProfCounts1 >= 2 & ProfCounts2 >= 2);
    if isempty(PlotProfSet1)
        PlotProfSet1 = ProfSet1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotProfSet2 = ProfSet2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotDeltaFCs1 = DeltaFCs1(ProfCounts1 >= 1 & ProfCounts2 >= 1);
        PlotDeltaFCs2 = DeltaFCs2(ProfCounts1 >= 1 & ProfCounts2 >= 1);
    end
        
    MaxValues = NaN(1, NumSets);
    close all
    DeltaFCFixCorrectedFig = figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    map = colormap(colors);
h = colorbar;
% %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
hold off
colorTitleHandle = get(h,'Title');
titleString = '\delta_{FC} (\mum)';
set(colorTitleHandle ,'String',titleString);
h.Ticks =  FractionalDeltaFCRange(6:5:41)+0.5*(FractionalDeltaFCRange(2)-FractionalDeltaFCRange(1)); %Create 8 ticks from zero to 1
h.TickLabels = {'5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0'} ;
hold on

    for plot_index = 1:length(PlotProfSet1)
    scatter(PlotProfSet1(plot_index),PlotProfSet2(plot_index),...
        75, 'MarkerFaceColor', colorsDeltaFC(PlotDeltaFCs1(plot_index)+1,:),...
        'MarkerEdgeColor', colorsDeltaFC(PlotDeltaFCs1(plot_index)+1,:));
    hold on 
    end

    ymax = ceil(max([PlotProfSet1 ; PlotProfSet2])/.5)*.5;
    plot([0, ymax], [0, ymax], 'k');
     %ymax = max(max(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,:,ch_index)))*1.1;
%     plot([y_positions(ch_index), y_positions(ch_index)], [0, ymax], 'k:','LineWidth', 2.0);
    grid on 
    hold off
    %xlabel('Time (Dubuis) into cycle 14 (min)', 'FontSize', 16)
xlim([0, ymax])
    xlab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp1))];
ylab = [num2str(unique_temperatures(temp_index)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(exp2))];
xlabel(xlab, 'FontSize', 16)
ylabel(ylab, 'FontSize', 16)

ylim([0, ymax])


DeltaFCAx.YAxis.FontSize = 16;
DeltaFCAx.XAxis.FontSize = 16;

title(DeltaFCAx, ChannelNames{ch_index}, 'FontSize', 18)
%


DeltaFCAx.FontSize = 16;

outpath = [AllSetsProfFigPath,filesep, 'HernanChat' filesep, 'DeltaFCSetComps',filesep,'T', strrep(num2str( num2str(unique_temperatures(temp_index))),'.','_'),'Ctest_StandardReplicates_DeltaFCcolor_',ChannelNames{ch_index}, '.png'];
saveas(DeltaFCFixCorrectedFig,outpath);
end
end
