version = 3;
AllCompiledPath = ['S:/Gabriella/Dropbox/ProteinProfiles/V',num2str(version),'Profiles/'];
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/V3/Figures/';
load([AllCompiledPath, 'AllCompiledEmbryos.Mat'], 'AllCompiledEmbryos','Universal_Imins','Universal_Imaxs','Imins','Imaxs');
AllSetInfo = GetFixedSetPrefixInfo;

NumSets = length(AllSetInfo.Temperatures);

NChannels = 5;

exp_index = find(AllSetInfo.Temperatures == 20 & AllSetInfo.Replicates == 1);


APbins = 0:0.025:1;
NumAPbins = length(APbins);

NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);

for i = 1:NumSets
    AllCompiledEmbryos{i} = AddBinnedProfiles( AllCompiledEmbryos{i});
    AllCompiledEmbryos{i} = AddSmoothedProfiles( AllCompiledEmbryos{i});
end

%%
ChannelNames = {'', 'Hoechst', 'Bicoid', 'Knirps', 'Hunchback'};

gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 177 32 ;
    0 0 0 ]/255;%0.9290 0.6940 0.1250]/255;
y_positions = [0, 0.5, 0.3, 0.6, 0.8];
DubuisTimeRange = 0:70;
colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);

y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
BinsToFit = zeros(NChannels, NumAPbins, 'logical');
BinsToFit(2,9:33) = true;
BinsToFit(3,5:13) = true;
BinsToFit(4,22:37) = true;
BinsToFit(5,5:21) = true;

%%
close all
gap_colors = [0 0 0 ;
    66 87 168;
    100 149 93;
    212 84 66;
    237 176 32 ]/255;%0.9290 0.6940 0.1250]/255;
DubuisTimeRange = 0:70;
colorsDubuisTimes = hsv(length(DubuisTimeRange)); % Colormap "jet" is another option
FractonalDubuisTimeRange = DubuisTimeRange/max(DubuisTimeRange);

y_positions = [0, 0.5, 0.3, 0.6, 0.8];
unique_temperatures = unique(AllSetInfo.Temperatures);
BinsToFit = zeros(NChannels, NumAPbins, 'logical');
BinsToFit(2,9:33) = true;
BinsToFit(3,5:13) = true;
BinsToFit(4,22:37) = true;
BinsToFit(5,5:21) = true;
i_list = [1 2 4 5 7 8  11 13 14];
j_list = [2 4 5 7 8 11 13 14 1];

for ch_index =[5 3 2]% 2:5
    for i_index = 1:10
        i = i_list(i_index);
        j_index = j_list(i_index);
        for j = j_index
            outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'control','_T',...
                strrep(num2str(AllSetInfo.Temperatures(j)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(j)), 'control_',ChannelNames{ch_index}];
            
            
            RefBin = find(round(APbins,6) == y_positions(ch_index));
            close all
            DeltaFCFixCorrectedFig = figure(2);
            set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
            set(gcf,'color','w');
            DeltaFCAx = axes(DeltaFCFixCorrectedFig);
            
            map = colormap(colorsDubuisTimes);
            h = colorbar;
            % %set(h, 'ylim', [min(Prefix_temp_obs) max(Prefix_temp_obs)])
            hold off
            colorTitleHandle = get(h,'Title');
            titleString = '\del (\mum)';
            titleString = 'time (min)';
            set(colorTitleHandle ,'String',titleString);
            h.Ticks =  FractonalDubuisTimeRange(1:5:71)+0.5*(FractonalDubuisTimeRange(2)-FractonalDubuisTimeRange(1)); %Create 8 ticks from zero to 1
            h.TickLabels = {'0.0', '5.0','10.0','15.0','20.0','25.0', '30.0', '35.0', '40.0','45.0','50.0','55.0', '60.0', '65.0', '70.0'};
            hold on
            
            
            AllDubuisTimes = AllCompiledEmbryos{i}.DubuisSmoothedTimes;
            
            UseTimes1 = AllSetInfo.ControlDubuisTimesToUse{i}(ch_index,:);
            AllowedDubuisTimes1 = AllDubuisTimes(UseTimes1);
            AllAPProfileExp1 = AllCompiledEmbryos{i}.DubuisTimeSmoothedAvgAPProfiles.Control(:,BinsToFit(ch_index,:),ch_index);
            APbinsMat1 = repmat(APbins, 1, size(AllAPProfileExp1,1));
            DubuisTimesMat1 = repmat(AllCompiledEmbryos{i}.DubuisSmoothedTimes.', 1, sum(BinsToFit(ch_index,:)));
            APbinsFlat1 = APbinsMat1(:).';
            DubuisTimesFlat1 = DubuisTimesMat1(:).';
            AllAPProfileFlatExp1 =  AllAPProfileExp1(:).';
            
            UseTimes2 = AllSetInfo.ControlDubuisTimesToUse{j}(ch_index,:);
            AllowedDubuisTimes2 = AllDubuisTimes(UseTimes2);
            AllAPProfileExp2 = AllCompiledEmbryos{j}.DubuisTimeSmoothedAvgAPProfiles.Control(:,BinsToFit(ch_index,:),ch_index);
            APbinsMat2 = repmat(APbins, 1, size(AllAPProfileExp2,1));
            DubuisTimesMat2 = repmat(AllCompiledEmbryos{j}.DubuisSmoothedTimes.', 1, sum(BinsToFit(ch_index,:)));
            APbinsFlat2 = APbinsMat2(:).';
            DubuisTimesFlat2 = DubuisTimesMat2(:).';
            AllAPProfileFlatExp2 =  AllAPProfileExp2(:).';
            
            
            TFNotNaN = ~isnan(AllAPProfileFlatExp1) & ~isnan(AllAPProfileFlatExp2) &...
                ismember(DubuisTimesFlat1,AllowedDubuisTimes1) & ismember(DubuisTimesFlat2, AllowedDubuisTimes2) ;
            AllAPProfileFlatExp1 = AllAPProfileFlatExp1(TFNotNaN);
            AllAPProfileFlatExp2 = AllAPProfileFlatExp2(TFNotNaN);
            PlotDubuisTimes = DubuisTimesFlat2(TFNotNaN);
            
            for plot_index = 1:length(AllAPProfileFlatExp1)
                scatter(AllAPProfileFlatExp1(plot_index),AllAPProfileFlatExp2(plot_index),...
                    75, 'MarkerFaceColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:),...
                    'MarkerEdgeColor', colorsDubuisTimes(PlotDubuisTimes(plot_index)+1,:));
                hold on
            end
            
            
            hold on
            
            
            ymax = ceil(max([AllAPProfileFlatExp1 ; AllAPProfileFlatExp2])/.5)*.5;
            xmax = ceil(max([AllAPProfileFlatExp1 ; AllAPProfileFlatExp2])/.5)*.5;
            ymax = ceil(max(AllAPProfileFlatExp2)/5)*5;
            xmax = ceil(max(AllAPProfileFlatExp1)/5)*5;
            plot([0, ymax], [0, ymax], 'k');
            
            
            grid on
            hold off
            xlab = ['Control for ', num2str(AllSetInfo.Temperatures(i)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(i))];
            ylab = ['Control for ', num2str(AllSetInfo.Temperatures(j)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(j))];
            %xlab = [num2str(AllSetInfo.Temperatures(i)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(i))];
            %ylab = [num2str(AllSetInfo.Temperatures(j)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(j))];
            
            if ~isempty(ymax)
                ylim([0, ymax])
            end
            if ~isempty(xmax)
                xlim([0, xmax])
            end
            xlabel(xlab, 'FontSize', 16)
            ylabel(ylab, 'FontSize', 16)
            DeltaFCAx.YAxis.FontSize = 16;
            DeltaFCAx.XAxis.FontSize = 16;
            
            PlotTitle = [num2str(AllSetInfo.Temperatures(i)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(i)), ' ',ChannelNames{ch_index}];
            if AllSetInfo.Replicates(i) == 0
                PlotTitle =  [num2str(AllSetInfo.Temperatures(i)), 'ºC Flipped Stain ', ' ',ChannelNames{ch_index}];
            end
            title(DeltaFCAx, PlotTitle, 'FontSize', 18)
            %
            h.FontSize = 16;
            
            DeltaFCAx.FontSize = 16;
            %set(DeltaFCAx, 'YScale', 'log')
            %set(DeltaFCAx, 'XScale', 'log')
            outpath = [AllSetsProfFigPath, filesep,'SetComps' filesep outstring,'_DeltaFCcolor.png'];
            saveas(DeltaFCFixCorrectedFig,outpath);
            
        end
    end
end

%%
NControls = cell(1, 15);
NTests = cell(1, 15);
for i = 1:15
    NControls{i} = zeros(1, length(AllSetInfo.Prefixes{i}));
    NTests{i} = zeros(1, length(AllSetInfo.Prefixes{i}));
    for j = 1:length(AllSetInfo.Prefixes{i})
        Prefix = AllSetInfo.Prefixes{i}{j};
        CE = AllCompiledEmbryos{i};
        NControls{i}(j) = sum(CE.IsNC14 & CE.ControlSetEmbryos & CE.SlideIDs == j);
        NTests{i}(j) = sum(CE.IsNC14 & CE.TestSetEmbryos & CE.SlideIDs == j);
        disp(['i = ', num2str(i), ', j = ', num2str(j), ', NTests = ', num2str(NTests{i}(j)), ', NControls = ', num2str(NControls{i}(j))]);
    end
end

%%
for i = 1:15
    RefBin = 13;
    CompiledEmbryos = AllCompiledEmbryos{i};
    close all
    DeltaFCFixCorrectedFig = figure(2);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    
    %
    %             AllDubuisTimes = AllCompiledEmbryos{i}.DubuisSmoothedTimes;
    ch2vals = CompiledEmbryos.DorsalAvgAPProfiles(:,RefBin,2).';
    ch3vals = CompiledEmbryos.DorsalAvgAPProfiles(:,RefBin,3).';
    SlideIDs = CompiledEmbryos.SlideIDs;
    IsControl = CompiledEmbryos.ControlSetEmbryos;
    IsTest = CompiledEmbryos.TestSetEmbryos;
    DubuisTimes = CompiledEmbryos.DubuisEmbryoTimes;
    
    TFNotNaN = ~isnan(ch2vals) & ~isnan(ch3vals);
    ch2vals = ch2vals(TFNotNaN);
    ch3vals = ch3vals(TFNotNaN);
    SlideIDs = SlideIDs(TFNotNaN);
    DubuisTimes = DubuisTimes(TFNotNaN);
    IsControl = IsControl(TFNotNaN);
    IsTest = IsTest(TFNotNaN);
    
    for slide_index = 1:length(AllSetInfo.Prefixes{i})
        scatter(ch2vals(SlideIDs == slide_index & IsTest),ch3vals(SlideIDs == slide_index & IsTest),...
            75, 'o', 'MarkerFaceColor', gap_colors(slide_index+1,:),...
            'MarkerEdgeColor', gap_colors(slide_index+1,:));
        hold on
        
        scatter(ch2vals(SlideIDs == slide_index & IsControl),ch3vals(SlideIDs == slide_index & IsControl),...
            75, 's', 'MarkerFaceColor', gap_colors(slide_index+1,:),...
            'MarkerEdgeColor','k');
    end
    
    PlotTitle = [num2str(AllSetInfo.Temperatures(i)), 'ºC Replicate ', num2str(AllSetInfo.Replicates(i))];%, ' ',ChannelNames{ch_index}];
    if AllSetInfo.Replicates(i) == 0
        PlotTitle =  [num2str(AllSetInfo.Temperatures(i)), 'ºC Flipped Stain '];%, ' ',ChannelNames{ch_index}];
    end
    title(DeltaFCAx, PlotTitle, 'FontSize', 18)
    xlabel('Hoechst Fluo', 'FontSize', 16)
            ylabel('Bicoid Fluo', 'FontSize', 16)
    %

    outstring = ['T', strrep(num2str(AllSetInfo.Temperatures(i)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(i)), 'test'];%,'_T',...
                %strrep(num2str(AllSetInfo.Temperatures(j)), '.', '_'),'CRep', num2str(AllSetInfo.Replicates(j)), 'test_',ChannelNames{ch_index}];
            
            
    DeltaFCAx.FontSize = 16;
    %set(DeltaFCAx, 'YScale', 'log')
    %set(DeltaFCAx, 'XScale', 'log')
    outpath = [AllSetsProfFigPath, filesep,'SetComps' filesep outstring,'_HoechstHunchbackFluoComps.png'];
    saveas(DeltaFCFixCorrectedFig,outpath);
    
end


%% First Fit 25ºC  Replicate 2 Test and  25ºC  Flipped Test
i = 2;
j = 3;
Set1Test = true;
Set2Test = true;

AllDubuisTimes = AllCompiledEmbryos{i}.DubuisSmoothedTimes;

if Set1Test
    UseTimes1 = AllSetInfo.TestDubuisTimesToUse{i}(ch_index,:);
    AllAPProfileExp1 = AllCompiledEmbryos{i}.DubuisTimeSmoothedAvgAPProfiles.Test(:,BinsToFit(ch_index,:),ch_index);
else
    UseTimes1 = AllSetInfo.ControlDubuisTimesToUse{i}(ch_index,:);
    AllAPProfileExp1 = AllCompiledEmbryos{i}.DubuisTimeSmoothedAvgAPProfiles.Control(:,BinsToFit(ch_index,:),ch_index);
end

AllowedDubuisTimes1 = AllDubuisTimes(UseTimes1);
APbinsMat1 = repmat(APbins, 1, size(AllAPProfileExp1,1));
DubuisTimesMat1 = repmat(AllCompiledEmbryos{i}.DubuisSmoothedTimes.', 1, sum(BinsToFit(ch_index,:)));
APbinsFlat1 = APbinsMat1(:).';
DubuisTimesFlat1 = DubuisTimesMat1(:).';
FitSet1 =  AllAPProfileExp1(:).';

if Set1Test
    UseTimes2 = AllSetInfo.TestDubuisTimesToUse{j}(ch_index,:);
    AllAPProfileExp2 = AllCompiledEmbryos{j}.DubuisTimeSmoothedAvgAPProfiles.Test(:,BinsToFit(ch_index,:),ch_index);
else
    UseTimes2 = AllSetInfo.ControlDubuisTimesToUse{j}(ch_index,:);
    AllAPProfileExp2 = AllCompiledEmbryos{j}.DubuisTimeSmoothedAvgAPProfiles.Control(:,BinsToFit(ch_index,:),ch_index);
end

AllowedDubuisTimes2 = AllDubuisTimes(UseTimes2);
APbinsMat2 = repmat(APbins, 1, size(AllAPProfileExp2,1));
DubuisTimesMat2 = repmat(AllCompiledEmbryos{j}.DubuisSmoothedTimes.', 1, sum(BinsToFit(ch_index,:)));
APbinsFlat2 = APbinsMat2(:).';
DubuisTimesFlat2 = DubuisTimesMat2(:).';
FitSet2 =  AllAPProfileExp2(:).';


dlm = fitlm(FitSet2, FitSet1);%,'Weights',  1./(FitSD1.^2));
%                 OkFitInfo{exp2, exp1, ch_index} = dlm;
%                 OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
%                 OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
%                 OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
%                 OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);

%% DUBUIS Smoothed PROGRAM
% The Fit 25ºC  Flipped  Test and 25ºC  Flipped Control
% The Fit 25ºC  Flipped  Test and 27.5ºC  Replicate 2 Control
% The Fit 25ºC  Flipped  Test and 27.5ºC  Flipped Control
% The Fit 25ºC  Flipped  Test and 22.5ºC  Flipped Control
% The Fit 25ºC  Flipped  Test and 20ºC  Flipped Control
% The Fit 25ºC  Flipped  Test and 17.5ºC  Replicate 2 Control
% The Fit 25ºC  Flipped  Test and 25ºC  Replicate 2 Control
% The Fit 25ºC  Flipped  Test and 22.5ºC  Replicate 2 Control


% Then fit 27.5ºC  Replicate 1 Test to  27.5ºC  Replicate 2 Test
% fit 27.5ºC  Replicate 1 Test to  27.5ºC  Flipped Test, and  27.5ºC  Replicate 2 Test to  27.5ºC  Flipped Test
% fit 22.5ºC  Replicate 2 Test to  22.5ºC  Flipped Test
% fit 20ºC  Replicate 2 Test to  20ºC  Flipped Test
% fit 17.5ºC  Replicate 2 Test to  17.5ºC  Flipped Test
% fit 17.5ºC  Replicate 1 Test to  17.5ºC  Flipped Test

%% DUBUIS Binned Program 
% 25ºC Flipped Control  27.5ºC Flipped Control & 22.5ºC Flipped Control &
% 20ºC Flipped Control & 17.5ºC Flipped Control ( can also include 25ºC Rep
% 1 Test & 25C Rep 2 Test)
FitBins = zeros(NChannels, NumAPbins, 'logical');
FitBins(2,:) = APbins >= 0.1 &  APbins <= 0.9 ;
FitBins(3,:) = APbins >= 0.1 &  APbins <= 0.9 ;
FitBins(4,:) = APbins >= 0.1 &  APbins <= 0.9 ;
FitBins(5,:) = (APbins >= 0.1 &  APbins <= 0.9) | (APbins >= 0.7 * APbins <= 0.9) ;
Ts = [25, 27.5, 22.5, 20, 17.5, 25, 25];
Reps = [0, 0, 0, 0, 0, 1, 2];
CTstrings = {'Control', 'Control', 'Control', 'Control', 'Control', 'Test', 'Test'};
NSubsets = length(Ts);
DubuisTimes = AllCompiledEmbryos{1}.WindowedProfiles.DubuisTime.x;
Ndt= length(DubuisTimes);
means = NaN(Ndt, NumAPbins, NChannels, NSubsets);
ses = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts_below = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts_above = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts_balance =  NaN(Ndt, NumAPbins, NChannels, NSubsets);
for i = 1:NSubsets
    set_index = find(AllSetInfo.Temperatures == Ts(i) & AllSetInfo.Replicates == Reps(i));
    means(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).mean;
    ses(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).se;
    counts(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count;
    counts_below(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count_belowcenter;
    counts_above(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count_abovecenter;
    counts_balance(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count_balance;
end


StrictConditions = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
LooseConditions = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
StrictConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets,NSubsets, 'logical');
LooseConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets,NSubsets, 'logical');
StrictCompactComps = zeros(NChannels, NSubsets,NSubsets);
LooseCompactComps = zeros(NChannels, NSubsets,NSubsets);

SubsetsIncluded = cell(1, NChannels);
Fits = cell(1, NChannels);

Slopes = NaN(NSubsets-1, 2, NChannels);
Intercepts = NaN(NSubsets-1, 2, NChannels);
IncludedMeans = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
IncludedSEs = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
IncludedCounts = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
CombinedMean = NaN(length(DubuisTimes), NumAPbins, NChannels);
CombinedSE = NaN(length(DubuisTimes), NumAPbins, NChannels);
CombinedCounts = NaN(length(DubuisTimes), NumAPbins, NChannels);

%%

for ch_index = 2:5


for i = 1:NSubsets
    for apbin = 1:NumAPbins
            counts_i = counts(:,apbin,ch_index,i);
            counts_above_i = counts_above(:,apbin,ch_index,i);
            counts_below_i = counts_below(:,apbin,ch_index,i);
            counts_balance_i = counts_balance(:,apbin,ch_index,i);
            StrictConditions(:,apbin,ch_index,i) = (counts_i >= 5 & counts_below_i  >= 2 ...
                & counts_above_i >= 2 & counts_balance_i>= 0.3 & counts_balance_i <= 0.7);
            LooseConditions(:,apbin,ch_index,i) = (counts_i >= 3 & counts_below_i  >= 1 ...
                & counts_above_i >= 1 & counts_balance_i >= 0.2 & counts_balance_i <= 0.8);
    end
    
end


AllSubsets = 1:NSubsets;
for i = 1:NSubsets
for j = AllSubsets(~ismember(AllSubsets, i))
        for apbin = 1:NumAPbins
   
            StrictConditionsComps(:,apbin,ch_index,i,j) = StrictConditions(:,apbin,ch_index,i) &...
                StrictConditions(:,apbin,ch_index,j);
            LooseConditionsComps(:,apbin,ch_index,i,j) = LooseConditions(:,apbin,ch_index,i) &...
                LooseConditions(:,apbin,ch_index,j);
        end
        SumStrictTimeBins = sum(StrictConditionsComps(:,:,ch_index,i,j), 1);
        StrictCompactComps(ch_index, i, j) = max(SumStrictTimeBins);
        SumLooseTimeBins = sum(LooseConditionsComps(:,:,ch_index,i,j), 1);
        LooseCompactComps(ch_index, i, j) = max(SumLooseTimeBins);
 
    
end
end

DubuisMats = repmat(DubuisTimes.', 1, NumAPbins);
APMats = repmat(APbins, length(DubuisTimes), 1);
RowLabels = repmat((1:NSubsets).', 1, NSubsets);
ColLabels = repmat(1:NSubsets, NSubsets, 1);
SubRowLabels = repmat((1:(NSubsets-2)).', 1, NSubsets-2);
SubColLabels = repmat(1:(NSubsets-2), NSubsets-2, 1);
SubsetsToFit = 1:NSubsets;
SubsetsIncluded{ch_index} = [];
Fits{ch_index} = {};

[maxval, max_idx] = max(squeeze(LooseCompactComps(ch_index, 1:5, 1:5)),[], 'all', 'Linear');
set1 = SubRowLabels(max_idx);
set2 = SubColLabels(max_idx);
SubsetsIncluded{ch_index} = [SubsetsIncluded{ch_index} set1 set2];
SubsetsToFit = SubsetsToFit(~ismember(SubsetsToFit, SubsetsIncluded{ch_index}));


%%

means1 = means(:,FitBins(ch_index,:),ch_index, set1);
means1 = means1(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
dubuis1 = DubuisMats(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
counts1 =  counts(:,FitBins(ch_index,:),ch_index, set1);
counts1 = counts1(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
ses1 = ses(:,FitBins(ch_index,:),ch_index, set1);
ses1 = ses1(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
means2 = means(:,FitBins(ch_index,:),ch_index, set2);
means2 = means2(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
dubuis2 = DubuisMats(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
ses2 = ses(:,FitBins(ch_index,:),ch_index, set2);
ses2 = ses2(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
counts2 =  counts(:,FitBins(ch_index,:),ch_index, set2);
counts2 = counts2(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
IncludedMeans(:,:,ch_index,1) =  means(:,:,ch_index, set1);
IncludedSEs(:,:,ch_index,1) = ses(:,:,ch_index, set1);
IncludedCounts(:,:,ch_index,1) = counts(:,:,ch_index, set1);

dlm = fitlm(means1, means2);
Fits{ch_index}{1} = dlm;
Slopes(size(Slopes, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(2) dlm.Coefficients.SE(2) ];
Intercepts(size(Intercepts, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(1) dlm.Coefficients.SE(1) ];
IncludedMeans(:,:,ch_index,2) =  means(:,:,ch_index, set2)*dlm.Coefficients.Estimate(2)+dlm.Coefficients.Estimate(1);
IncludedSEs(:,:,ch_index,2) =  ses(:,:,ch_index, set2)*dlm.Coefficients.Estimate(2);
IncludedCounts(:,:,ch_index,2) = counts(:,:,ch_index, set2);

for t_index = 1:length(DubuisTimes)
    for ap_index = 1:NumAPbins
        MeansToUse = ~isnan(squeeze(IncludedMeans(t_index, ap_index,ch_index, :))).';
        if sum(MeansToUse) > 0
            Weights = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).'/(sum(squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse))));
            CountsToUse = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse = squeeze(IncludedSEs(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse(isnan(SEsToUse)) = 0;
            CombinedMean(t_index,ap_index, ch_index)  = sum((squeeze(IncludedMeans(t_index, ap_index,ch_index, MeansToUse)).').*Weights);
            CombinedSE(t_index,ap_index, ch_index)  = sqrt(sum(SEsToUse.*SEsToUse.*Weights.*Weights));
            CombinedCounts(t_index,ap_index, ch_index)   = sum(CountsToUse);
        else
            CombinedMean(t_index,ap_index, ch_index) = NaN;
            CombinedSE(t_index,ap_index, ch_index) = NaN;
            CombinedCounts(t_index,ap_index, ch_index) = 0;
        end
        
    end
end

CombinedMean(CombinedMean == 0) = NaN;
plot(APbins, CombinedMean(:,:,ch_index).')
disp('pause')


%%

while ~isempty(SubsetsToFit) 
TempStrictConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
TempLooseConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
TempStrictCompactComps = zeros(NChannels,NSubsets);
TempLooseCompactComps = zeros(NChannels, NSubsets);

for i = SubsetsToFit
    
    for apbin = 1:NumAPbins
        
        TempStrictConditionsComps(:,apbin,ch_index,i) = StrictConditions(:,apbin,ch_index,SubsetsIncluded{ch_index}(1)); 
        TempLooseConditionsComps(:,apbin,ch_index,i) = LooseConditions(:,apbin,ch_index,SubsetsIncluded{ch_index}(1)); 
        for j = SubsetsIncluded{ch_index}(2:end)
            TempStrictConditionsComps(:,apbin,ch_index,i) = TempStrictConditionsComps(:,apbin,ch_index,i) |  StrictConditions(:,apbin,ch_index,j);
            TempLooseConditionsComps(:,apbin,ch_index,i) = TempLooseConditionsComps(:,apbin,ch_index,i) |  LooseConditions(:,apbin,ch_index,j);
        end
         TempStrictConditionsComps(:,apbin,ch_index,i) =  TempStrictConditionsComps(:,apbin,ch_index,i) &  StrictConditions(:,apbin,ch_index,i);
         TempLooseConditionsComps(:,apbin,ch_index,i) =  TempLooseConditionsComps(:,apbin,ch_index,i) &  LooseConditions(:,apbin,ch_index,i);
        
    end
    SumStrictTimeBins = sum(TempStrictConditionsComps(:,:,ch_index,i), 1);
    TempStrictCompactComps(ch_index, i) = max(SumStrictTimeBins);
    SumLooseTimeBins = sum(TempLooseConditionsComps(:,:,ch_index,i), 1);
    TempLooseCompactComps(ch_index, i) = max(SumLooseTimeBins);
    
    
end

if length(SubsetsIncluded{ch_index}) < NSubsets - 2
    [maxval, max_idx] = max(squeeze(TempLooseCompactComps(ch_index, 1:(NSubsets-2))),[], 'all', 'Linear');
    tempset = SubRowLabels(max_idx);
else
    [maxval, max_idx] = max(squeeze(TempLooseCompactComps(ch_index, 1:(NSubsets))),[], 'all', 'Linear');
    tempset = RowLabels(max_idx);
end
SubsetsIncluded{ch_index} = [SubsetsIncluded{ch_index}  tempset];
SubsetsToFit = SubsetsToFit(~ismember(SubsetsToFit, SubsetsIncluded{ch_index} ));


%%

means1 = means(:,FitBins(ch_index,:),ch_index, tempset);
means1 = means1(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
dubuis1 = DubuisMats(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
counts1 =  counts(:,FitBins(ch_index,:),ch_index, tempset);
counts1 = counts1(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
ses1 = ses(:,FitBins(ch_index,:),ch_index, tempset);
ses1 = ses1(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 



means2 = CombinedMean(:,FitBins(ch_index,:),ch_index);
means2 = means2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
dubuis2 = DubuisMats(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
ses2 = CombinedSE(:,FitBins(ch_index,:),ch_index);
ses2 = ses2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 

counts2 =  CombinedCounts(:,FitBins(ch_index,:),ch_index);
counts2 = counts2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 



dlm = fitlm(means1, means2)
Fits{ch_index}{length(SubsetsIncluded{ch_index})-1} = dlm;
Slopes(size(Slopes, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(2) dlm.Coefficients.SE(2) ];
Intercepts(size(Intercepts, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(1) dlm.Coefficients.SE(1) ];
IncludedMeans(:,:,ch_index,length(SubsetsIncluded{ch_index}))=  means(:,:,ch_index, tempset)*dlm.Coefficients.Estimate(2)+dlm.Coefficients.Estimate(1);
IncludedSEs(:,:,ch_index,length(SubsetsIncluded{ch_index})) =  ses(:,:,ch_index, tempset)*dlm.Coefficients.Estimate(2);
IncludedCounts(:,:,ch_index,length(SubsetsIncluded{ch_index})) = counts(:,:,ch_index, tempset);

for t_index = 1:length(DubuisTimes)
    for ap_index = 1:NumAPbins
        MeansToUse = ~isnan(squeeze(IncludedMeans(t_index, ap_index,ch_index, :))).';
        if sum(MeansToUse) > -1
            Weights = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).'/(sum(squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse))));
            CountsToUse = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse = squeeze(IncludedSEs(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse(isnan(SEsToUse)) = 0;
            CombinedMean(t_index,ap_index, ch_index)  = sum((squeeze(IncludedMeans(t_index, ap_index,ch_index, MeansToUse)).').*Weights);
            CombinedSE(t_index,ap_index, ch_index)  = sqrt(sum(SEsToUse.*SEsToUse.*Weights.*Weights));
            CombinedCounts(t_index,ap_index, ch_index)   = sum(CountsToUse);
        else
            CombinedMean(t_index,ap_index, ch_index) = NaN;
            CombinedSE(t_index,ap_index, ch_index) = NaN;
            CombinedCounts(t_index,ap_index, ch_index) = 0;
        end
    end
end

CombinedMean(CombinedMean == 0) = NaN;
plot(APbins, CombinedMean(:,:,ch_index).')

disp('pause')
end
end
%%
AllCompiledPath = 'S:/Gabriella/Dropbox/ProteinProfiles/25CMasterSets.mat';
save(AllCompiledPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');



%% DUBUIS Smoothed Profile 
% 25ºC Flipped Control  27.5ºC Flipped Control & 22.5ºC Flipped Control &
% 20ºC Flipped Control & 17.5ºC Flipped Control ( can also include 25ºC Rep
% 1 Test & 25C Rep 2 Test)
Ts = [25, 25,25, 27.5,27.5, 22.5,22.5,20, 17.5, 17.5];
Reps = [0, 1, 2, 1, 2, 1, 2, 2, 1, 2];
CTstrings = {'Test', 'Control', 'Control', 'Control', 'Control', 'Control',  'Control', 'Control', 'Control','Control'};
NSubsets = length(Ts);
DubuisTimes = AllCompiledEmbryos{1}.WindowedProfiles.DubuisTime.x;
Ndt= length(DubuisTimes);
means = NaN(Ndt, NumAPbins, NChannels, NSubsets);
ses = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts_below = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts_above = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts_balance =  NaN(Ndt, NumAPbins, NChannels, NSubsets);
for i = 1:NSubsets
    set_index = find(AllSetInfo.Temperatures == Ts(i) & AllSetInfo.Replicates == Reps(i));
    means(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).mean;
    ses(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).se;
    counts(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count;
    counts_below(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count_belowcenter;
    counts_above(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count_abovecenter;
    counts_balance(:,:,:,i) = AllCompiledEmbryos{set_index}.WindowedProfiles.DubuisTime.(CTstrings{i}).count_balance;
end


StrictConditions = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
LooseConditions = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
StrictConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets,NSubsets, 'logical');
LooseConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets,NSubsets, 'logical');
StrictCompactComps = zeros(NChannels, NSubsets,NSubsets);
LooseCompactComps = zeros(NChannels, NSubsets,NSubsets);

SubsetsIncluded = cell(1, NChannels);
Fits = cell(1, NChannels);

Slopes = NaN(NSubsets-1, 2, NChannels);
Intercepts = NaN(NSubsets-1, 2, NChannels);
IncludedMeans = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
IncludedSEs = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
IncludedCounts = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
CombinedMean = NaN(length(DubuisTimes), NumAPbins, NChannels);
CombinedSE = NaN(length(DubuisTimes), NumAPbins, NChannels);
CombinedCounts = NaN(length(DubuisTimes), NumAPbins, NChannels);

%%

for ch_index = 2:5

for i = 1:NSubsets
    for apbin = 1:NumAPbins
            counts_i = counts(:,apbin,ch_index,i);
            counts_above_i = counts_above(:,apbin,ch_index,i);
            counts_below_i = counts_below(:,apbin,ch_index,i);
            counts_balance_i = counts_balance(:,apbin,ch_index,i);
            StrictConditions(:,apbin,ch_index,i) = (counts_i >= 5 & counts_below_i  >= 2 ...
                & counts_above_i >= 2 & counts_balance_i>= 0.3 & counts_balance_i <= 0.7);
            LooseConditions(:,apbin,ch_index,i) = (counts_i >= 2 & counts_below_i  >= 1 ...
                & counts_above_i >= 1 & counts_balance_i >= 0.2 & counts_balance_i <= 0.8);
    end
    
end


AllSubsets = 1:NSubsets;
for i = 1:NSubsets
for j = AllSubsets(~ismember(AllSubsets, i))
        for apbin = 1:NumAPbins
   
            StrictConditionsComps(:,apbin,ch_index,i,j) = StrictConditions(:,apbin,ch_index,i) &...
                StrictConditions(:,apbin,ch_index,j);
            LooseConditionsComps(:,apbin,ch_index,i,j) = LooseConditions(:,apbin,ch_index,i) &...
                LooseConditions(:,apbin,ch_index,j);
        end
        SumStrictTimeBins = sum(StrictConditionsComps(:,:,ch_index,i,j), 1);
        StrictCompactComps(ch_index, i, j) = max(SumStrictTimeBins);
        SumLooseTimeBins = sum(LooseConditionsComps(:,:,ch_index,i,j), 1);
        LooseCompactComps(ch_index, i, j) = max(SumLooseTimeBins);
 
    
end
end

DubuisMats = repmat(DubuisTimes.', 1, NumAPbins);

RowLabels = repmat((1:NSubsets).', 1, NSubsets);
ColLabels = repmat(1:NSubsets, NSubsets, 1);
SubRowLabels = repmat((1:(NSubsets-2)).', 1, NSubsets-2);
SubColLabels = repmat(1:(NSubsets-2), NSubsets-2, 1);
SubsetsToFit = 1:NSubsets;
SubsetsIncluded{ch_index} = [];
Fits{ch_index} = {};

[maxval, max_idx] = max(squeeze(LooseCompactComps(ch_index, 1:(NSubsets), 1:(NSubsets))),[], 'all', 'Linear');
set1 = RowLabels(max_idx);
set2 = ColLabels(max_idx);
SubsetsIncluded{ch_index} = [SubsetsIncluded{ch_index} set1 set2];
SubsetsToFit = SubsetsToFit(~ismember(SubsetsToFit, SubsetsIncluded{ch_index}));


%%

means1 = means(:,FitBins(ch_index,:),ch_index, set1);
means1 = means1(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
dubuis1 = DubuisMats(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
counts1 =  counts(:,FitBins(ch_index,:),ch_index, set1);
counts1 = counts1(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
ses1 = ses(:,FitBins(ch_index,:),ch_index, set1);
ses1 = ses1(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
means2 = means(:,FitBins(ch_index,:),ch_index, set2);
means2 = means2(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
dubuis2 = DubuisMats(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
ses2 = ses(:,FitBins(ch_index,:),ch_index, set2);
ses2 = ses2(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
counts2 =  counts(:,FitBins(ch_index,:),ch_index, set2);
counts2 = counts2(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
IncludedMeans(:,:,ch_index,1) =  means(:,:,ch_index, set1);
IncludedSEs(:,:,ch_index,1) = ses(:,:,ch_index, set1);
IncludedCounts(:,:,ch_index,1) = counts(:,:,ch_index, set1);

dlm = fitlm(means1, means2);
Fits{ch_index}{1} = dlm;
Slopes(size(Slopes, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(2) dlm.Coefficients.SE(2) ];
Intercepts(size(Intercepts, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(1) dlm.Coefficients.SE(1) ];
IncludedMeans(:,:,ch_index,2) =  means(:,:,ch_index, set2)*dlm.Coefficients.Estimate(2)+dlm.Coefficients.Estimate(1);
IncludedSEs(:,:,ch_index,2) =  ses(:,:,ch_index, set2)*dlm.Coefficients.Estimate(2);
IncludedCounts(:,:,ch_index,2) = counts(:,:,ch_index, set2);

for t_index = 1:length(DubuisTimes)
    for ap_index = 1:NumAPbins
        MeansToUse = ~isnan(squeeze(IncludedMeans(t_index, ap_index,ch_index, :))).';
        if sum(MeansToUse) > 0
            Weights = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).'/(sum(squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse))));
            CountsToUse = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse = squeeze(IncludedSEs(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse(isnan(SEsToUse)) = 0;
            CombinedMean(t_index,ap_index, ch_index)  = sum((squeeze(IncludedMeans(t_index, ap_index,ch_index, MeansToUse)).').*Weights);
            CombinedSE(t_index,ap_index, ch_index)  = sqrt(sum(SEsToUse.*SEsToUse.*Weights.*Weights));
            CombinedCounts(t_index,ap_index, ch_index)   = sum(CountsToUse);
        else
            CombinedMean(t_index,ap_index, ch_index) = NaN;
            CombinedSE(t_index,ap_index, ch_index) = NaN;
            CombinedCounts(t_index,ap_index, ch_index) = 0;
        end
        
    end
end

CombinedMean(CombinedMean == 0) = NaN;
plot(APbins, CombinedMean(:,:,ch_index).')
disp('pause')


%%

while ~isempty(SubsetsToFit) 
TempStrictConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
TempLooseConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
TempStrictCompactComps = zeros(NChannels,NSubsets);
TempLooseCompactComps = zeros(NChannels, NSubsets);

for i = SubsetsToFit
    
    for apbin = 1:NumAPbins
        
        TempStrictConditionsComps(:,apbin,ch_index,i) = StrictConditions(:,apbin,ch_index,SubsetsIncluded{ch_index}(1)); 
        TempLooseConditionsComps(:,apbin,ch_index,i) = LooseConditions(:,apbin,ch_index,SubsetsIncluded{ch_index}(1)); 
        for j = SubsetsIncluded{ch_index}(2:end)
            TempStrictConditionsComps(:,apbin,ch_index,i) = TempStrictConditionsComps(:,apbin,ch_index,i) |  StrictConditions(:,apbin,ch_index,j);
            TempLooseConditionsComps(:,apbin,ch_index,i) = TempLooseConditionsComps(:,apbin,ch_index,i) |  LooseConditions(:,apbin,ch_index,j);
        end
         TempStrictConditionsComps(:,apbin,ch_index,i) =  TempStrictConditionsComps(:,apbin,ch_index,i) &  StrictConditions(:,apbin,ch_index,i);
         TempLooseConditionsComps(:,apbin,ch_index,i) =  TempLooseConditionsComps(:,apbin,ch_index,i) &  LooseConditions(:,apbin,ch_index,i);
        
    end
    SumStrictTimeBins = sum(TempStrictConditionsComps(:,:,ch_index,i), 1);
    TempStrictCompactComps(ch_index, i) = max(SumStrictTimeBins);
    SumLooseTimeBins = sum(TempLooseConditionsComps(:,:,ch_index,i), 1);
    TempLooseCompactComps(ch_index, i) = max(SumLooseTimeBins);
    
    
end


[maxval, max_idx] = max(squeeze(TempLooseCompactComps(ch_index, 1:(NSubsets))),[], 'all', 'Linear');
tempset = RowLabels(max_idx);

SubsetsIncluded{ch_index} = [SubsetsIncluded{ch_index}  tempset];
SubsetsToFit = SubsetsToFit(~ismember(SubsetsToFit, SubsetsIncluded{ch_index} ));


%%

means1 = means(:,FitBins(ch_index,:),ch_index, tempset);
means1 = means1(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
dubuis1 = DubuisMats(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
counts1 =  counts(:,FitBins(ch_index,:),ch_index, tempset);
counts1 = counts1(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
ses1 = ses(:,FitBins(ch_index,:),ch_index, tempset);
ses1 = ses1(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 



means2 = CombinedMean(:,FitBins(ch_index,:),ch_index);
means2 = means2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
dubuis2 = DubuisMats(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
ses2 = CombinedSE(:,FitBins(ch_index,:),ch_index);
ses2 = ses2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 

counts2 =  CombinedCounts(:,FitBins(ch_index,:),ch_index);
counts2 = counts2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 



dlm = fitlm(means1, means2)
Fits{ch_index}{length(SubsetsIncluded{ch_index})-1} = dlm;
Slopes(size(Slopes, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(2) dlm.Coefficients.SE(2) ];
Intercepts(size(Intercepts, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(1) dlm.Coefficients.SE(1) ];
IncludedMeans(:,:,ch_index,length(SubsetsIncluded{ch_index}))=  means(:,:,ch_index, tempset)*dlm.Coefficients.Estimate(2)+dlm.Coefficients.Estimate(1);
IncludedSEs(:,:,ch_index,length(SubsetsIncluded{ch_index})) =  ses(:,:,ch_index, tempset)*dlm.Coefficients.Estimate(2);
IncludedCounts(:,:,ch_index,length(SubsetsIncluded{ch_index})) = counts(:,:,ch_index, tempset);

for t_index = 1:length(DubuisTimes)
    for ap_index = 1:NumAPbins
        MeansToUse = ~isnan(squeeze(IncludedMeans(t_index, ap_index,ch_index, :))).';
        if sum(MeansToUse) > -1
            Weights = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).'/(sum(squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse))));
            CountsToUse = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse = squeeze(IncludedSEs(t_index, ap_index,ch_index, MeansToUse)).';
            SEsToUse(isnan(SEsToUse)) = 0;
            CombinedMean(t_index,ap_index, ch_index)  = sum((squeeze(IncludedMeans(t_index, ap_index,ch_index, MeansToUse)).').*Weights);
            CombinedSE(t_index,ap_index, ch_index)  = sqrt(sum(SEsToUse.*SEsToUse.*Weights.*Weights));
            CombinedCounts(t_index,ap_index, ch_index)   = sum(CountsToUse);
        else
            CombinedMean(t_index,ap_index, ch_index) = NaN;
            CombinedSE(t_index,ap_index, ch_index) = NaN;
            CombinedCounts(t_index,ap_index, ch_index) = 0;
        end
    end
end

CombinedMean(CombinedMean == 0) = NaN;
plot(APbins, CombinedMean(:,:,ch_index).')

disp('pause')
end
end
%%
AllFlippedCompiledPath = 'S:/Gabriella/Dropbox/ProteinProfiles/25CUnflippedMasterSets.mat';
save(AllFlippedCompiledPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');




%% DUBUIS Smoothed Program 
% 25ºC Flipped Control  27.5ºC Flipped Control & 22.5ºC Flipped Control &
% 20ºC Flipped Control & 17.5ºC Flipped Control ( can also include 25ºC Rep
% 1 Test & 25C Rep 2 Test)
FitBins = zeros(NChannels, NumAPbins, 'logical');
FitBins(2,:) = APbins >= 0.2 &  APbins <= 0.8 ;
FitBins(3,:) = APbins >= 0.1 &  APbins <= 0.5 ;
FitBins(4,:) = APbins >= 0.5 &  APbins <= 0.9 ;
FitBins(5,:) = (APbins >= 0.2 &  APbins <= 0.5) | (APbins >= 0.7 * APbins <= 0.9) ;
Ts = [25, 27.5, 22.5, 20, 17.5, 25, 25];
Reps = [0, 0, 0, 0, 0, 1, 2];
CTstrings = {'Control', 'Control', 'Control', 'Control', 'Control', 'Test', 'Test'};
NSubsets = length(Ts);
DubuisTimes = AllCompiledEmbryos{1}.DubuisTimesSmoothedProfiles.ControlCounts;
DubuisTimes = AllCompiledEmbryos{1}.DubuisSmoothedTimes;
Ndt= length(DubuisTimes);
means = NaN(Ndt, NumAPbins, NChannels, NSubsets);
counts = NaN(Ndt,  NSubsets);
counts_below = NaN(Ndt, NSubsets);
counts_above = NaN(Ndt, NSubsets);
counts_balance =  NaN(Ndt,  NSubsets);
for i = 1:NSubsets
    set_index = find(AllSetInfo.Temperatures == Ts(i) & AllSetInfo.Replicates == Reps(i));
    means(:,:,:,i) = AllCompiledEmbryos{set_index}.DubuisTimeSmoothedAvgAPProfiles.(CTstrings{i});
    counts(:,i) = AllCompiledEmbryos{set_index}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_1sigma+...
        AllCompiledEmbryos{set_index}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma;
    counts_below(:,i) = AllCompiledEmbryos{set_index}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma_belowcenter;
    counts_above(:,i) = AllCompiledEmbryos{set_index}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma_abovecenter;
    counts_balance(:,i) = AllCompiledEmbryos{set_index}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2balance;
end


StrictConditions = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
LooseConditions = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
StrictConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets,NSubsets, 'logical');
LooseConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets,NSubsets, 'logical');
StrictCompactComps = zeros(NChannels, NSubsets,NSubsets);
LooseCompactComps = zeros(NChannels, NSubsets,NSubsets);

SubsetsIncluded = cell(1, NChannels);
Fits = cell(1, NChannels);

Slopes = NaN(NSubsets-1, 2, NChannels);
Intercepts = NaN(NSubsets-1, 2, NChannels);
IncludedMeans = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
IncludedCounts = NaN(length(DubuisTimes), NumAPbins, NChannels, NSubsets);
CombinedMean = NaN(length(DubuisTimes), NumAPbins, NChannels);
CombinedCounts = NaN(length(DubuisTimes), NumAPbins, NChannels);

%%

for ch_index = 2:5



for i = 1:NSubsets
    for apbin = 1:NumAPbins
            counts_i = counts(:,i);
            counts_above_i = counts_above(:,i);
            counts_below_i = counts_below(:,i);
            counts_balance_i = counts_balance(:,i);
            StrictConditions(:,apbin, ch_index, i) = (counts_i >= 5 & counts_below_i  >= 2 ...
                & counts_above_i >= 2 & counts_balance_i>= 0.3 & counts_balance_i <= 0.7);
            LooseConditions(:,apbin, ch_index, i)  = (counts_i >= 3 & counts_below_i  >= 1 ...
                & counts_above_i >= 1 & counts_balance_i >= 0.2 & counts_balance_i <= 0.8);
    end
    
end


AllSubsets = 1:NSubsets;
for i = 1:NSubsets
for j = AllSubsets(~ismember(AllSubsets, i))
        for apbin = 1:NumAPbins
   
            StrictConditionsComps(:,apbin,ch_index,i,j) = StrictConditions(:,apbin,ch_index,i) &...
                StrictConditions(:,apbin,ch_index,j);
            LooseConditionsComps(:,apbin,ch_index,i,j) = LooseConditions(:,apbin,ch_index,i) &...
                LooseConditions(:,apbin,ch_index,j);
        end
        SumStrictTimeBins = sum(StrictConditionsComps(:,:,ch_index,i,j), 1);
        StrictCompactComps(ch_index, i, j) = max(SumStrictTimeBins);
        SumLooseTimeBins = sum(LooseConditionsComps(:,:,ch_index,i,j), 1);
        LooseCompactComps(ch_index, i, j) = max(SumLooseTimeBins);
 
    
end
end

DubuisMats = repmat(DubuisTimes.', 1, NumAPbins);
APMats = repmat(APbins, length(DubuisTimes), 1);
RowLabels = repmat((1:NSubsets).', 1, NSubsets);
ColLabels = repmat(1:NSubsets, NSubsets, 1);
SubRowLabels = repmat((1:(NSubsets-2)).', 1, NSubsets-2);
SubColLabels = repmat(1:(NSubsets-2), NSubsets-2, 1);
SubsetsToFit = 1:NSubsets;
SubsetsIncluded{ch_index} = [];
Fits{ch_index} = {};

[maxval, max_idx] = max(squeeze(LooseCompactComps(ch_index, 1:5, 1:5)),[], 'all', 'Linear');
set1 = SubRowLabels(max_idx);
set2 = SubColLabels(max_idx);
SubsetsIncluded{ch_index} = [SubsetsIncluded{ch_index} set1 set2];
SubsetsToFit = SubsetsToFit(~ismember(SubsetsToFit, SubsetsIncluded{ch_index}));


%%

means1 = means(:,FitBins(ch_index,:),ch_index, set1);
means1 = means1(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
dubuis1 = DubuisMats(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
counts1 =  counts(:,set1);
counts1 = counts1(LooseConditions(:,21,ch_index, set1)& LooseConditions(:,21,ch_index, set2)); 
means2 = means(:,FitBins(ch_index,:),ch_index, set2);
means2 = means2(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
dubuis2 = DubuisMats(LooseConditions(:,FitBins(ch_index,:),ch_index, set1) & LooseConditions(:,FitBins(ch_index,:),ch_index, set2)); 
counts2 =  counts(:,set2);
counts2 = counts2(LooseConditions(:,21,ch_index, set1) & LooseConditions(:,21,ch_index, set2)); 
IncludedMeans(:,:,ch_index,1) =  means(:,:,ch_index, set1);
IncludedCounts(:,:,ch_index,1) = repmat(counts(:,set1), 1, NumAPbins);

dlm = fitlm(means1, means2);
Fits{ch_index}{1} = dlm;
Slopes(size(Slopes, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(2) dlm.Coefficients.SE(2) ];
Intercepts(size(Intercepts, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(1) dlm.Coefficients.SE(1) ];
IncludedMeans(:,:,ch_index,2) =  means(:,:,ch_index, set2)*dlm.Coefficients.Estimate(2)+dlm.Coefficients.Estimate(1);
IncludedCounts(:,:,ch_index,2) = repmat(counts(:,set2), 1, NumAPbins);

for t_index = 1:length(DubuisTimes)
    for ap_index = 1:NumAPbins
        MeansToUse = ~isnan(squeeze(IncludedMeans(t_index, ap_index,ch_index, :))).';
        if sum(MeansToUse) > 0
            Weights = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).'/(sum(squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse))));
            CountsToUse = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).';
            CombinedMean(t_index,ap_index, ch_index)  = sum((squeeze(IncludedMeans(t_index, ap_index,ch_index, MeansToUse)).').*Weights);
            CombinedCounts(t_index,ap_index, ch_index)   = sum(CountsToUse);
        else
            CombinedMean(t_index,ap_index, ch_index) = NaN;
            CombinedCounts(t_index,ap_index, ch_index) = 0;
        end
        
    end
end

CombinedMean(CombinedMean == 0) = NaN;
plot(APbins, CombinedMean(:,:,ch_index).')
disp('pause')


%%

while ~isempty(SubsetsToFit) 
TempStrictConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
TempLooseConditionsComps = zeros(Ndt, NumAPbins, NChannels, NSubsets, 'logical');
TempStrictCompactComps = zeros(NChannels,NSubsets);
TempLooseCompactComps = zeros(NChannels, NSubsets);

for i = SubsetsToFit
    
    for apbin = 1:NumAPbins
        
        TempStrictConditionsComps(:,apbin,ch_index,i) = StrictConditions(:,apbin,ch_index,SubsetsIncluded{ch_index}(1)); 
        TempLooseConditionsComps(:,apbin,ch_index,i) = LooseConditions(:,apbin,ch_index,SubsetsIncluded{ch_index}(1)); 
        for j = SubsetsIncluded{ch_index}(2:end)
            TempStrictConditionsComps(:,apbin,ch_index,i) = TempStrictConditionsComps(:,apbin,ch_index,i) |  StrictConditions(:,apbin,ch_index,j);
            TempLooseConditionsComps(:,apbin,ch_index,i) = TempLooseConditionsComps(:,apbin,ch_index,i) |  LooseConditions(:,apbin,ch_index,j);
        end
         TempStrictConditionsComps(:,apbin,ch_index,i) =  TempStrictConditionsComps(:,apbin,ch_index,i) &  StrictConditions(:,apbin,ch_index,i);
         TempLooseConditionsComps(:,apbin,ch_index,i) =  TempLooseConditionsComps(:,apbin,ch_index,i) &  LooseConditions(:,apbin,ch_index,i);
        
    end
    SumStrictTimeBins = sum(TempStrictConditionsComps(:,:,ch_index,i), 1);
    TempStrictCompactComps(ch_index, i) = max(SumStrictTimeBins);
    SumLooseTimeBins = sum(TempLooseConditionsComps(:,:,ch_index,i), 1);
    TempLooseCompactComps(ch_index, i) = max(SumLooseTimeBins);
    
    
end

if length(SubsetsIncluded{ch_index}) < NSubsets - 2
    [maxval, max_idx] = max(squeeze(TempLooseCompactComps(ch_index, 1:(NSubsets-2))),[], 'all', 'Linear');
    tempset = SubRowLabels(max_idx);
else
    [maxval, max_idx] = max(squeeze(TempLooseCompactComps(ch_index, 1:(NSubsets))),[], 'all', 'Linear');
    tempset = RowLabels(max_idx);
end
SubsetsIncluded{ch_index} = [SubsetsIncluded{ch_index}  tempset];
SubsetsToFit = SubsetsToFit(~ismember(SubsetsToFit, SubsetsIncluded{ch_index} ));


%%

means1 = means(:,FitBins(ch_index,:),ch_index, tempset);
means1 = means1(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
dubuis1 = DubuisMats(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
counts1 =  counts(:,tempset);
counts1 = counts1(TempLooseConditionsComps(:,21,ch_index, tempset)); 

means2 = CombinedMean(:,FitBins(ch_index,:),ch_index);
means2 = means2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 
dubuis2 = DubuisMats(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 

counts2 =  CombinedCounts(:,FitBins(ch_index,:),ch_index);
counts2 = counts2(TempLooseConditionsComps(:,FitBins(ch_index,:),ch_index, tempset)); 



dlm = fitlm(means1, means2)
Fits{ch_index}{length(SubsetsIncluded{ch_index})-1} = dlm;
Slopes(size(Slopes, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(2) dlm.Coefficients.SE(2) ];
Intercepts(size(Intercepts, 1)+1, :, ch_index) = [dlm.Coefficients.Estimate(1) dlm.Coefficients.SE(1) ];
IncludedMeans(:,:,ch_index,length(SubsetsIncluded{ch_index}))=  means(:,:,ch_index, tempset)*dlm.Coefficients.Estimate(2)+dlm.Coefficients.Estimate(1);
IncludedCounts(:,:,ch_index,length(SubsetsIncluded{ch_index})) = repmat(counts(:,tempset), 1, NumAPbins);

for t_index = 1:length(DubuisTimes)
    for ap_index = 1:NumAPbins
        MeansToUse = ~isnan(squeeze(IncludedMeans(t_index, ap_index,ch_index, :))).';
        if sum(MeansToUse) > -1
            Weights = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).'/(sum(squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse))));
            CountsToUse = squeeze(IncludedCounts(t_index, ap_index,ch_index, MeansToUse)).';
            CombinedMean(t_index,ap_index, ch_index)  = sum((squeeze(IncludedMeans(t_index, ap_index,ch_index, MeansToUse)).').*Weights);
  
            CombinedCounts(t_index,ap_index, ch_index)   = sum(CountsToUse);
        else
            CombinedMean(t_index,ap_index, ch_index) = NaN;
            CombinedCounts(t_index,ap_index, ch_index) = 0;
        end
    end
end

CombinedMean(CombinedMean == 0) = NaN;
plot(APbins, CombinedMean(:,:,ch_index).')

disp('pause')
end
end
%%
AllSmoothCompiledPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Smoothed25CMasterSets.mat';
save(AllSmoothCompiledPath, 'CombinedMean', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');













