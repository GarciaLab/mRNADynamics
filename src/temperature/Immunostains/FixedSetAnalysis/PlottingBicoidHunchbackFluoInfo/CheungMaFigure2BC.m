function CheungMaFigure2B(AllCompiledEmbryos)
%%
APbins = 0:0.025:1;
NumAPbins = length(APbins);

NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);
NumAPbins = length(APbins);

AllSetInfo = GetFixedSetPrefixInfo;

NumSets = length(AllSetInfo.Temperatures);

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


temps = [17.5 20 22.5 25 27.5];

%%
f = @(b,x) b(1).*exp(b(2).*x) + b(3);
UseForFit = APbins >= 0.1 & APbins <= 0.5;
x = APbins(UseForFit);
master_index = 1;
ch_index = 3;
for exp_index = 1:NumSets
    disp(['Set: ', num2str(exp_index)])
    if ~isfield(AllCompiledEmbryos{exp_index}, 'NormalizedBicoidFitInfo')
        AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo = {};
    end
    if ~isfield(AllCompiledEmbryos{exp_index}, 'BgdSubNormalizedBicoidFitInfo')
        AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo = {};
    end
    NEmbryos = size(AllCompiledEmbryos{exp_index}.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles, 1);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.fits = cell(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.Bmax = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.BmaxSE = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.RealBmax = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.lambda = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.lambdaSE = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.Intercept = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.InterceptSE = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.physlambda = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.physlambdaSE = NaN(1,NEmbryos);
    
    for i = 1:NEmbryos
        AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.RealBmax(i) = max(AllCompiledEmbryos{exp_index}.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(i,:,ch_index, master_index));
        BgdGuess = min(AllCompiledEmbryos{exp_index}.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(i,:,ch_index, master_index));
        APLength =  AllCompiledEmbryos{exp_index}.APLengths(i);
        y = AllCompiledEmbryos{exp_index}.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(i,UseForFit,ch_index, master_index);
        
        if ~isempty(y) & sum(~isnan(y)) > 7
            beta0 =[ AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.RealBmax(i) , -3, BgdGuess];
            nlm = fitnlm(x, y,f,beta0);
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.fits{i} = nlm;
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.Bmax(i) = nlm.Coefficients.Estimate(1);
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.BmaxSE(i) = nlm.Coefficients.SE(1);
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.lambda(i) = -1/nlm.Coefficients.Estimate(2);
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.lambdaSE(i) = nlm.Coefficients.SE(2)/nlm.Coefficients.Estimate(2);
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.Intercept(i) = nlm.Coefficients.Estimate(3);
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.InterceptSE(i) = nlm.Coefficients.SE(3);
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.physlambda(i) = -APLength/(nlm.Coefficients.Estimate(2));
            AllCompiledEmbryos{exp_index}.NormalizedBicoidFitInfo.physlambdaSE(i) = APLength*(nlm.Coefficients.SE(2)/nlm.Coefficients.Estimate(2));
            
        end
    end
    
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.fits = cell(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.Bmax = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.BmaxSE = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.RealBmax = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.lambda = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.lambdaSE = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.Intercept = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.InterceptSE = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.physlambda = NaN(1,NEmbryos);
    AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.physlambdaSE = NaN(1,NEmbryos);
    
    for i = 1:NEmbryos
        AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.RealBmax(i) = max(AllCompiledEmbryos{exp_index}.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(i,:,ch_index, master_index));
        BgdGuess = min(AllCompiledEmbryos{exp_index}.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(i,:,ch_index, master_index));
        APLength =  AllCompiledEmbryos{exp_index}.APLengths(i);
        y = AllCompiledEmbryos{exp_index}.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(i,UseForFit,ch_index, master_index);
        
        if ~isempty(y)& sum(~isnan(y)) > 7
            beta0 =[ AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.RealBmax(i) , -3, BgdGuess];
            nlm = fitnlm(x, y,f,beta0);
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.fits{i} = nlm;
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.Bmax(i) = nlm.Coefficients.Estimate(1);
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.BmaxSE(i) = nlm.Coefficients.SE(1);
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.lambda(i) = -1/nlm.Coefficients.Estimate(2);
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.lambdaSE(i) = nlm.Coefficients.SE(2)/nlm.Coefficients.Estimate(2);
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.Intercept(i) = nlm.Coefficients.Estimate(3);
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.InterceptSE(i) = nlm.Coefficients.SE(3);
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.physlambda(i) = -APLength/(nlm.Coefficients.Estimate(2));
            AllCompiledEmbryos{exp_index}.BgdSubNormalizedBicoidFitInfo.physlambdaSE(i) = APLength*(nlm.Coefficients.SE(2)/nlm.Coefficients.Estimate(2));
            
        end
    end
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.lambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 0.7);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector]; %#ok<AGROW>
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
MarkerStyles = {'o', 's', 'd'};

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.lambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 0.7);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
   
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.physlambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 350);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mu m)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.physlambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 350);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mum)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end

%%
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
   
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13; 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.lambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 0.7);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, ['NC13'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%
for rep_index = 1:3
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.lambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 0.7);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end


%%
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.physlambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 350);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mu m)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%
for rep_index = 1:3
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
   
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.physlambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 350);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mum)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_NC13.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end

%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
   
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.lambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 0.7);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.lambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 0.7);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.physlambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 350);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mu m)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
   
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.physlambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 350);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mum)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end

%%
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13; 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.lambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 0.7);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, ['NC13'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%


for rep_index = 1:3

    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.lambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 0.7);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/0.2))*0.2;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Lambda_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end


%%
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
 
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.physlambda(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 350);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mu m)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%


for rep_index = 1:3
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
 
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.physlambda(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 350);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/50))*50;
    plot_ymin = floor(min(ymins)/50)*50;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\lambda (\mum)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_NC13.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2C', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'PhysLambda_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end

%%


times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.Bmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    map = colormap(colorsTemperatures);
%  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.Bmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 5);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/0.1)*0.1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.RealBmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.RealBmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 5);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end

%%

    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13; 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.Bmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, ['NC13'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%


for rep_index = 1:3

    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.Bmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 5);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('\Lambda', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end


%%
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    map = colormap(colorsTemperatures);

    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.RealBmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%
for rep_index = 1:3
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.NormalizedBicoidFitInfo.RealBmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 5);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_NC13.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'NormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end

%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
   hold on
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.Bmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.Bmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 5);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
                (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2); 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.RealBmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max}(AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_AllSets_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);
end


%%
times_to_plot = 10:10:50;
time_diff = times_to_plot(2)-times_to_plot(1);
colorsDubuisTimes = (hsv(round(9*1.5))); % Colormap "jet" is another option
colorsDubuisTimes = colorsDubuisTimes(1:9,:);
colorsTemperatures = colorsDubuisTimes([9 8 5 3 1],:);
FractionalTemperatures = (1:2:(2*length(temps)-1))/(2*length(temps));
TempTickLabels = {};
MarkerStyles = {'o', 's', 'd'};
for temp_index = 1:length(temps)
    TempTickLabels{temp_index} = num2str(round(temps(temp_index),1));
end

for rep_index = 1:3
for time_index = 1:length(times_to_plot)
    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC14 & ~isnan(AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes) &...
            (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes >= times_to_plot(time_index)-time_diff/2)& (AllCompiledEmbryos{TempMatches(rep_index)}.DubuisEmbryoTimes < times_to_plot(time_index)+time_diff/2);
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.RealBmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 350);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, [num2str(times_to_plot(time_index)), 'min into NC14'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_FlippedStain_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);
end
end

%%
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    map = colormap(colorsTemperatures);
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13; 
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.Bmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx, ['NC13'], 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%


for rep_index = 1:3

    
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.Bmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 5);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_TimeBin', strrep(num2str(times_to_plot(time_index)), '.', '_'), 'm.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'Bmax_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end


%%
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
    
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];
        for rep_index = 1:length(TempMatches)
            UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
            LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.RealBmax(UseTF);
            LambdaVector = LambdaVector(LambdaVector < 5);
            TempVector = (temps(temp_index)-0.5+rep_index*0.25)*ones(1, length(LambdaVector));
           
            if ~isempty(LambdaVector)
                TempLambdas = [TempLambdas LambdaVector];
                ymaxes(TempMatches(rep_index)) = max(LambdaVector);
                ymins(TempMatches(rep_index)) = min(LambdaVector);
                scatter(TempVector, LambdaVector,100,MarkerStyles{rep_index}, 'MarkerFaceColor',...
                    colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.2);
                hold on
            end
        end
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_AllSets_MasterSet',num2str(master_index),'_NC13.png' ];
    saveas(DeltaFCFixCorrectedFig,outpath2);



%%


for rep_index = 1:3
    close all
    DeltaFCFixCorrectedFig =figure(1);
    set(DeltaFCFixCorrectedFig,'units', 'normalized', 'position',[0.01, 0.05, .6, .6]);
    set(gcf,'color','w');
    DeltaFCAx = axes(DeltaFCFixCorrectedFig);
  
    ymaxes = NaN(1, NumSets);
    ymins = NaN(1, NumSets);
    LambdaVectors = cell(1, 5);
    for temp_index = 1:length(temps)
        TempMatches = find(AllSetInfo.Temperatures == temps(temp_index));
        TempLambdas = [];

        UseTF = AllCompiledEmbryos{TempMatches(rep_index)}.Approved & AllCompiledEmbryos{TempMatches(rep_index)}.IsNC13;
        LambdaVector = AllCompiledEmbryos{TempMatches(rep_index)}.BgdSubNormalizedBicoidFitInfo.RealBmax(UseTF);
        LambdaVector = LambdaVector(LambdaVector < 5);
        TempVector = (temps(temp_index))*ones(1, length(LambdaVector));
        
        if ~isempty(LambdaVector)
            TempLambdas = [TempLambdas LambdaVector];
            ymaxes(TempMatches(rep_index)) = max(LambdaVector);
            ymins(TempMatches(rep_index)) = min(LambdaVector);
            scatter(TempVector, LambdaVector,100,'o', 'MarkerFaceColor',...
                colorsTemperatures(temp_index,:), 'MarkerEdgeColor', 'k', 'jitter', 'on', 'jitterAmount', 0.3);
            hold on
        end
        
        LambdaVectors{temp_index} = TempLambdas;
        if ~isempty(TempLambdas)
        MeanLambda = mean(TempLambdas, 'omitnan');
        sdLambda = std(TempLambdas, 'omitnan');
        plot([temps(temp_index)-1.15,temps(temp_index)-0.85],[MeanLambda MeanLambda], 'k-', 'LineWidth', 5.0);
        errorbar(temps(temp_index)-1,MeanLambda,  sdLambda, '-', 'Color', 'k');
        end
        
    end
    plot_ymax = (ceil(max(ymaxes)/1))*1;
    plot_ymin = floor(min(ymins)/1)*1;
    grid on
    hold off
    xlabel('Temperature (ºC)', 'FontSize', 20)
    ylabel('Obs. B_{max} (AU)', 'FontSize', 20)
    ylim([0,  plot_ymax])
    xlim([15 30])
    xticks([17.5 20 22.5 25 27.5])
    DeltaFCAx.YAxis.FontSize = 20;
    DeltaFCAx.XAxis.FontSize = 20;
    title(DeltaFCAx,  'NC13', 'FontSize', 24)
    DeltaFCAx.FontSize = 20;
    if ~isdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
        mkdir([AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles'])
    end
    if rep_index < 3
    outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_Replicate',num2str(rep_index),'_MasterSet',num2str(master_index),'_NC13.png' ];
    else
        outpath2 = [AllSetsProfFigPath, filesep, 'CheungMaFig2B', filesep, 'BgdSubNormalizedBootstrappedTestProfiles', filesep,...
        'RealBmax_FlippedStain_MasterSet',num2str(master_index),'_NC13.png' ];
    end
    saveas(DeltaFCFixCorrectedFig,outpath2);

end
