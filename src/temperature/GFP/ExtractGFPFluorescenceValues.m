
[plot_colors, ~, ~] = getColorPalettes();

ddH2OPrefix = {'2022-04-19-ddH2O_25C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-ddH2O_25C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-ddH2O_25C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-19-ddH2O_22_5C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-ddH2O_22_5C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-ddH2O_22_5C_PetriWithGlueNoEmbryos_50uW_Dish1'};

Temperatures = [17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 20, 20, 20, 20, 20, 20,...
    22.5, 22.5, 22.5, 22.5, 22.5, 22.5,...
    25, 25, 25, 25, 25, 25,27.5, 27.5, 27.5, 27.5, 27.5, 27.5,...
     20, 20, 20,22.5, 22.5, 22.5, 25, 25, 25, 27.5, 27.5, 27.5];
Dilutions = [1000, 1000, 1000, 2000, 2000, 2000, 1000, 1000, 1000, 2000, 2000, 2000,...
    1000, 1000, 1000, 2000, 2000, 2000, 1000, 1000, 1000, 2000, 2000, 2000,...
    1000, 1000, 1000, 2000, 2000, 2000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000,...
    1000000, 1000000, 1000000,  1000000, 1000000, 1000000];
Concentrations = [100, 100, 100, 50 , 50 , 50, 100, 100, 100, 50, 50, 50,...
    100, 100, 100, 50 , 50 , 50, 100, 100, 100, 50, 50, 50,...
    100, 100, 100, 50 , 50 , 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
Powers = [25, 35, 50, 25, 35, 50, 25, 35, 50, 25, 35, 50,...
    25, 35, 50, 25, 35, 50, 25, 35, 50, 25, 35, 50,...
    25, 35, 50, 25, 35, 50, 25, 35, 50, 25, 35, 50, 25, 35, 50, 25, 35, 50];

Prefixes = {'2022-04-16-EGFP_1_1000Dilution_17_5C_PetriWithGlueWithOutsideGreaseNoEmbryos_25uW_Dish1',...
'2022-04-16-EGFP_1_1000Dilution_17_5C_PetriWithGlueWithOutsideGreaseNoEmbryos_35uW_Dish1',...
'2022-04-16-EGFP_1_1000Dilution_17_5C_PetriWithGlueWithOutsideGreaseNoEmbryos_50uW_Dish1',...
'2022-04-16-EGFP_1_2000Dilution_17_5C_PetriWithGlueWithOutsideGreaseNoEmbryos_25uW_Dish1',...
'2022-04-16-EGFP_1_2000Dilution_17_5C_PetriWithGlueWithOutsideGreaseNoEmbryos_35uW_Dish1',...
'2022-04-16-EGFP_1_2000Dilution_17_5C_PetriWithGlueWithOutsideGreaseNoEmbryos_50uW_Dish1',...
'2022-04-17-EGFP_1_1000Dilution_20C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-17-EGFP_1_1000Dilution_20C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-17-EGFP_1_1000Dilution_20C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-17-EGFP_1_2000Dilution_20C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-17-EGFP_1_2000Dilution_20C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-17-EGFP_1_2000Dilution_20C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-19-EGFP_1_1000Dilution_22_5C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-EGFP_1_1000Dilution_22_5C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-EGFP_1_1000Dilution_22_5C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-19-EGFP_1_2000Dilution_22_5C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-EGFP_1_2000Dilution_22_5C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-EGFP_1_2000Dilution_22_5C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-19-EGFP_1_1000Dilution_25C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-EGFP_1_1000Dilution_25C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-EGFP_1_1000Dilution_25C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-19-EGFP_1_2000Dilution_25C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-EGFP_1_2000Dilution_25C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-EGFP_1_2000Dilution_25C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-18-EGFP_1_1000Dilution_27_5C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-18-EGFP_1_1000Dilution_27_5C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-18-EGFP_1_1000Dilution_27_5C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-18-EGFP_1_2000Dilution_27_5C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-18-EGFP_1_2000Dilution_27_5C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-18-EGFP_1_2000Dilution_27_5C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-17-ddH2O_20C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-17-ddH2O_20C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-17-ddH2O_20C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-19-ddH2O_22_5C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-ddH2O_22_5C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-ddH2O_22_5C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-19-ddH2O_25C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-19-ddH2O_25C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-19-ddH2O_25C_PetriWithGlueNoEmbryos_50uW_Dish1',...
'2022-04-18-ddH2O_27_5C_PetriWithGlueNoEmbryos_25uW_Dish1',...
'2022-04-18-ddH2O_27_5C_PetriWithGlueNoEmbryos_35uW_Dish1',...
'2022-04-18-ddH2O_27_5C_PetriWithGlueNoEmbryos_50uW_Dish1'...
};
%%
AllTIFs = {};
AllMaxTIFs = {};
MeanValues = NaN(36, 3);
MeanValues2 = NaN(36, 3);
MedianValues = NaN(36, 3);
MeanValueSlices = NaN(36, 3, 30);
for i = 1:length(Prefixes)
    disp(['i = ', num2str(i)])
    Prefix = Prefixes{i};
    rootFolder = '';
  [rawDataPath, ProcPath, DropboxFolder, ~, PreProcPath, rawDataFolder, Prefix, ExperimentType, Channel1, Channel2, ~,...
      Channel3] = readMovieDatabase(Prefix,'rootFolder', rootFolder);
  
  [LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);
  ySize = size(LIFImages{1}{1,1}, 1);
  xSize = size(LIFImages{1}{1,1}, 2);
  MovieTIFs = cell(1, 3);
  MaxTIFs = cell(1, 3);
  % loop through each series
  for seriesIndex = 1:size(LIFImages, 1)
      NumSlices = size(LIFImages{seriesIndex,1}, 1);
      TIFstack = zeros(ySize, xSize,NumSlices, 'uint16');
      for imageIndex = 1:NumSlices
          
          TIFstack(:,:,imageIndex) = LIFImages{seriesIndex,1}{imageIndex,1};
          MeanValueSlices(i, seriesIndex,  imageIndex) = mean(TIFstack(:,:,imageIndex), 'all');
      end
      MovieTIFs{seriesIndex} = TIFstack;
      MaxTIFs{seriesIndex} = max(TIFstack,[], 3);
      MeanValues(i, seriesIndex) = mean(MaxTIFs{seriesIndex} , 'all');
      MedianValues(i, seriesIndex) = median(TIFstack, 'all');
      MeanValues2(i, seriesIndex) = max(TIFstack, [],'all');
  end
  AllTIFs{i} = MovieTIFs;
  AllMaxTIFs{i} = MaxTIFs;
      
end
%%

MeanFluos = mean(MeanValues, 2).';
StdFluos = std(MeanValues,0, 2).';
MeanFluos2 = [MeanValues(:,1).' MeanValues(:,2).' MeanValues(:,3).' ];
MeanFluos3 = mean(MeanValues2, 2).';
MeanTemperatures = [Temperatures Temperatures Temperatures];
MeanDilutions = [Dilutions Dilutions Dilutions];
MeanConcentrations = [Concentrations Concentrations Concentrations];
MeanPowers = [Powers Powers Powers];

%%
MeanTop1Percent = NaN(30, 3);
MeanTop5Percent = NaN(30, 3);
MeanMiddle50Percent = NaN(30, 3);
for i = 1:30
    for j = 1:3
        TIFvals = AllTIFs{i}{j}(:).';
       p1 = prctile(TIFvals,99);
       MeanTop1Percent(i, j) = mean(double(TIFvals(TIFvals >= p1)));
       p5 = prctile(TIFvals,95);
       MeanTop5Percent(i, j) = mean(double(TIFvals(TIFvals >= p5)));
       p25 = prctile(TIFvals,25);
       p75 = prctile(TIFvals,75);
       MeanMiddle50Percent(i, j) = mean(double(TIFvals(TIFvals >= p25 & TIFvals <= 75)));
    end
end

MeanFluos3 = [MeanMiddle50Percent(:,1).' MeanMiddle50Percent(:,2).' MeanMiddle50Percent(:,3).' ];

%% Calculate ddH2O background fluorescence 
power_levels =  [25, 35, 50];
BackgroundFluos = NaN(1, 3);
for i = 1:3
BackgroundFluos(i) = mean(MeanValues(Concentrations == 0 & Powers == power_levels(i), 2), 'all');
end
%%
UniqueTemperatures = [27.5, 25, 22.5, 20, 17.5];
fits = {};
PredictedValues = NaN(9, 5);
PredictedCIs = NaN(9, 5);
PredictedRatios = NaN(9, 5);
PredictedV2s = NaN(6, 5);
PredictedR2s = NaN(6, 5);
fits2 = {};
PredictedV3s = NaN(6, 5);
PredictedR3s = NaN(6, 5);
fits3 = {};
power_levels = [25, 35, 50];
counter = 1;
figure(1)
for conc_level = [50, 100]
    for  i = 1:length(power_levels)
        power_level = power_levels(i);
        fits{end+1} = fitlm(Temperatures(Concentrations == conc_level & Powers == power_level),...
            MeanFluos(Concentrations == conc_level & Powers == power_level));% - MeanFluos(Dilutions == 1000000 & Powers == power_level& Temperatures > 17.5));
        PredictedValues(counter,:) = predict(fits{end}, UniqueTemperatures.').';
        PredictedRatios(counter,:) = PredictedValues(counter,:)/PredictedValues(counter,2);
%         if dilution_level <= 2000
%             fits2{end+1} = fitlm(Temperatures(Concentrations == conc_level & Powers == power_level & Temperatures > 17.5),...
%             MeanFluos(Concentrations == conc_level & Powers == power_level & Temperatures > 17.5)-MeanFluos(Concentrations == 0 & Powers == power_level & Temperatures > 17.5));
%             PredictedV2s(counter,:) = predict(fits2{end}, ltm.UniqueTemperatures.').';
%             PredictedR2s(counter,:) = PredictedV2s(counter,:)/PredictedV2s(counter,2);
%             fits3{end+1} = fitlm(Temperatures(Concentrations == conc_level & Powers == power_level),...
%             MeanFluos(Concentrations == conc_level & Powers == power_level)/MeanFluos(Concentrations == conc_level& Powers == power_level & Temperatures ==25));
%             PredictedV3s(counter,:) = predict(fits3{end}, ltm.UniqueTemperatures.').';
%             PredictedR3s(counter,:) = PredictedV3s(counter,:)/PredictedV3s(counter,2);
%    
%         end
                 

           scatter(Temperatures(Concentrations == conc_level & Powers == power_level),...
               MeanFluos(Concentrations == conc_level & Powers == power_level),...
               50, 'filled', 'MarkerFaceColor', plot_colors(counter,:));
            hold on
        counter = counter + 1;
    end
end
hold off
%%
counter = 1;
figure(2)
for conc_level = [0,  50, 100]
    for T_level = [17.5, 20, 22.5, 25, 27.5]
        
           scatter(Powers(Concentrations == conc_level & Temperatures == T_level),...
               MeanFluos(Concentrations == conc_level & Temperatures == T_level),...
               50, 'filled', 'MarkerFaceColor', plot_colors(mod(counter, 20),:));
            hold on
        counter = counter + 1;
        if mod(counter, 20) == 0
            counter = counter + 1;
        end
    end
end
hold off

counter = 1;
figure(3)
for T_level = [17.5, 20, 22.5, 25, 27.5]
    for power_level = [25 35 50]
        
           scatter(Dilutions(Temperatures == T_level & Powers == power_level & Concentrations > 0),...
               MeanFluos(Temperatures == T_level & Powers == power_level & Concentrations > 0),...
               50, 'filled', 'MarkerFaceColor', plot_colors(mod(counter, 20),:));
            hold on
        counter = counter + 1;
        if mod(counter, 20) == 0
            counter = counter + 1;
        end
    end
end
hold off

%%
tbl = table(Temperatures.',Powers.',Concentrations.', MeanFluos.',...
    'VariableNames',{'Temperature','Power','Concentration', 'MeanFluo'});
tbl2 = table(MeanTemperatures.',MeanPowers.',MeanConcentrations.', MeanFluos2.',...
    'VariableNames',{'Temperature','Power','Concentration', 'MeanFluo'});
T_levels = [17.5, 20, 22.5,25,27.5];
mdl3 = {};
for i = 1:5
mdl3{i}= fitlm([Powers(Temperatures == T_levels(i)).' Concentrations(Temperatures == T_levels(i)).'],  MeanFluos(Temperatures == T_levels(i)).',...
    'interactions', 'VarNames', {'P', 'C', 'Fluo'});
end
%mdl4 = removeTerms(mdl3, 'C+ C:C^2+T:T + T:T^2 + P:P^2 +T:P:P + P:P:C+T:C:C + P:C:C+T:C+P:C')

%% Make Plots for Hernan 
FluoParamPath = 'S:/Gabriella\Dropbox\TemperatureParameters\GFP\';
load([FluoParamPath filesep 'FluoCoefficients.mat']);
outdir = 'S:/Gabriella/Dropbox\GMPlots/TemperatureAnalysis/Hunchback/MS2V5/ThesisPlotting/GFPFluoPlots';
if ~isdir(outdir)
    mkdir(outdir)
end

%%

eb = cell(1, 3);
prof = cell(1,3);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
    eb{power_level_idx} = errorbar(Temperatures(Powers == power_level & Concentrations == 100),...
        MeanFluos(Powers == power_level & Concentrations == 100),...
        StdFluos(Powers == power_level & Concentrations == 100),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
        prof{power_level_idx} = plot(Temperatures(Powers == power_level & Concentrations == 100),...
            MeanFluos(Powers == power_level & Concentrations == 100),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end


hold off
xlabel('Temperature (ºC)', 'FontSize', 10)

xlim([16, 29])
ylim([0 100])
xticks([17.5, 20, 22.5, 25, 27.5])
yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('1000-fold Dilution')
legend({'25 \muW', '35 \muW', '50\muW'})
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '1000foldDilutionDataVsTemperature.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');

%%

eb = cell(1, 3);
prof = cell(1,3);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
    eb{power_level_idx} = errorbar(Temperatures(Powers == power_level & Concentrations == 50),...
        MeanFluos(Powers == power_level & Concentrations == 50),...
        StdFluos(Powers == power_level & Concentrations == 50),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
        prof{power_level_idx} = plot(Temperatures(Powers == power_level & Concentrations == 50),...
            MeanFluos(Powers == power_level & Concentrations == 50),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end


hold off
xlabel('Temperature (ºC)', 'FontSize', 10)

xlim([16, 29])
ylim([0 100])
xticks([17.5, 20, 22.5, 25, 27.5])
yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('2000-fold Dilution')
legend({'25 \muW', '35 \muW', '50\muW'})
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '2000foldDilutionDataVsTemperature.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');

%%

eb = cell(1, 5);
prof = cell(1,5);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:5
    power_level = UniqueTemperatures(power_level_idx);
    eb{power_level_idx} = errorbar(Powers(Temperatures == power_level & Concentrations == 100),...
        MeanFluos(Temperatures == power_level & Concentrations == 100),...
        StdFluos(Temperatures == power_level & Concentrations == 100),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:5
    power_level = UniqueTemperatures(power_level_idx);
        prof{power_level_idx} = plot(Powers(Temperatures == power_level & Concentrations == 100),...
            MeanFluos(Temperatures == power_level & Concentrations == 100),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', '-', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end

prof{6} = plot(power_levels,...
            BackgroundFluos,...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(6,:),...
            'MarkerSize', 5, 'LineStyle', '-', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);

hold off
xlabel('Power (\muW)', 'FontSize', 10)
xlim([20, 60])
ylim([0 100])
xticks([20 25 30 35 40 45 50, 55, 60])
yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('1000-fold Dilution')
legend({'27.5ºC', '25ºC', '22.5ºC', '20ºC', '17.5ºC', 'ddH2O'}, 'location', 'northwest')
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '1000foldDilutionDataVsPower.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');

%%

eb = cell(1, 5);
prof = cell(1,5);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:5
    power_level = UniqueTemperatures(power_level_idx);
    eb{power_level_idx} = errorbar(Powers(Temperatures == power_level & Concentrations == 50),...
        MeanFluos(Temperatures == power_level & Concentrations == 50),...
        StdFluos(Temperatures == power_level & Concentrations == 50),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:5
    power_level = UniqueTemperatures(power_level_idx);
        prof{power_level_idx} = plot(Powers(Temperatures == power_level & Concentrations == 50),...
            MeanFluos(Temperatures == power_level & Concentrations == 50),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', '-', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end

prof{6} = plot(power_levels,...
            BackgroundFluos,...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(6,:),...
            'MarkerSize', 5, 'LineStyle', '-', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
hold off
xlabel('Power (\muW)', 'FontSize', 10)
xlim([20, 60])
ylim([0 100])
xticks([20 25 30 35 40 45 50, 55, 60])
yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('2000-fold Dilution')
legend({'27.5ºC', '25ºC', '22.5ºC', '20ºC', '17.5ºC', 'ddH2O'}, 'location', 'northwest')
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '2000foldDilutionDataVsPower.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');


%%
concentrations = [50 100];
power_levels = [25 35 50];
errorbar_values = NaN(6, 5);
ratio_errorbar_values = NaN(6, 5);
for conc_index = 1:2

    for power_index = 1:3
        fit_index = (conc_index-1)*3 + power_index;
        errorbar_values(fit_index, :) = sqrt(fits{fit_index}.Coefficients.SE(1)^2+UniqueTemperatures*fits{1}.Coefficients.SE(2)^2);
        ratio_errorbar_values(fit_index,:) = sqrt(errorbar_values(fit_index, :).^2/(PredictedValues(fit_index, 2)^2)+ ...
            (errorbar_values(fit_index, 2)^2)*(PredictedValues(fit_index, :).^2)./(PredictedValues(fit_index, 2)^4));
    end
end

ratio_errorbar_values(:,2) = 0;


%%


eb = cell(1, 3);
prof = cell(1,3);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
    eb{power_level_idx} = errorbar(UniqueTemperatures,...
        PredictedRatios(power_level_idx,:),...
        ratio_errorbar_values(power_level_idx,:),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
        prof{power_level_idx} = plot(UniqueTemperatures,...
        PredictedRatios(power_level_idx,:),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end

prof{4} = plot(FluoInfo.Temperatures,...
        1./FluoInfo.FluoCoeffs,...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(4,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);


hold off
xlabel('Temperature (ºC)', 'FontSize', 10)

xlim([16, 29])
ylim([0 2.5])
xticks([17.5, 20, 22.5, 25, 27.5])
%yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('1000-fold Dilution')
legend({'25 \muW', '35 \muW', '50\muW', 'Zhang et al. 2009'})
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '1000foldDilutionDataRatios.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');


eb = cell(1, 3);
prof = cell(1,3);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
    eb{power_level_idx} = errorbar(UniqueTemperatures,...
        PredictedRatios(power_level_idx+3,:),...
        ratio_errorbar_values(power_level_idx+3,:),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
        prof{power_level_idx} = plot(UniqueTemperatures,...
        PredictedRatios(power_level_idx+3,:),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end

prof{4} = plot(FluoInfo.Temperatures,...
        1./FluoInfo.FluoCoeffs,...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(4,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);


hold off
xlabel('Temperature (ºC)', 'FontSize', 10)

xlim([16, 29])
ylim([0 2.5])
xticks([17.5, 20, 22.5, 25, 27.5])
%yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('2000-fold Dilution')
legend({'25 \muW', '35 \muW', '50\muW', 'Zhang et al. 2009'})
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '2000foldDilutionDataRatios.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');

%%
concentrations = [50 100];
power_levels = [25 35 50];
errorbar_values = NaN(6, 5);
ratio_errorbar_values = NaN(6, 5);
for conc_index = 1:2

    for power_index = 1:3
        fit_index = (conc_index-1)*3 + power_index;
        errorbar_values(fit_index, :) = sqrt(fits{fit_index}.Coefficients.SE(1)^2+UniqueTemperatures*fits{1}.Coefficients.SE(2)^2);
        ratio_errorbar_values(fit_index,:) = sqrt(errorbar_values(fit_index, :).^2/((PredictedValues(fit_index, 2)-BackgroundFluos(power_index))^2)+ ...
            (errorbar_values(fit_index, 2)^2)*((PredictedValues(fit_index, :)-BackgroundFluos(power_index)).^2)./((PredictedValues(fit_index, 2)-BackgroundFluos(power_index)).^4));
    end
end

ratio_errorbar_values(:,2) = 0;
PredictedRatios2 = NaN(size(PredictedRatios));
for i = 1:6
    bgd_index = i;
    if bgd_index > 3
        bgd_index = i-3;
    end
    PredictedRatios2(i,:) = (PredictedValues(i,:)-BackgroundFluos(bgd_index))/(PredictedValues(i,2)-BackgroundFluos(bgd_index));
end


%%


eb = cell(1, 3);
prof = cell(1,3);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
    eb{power_level_idx} = errorbar(UniqueTemperatures,...
        PredictedRatios2(power_level_idx,:),...
        ratio_errorbar_values(power_level_idx,:),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
        prof{power_level_idx} = plot(UniqueTemperatures,...
        PredictedRatios2(power_level_idx,:),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end

prof{4} = plot(FluoInfo.Temperatures,...
        1./FluoInfo.FluoCoeffs,...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(4,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);


hold off
xlabel('Temperature (ºC)', 'FontSize', 10)

xlim([16, 29])
ylim([0 3])
xticks([17.5, 20, 22.5, 25, 27.5])
%yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('1000-fold Dilution')
legend({'25 \muW', '35 \muW', '50\muW', 'Zhang et al. 2009'})
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '1000foldDilutionDataRatiosBackgroundSubtracted.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');

%%
eb = cell(1, 3);
prof = cell(1,3);
close all

FrameProfFig = figure(1);
set(FrameProfFig,'units', 'normalized', 'position',[0.05, 0.05, .3, .4]);
set(gcf,'color','w');
FrameProfAx = axes(FrameProfFig);

power_levels = [25 35 50];

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
    eb{power_level_idx} = errorbar(UniqueTemperatures,...
        PredictedRatios2(power_level_idx+3,:),...
        ratio_errorbar_values(power_level_idx+3,:),...
        'vertical', 'LineStyle', 'none',...
        'CapSize', 0, 'Color', 'k');
    hold on
    set(eb{power_level_idx}, 'color', 'k',...%colors(temp_idx,:),
        'LineWidth', 1);
    set(get(get(eb{power_level_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %set(eb{power_level_idx},'Visible','off'); %'off' or 'on'
end

for  power_level_idx = 1:3
    power_level = power_levels(power_level_idx);
        prof{power_level_idx} = plot(UniqueTemperatures,...
        PredictedRatios2(power_level_idx+3,:),...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(power_level_idx,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);
end

prof{4} = plot(FluoInfo.Temperatures,...
        1./FluoInfo.FluoCoeffs,...
            'o',...
            'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', plot_colors(4,:),...
            'MarkerSize', 5, 'LineStyle', 'none', 'Color', [0 0 0],...%colors(temp_idx,:),...
            'LineWidth', 1);


hold off
xlabel('Temperature (ºC)', 'FontSize', 10)

xlim([16, 29])
ylim([0 3])
xticks([17.5, 20, 22.5, 25, 27.5])
%yticks([0 20 40 60 80 100])
ylabel('Fluorescence (AU)', 'FontSize', 10)
title('2000-fold Dilution')
legend({'25 \muW', '35 \muW', '50\muW', 'Zhang et al. 2009'})
%

FrameProfAx.XAxis.FontSize = 12;
FrameProfAx.YAxis.FontSize = 12;

outpath = [outdir, filesep, '2000foldDilutionDataRatiosBackgroundSubtracted.pdf'];
exportgraphics(FrameProfFig,outpath, 'ContentType','vector');














 
 