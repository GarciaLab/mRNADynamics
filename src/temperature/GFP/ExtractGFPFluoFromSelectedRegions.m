Prefixes = {'2022-04-14-EGFP_1_1000Dilution_25C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide1',...
'2022-04-14-EGFP_1_1000Dilution_25C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide1Pass2',...
'2022-04-14-EGFP_1_1000Dilution_27_5C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide1',...
'2022-04-14-EGFP_1_1000Dilution_22_5C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide1',...
'2022-04-14-EGFP_1_1000Dilution_20C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide1',...
'2022-04-14-EGFP_1_1000Dilution_17_5C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide1',...
'2022-04-15-EGFP_1_1000Dilution_17_5C_HydrophilicWithGlueWithOutsideGrease_25uW_Slide1',...
'2022-04-15-EGFP_1_1000Dilution_17_5C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide1',...
'2022-04-15-EGFP_1_1000Dilution_17_5C_HydrophilicWithGlueWithOutsideGrease_50uW_Slide1',...
'2022-04-15-EGFP_1_1000Dilution_17_5C_HydrophilicWithGlueWithOutsideGrease_25uW_Slide2',...
'2022-04-15-EGFP_1_1000Dilution_17_5C_HydrophilicWithGlueWithOutsideGrease_35uW_Slide2',...
'2022-04-15-EGFP_1_1000Dilution_17_5C_HydrophilicWithGlueWithOutsideGrease_50uW_Slide2'};

Concentrations2 = 100*ones(1, length(Prefixes));
Powers2 = [35 35 35 35 35 35 35 35 25 35 50 25 35 50];
Tempertures2 = [25 25  27.5 22.5 20 17.5 17.5 17.5 17.5 17.5 17.5 17.5];
AllTIFsV2 = {};
AllMaxTIFsV2 = {};
MeanValuesV2 = NaN(length(Prefixes), 3);
MeanValueSlicesV2 = NaN(length(Prefixes), 3, 30);
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
  MaxTIFs = {};
  % loop through each series
  for seriesIndex = 1:size(LIFImages, 1)
      NumSlices = size(LIFImages{seriesIndex,1}, 1);
      TIFstack = zeros(ySize, xSize,NumSlices, 'uint16');
      for imageIndex = 1:NumSlices
          
          TIFstack(:,:,imageIndex) = LIFImages{seriesIndex,1}{imageIndex,1};
          MeanValueSlicesV2(i, seriesIndex,  imageIndex) = mean(TIFstack(:,:,imageIndex), 'all');
      end
      MovieTIFs{seriesIndex} = TIFstack;
      MaxTIFs{seriesIndex} = max(TIFstack,[], 3);
      MeanValuesV2(i, seriesIndex) = mean(MaxTIFs{seriesIndex} , 'all');
  end
  AllTIFsV2{i} = MovieTIFs;
  AllMaxTIFsV2{i} = MaxTIFs;
      
end


MeanFluosV2 = mean(MeanValuesV2, 2).';



%%
MeanValuesV2 = NaN(12, 5);
RegionsToUse = NaN(length(Prefixes), 3, 4);
RegionsToUse(1, 1, :) = [1 250 1 856];
RegionsToUse(2, 1, :) = [1 250 1 856];
RegionsToUse(3, 1, :) = [1 250 1 856];
RegionsToUse(4, 1, :) = [1 250 1 856];
RegionsToUse(5, 1, :) = [1 250 1 856];
RegionsToUse(6, 1, :) = [1 250 1 856];
RegionsToUse(7, 1, :) = [1 856 1 250];
RegionsToUse(7, 2, :) = [500 856 1 856];
RegionsToUse(7, 3, :) = [600 856 1 856];
RegionsToUse(7, 4, :) = [1 856 500 856];
RegionsToUse(7, 5, :) = [500 856 1 856];
RegionsToUse(8, 1, :) = [1 856 1 250];
RegionsToUse(8, 2, :) = [500 856 1 856];
RegionsToUse(8, 3, :) = [600 856 1 856];
RegionsToUse(8, 4, :) = [1 856 500 856];
RegionsToUse(8, 5, :) = [500 856 1 856];
RegionsToUse(9, 1, :) = [1 856 1 250];
RegionsToUse(9, 2, :) = [500 856 1 856];
RegionsToUse(9, 3, :) = [600 856 1 856];
RegionsToUse(9, 4, :) = [1 856 500 856];
RegionsToUse(9, 5, :) = [500 856 1 856];
RegionsToUse(10, 1, :) = [200 856 500 856];
RegionsToUse(10, 2, :) = [500 856 1 856];
RegionsToUse(10, 3, :) = [1 856 500 856];
RegionsToUse(10, 4, :) = [1 300 1 856];
RegionsToUse(10, 5, :) = [1 300 1 856];
RegionsToUse(11, 1, :) = [200 856 500 856];
RegionsToUse(11, 2, :) = [500 856 1 856];
RegionsToUse(11, 3, :) = [1 856 500 856];
RegionsToUse(11, 4, :) = [1 300 1 856];
RegionsToUse(11, 5, :) = [1 300 1 856];
RegionsToUse(12, 1, :) = [200 856 500 856];
RegionsToUse(12, 2, :) = [500 856 1 856];
RegionsToUse(12, 3, :) = [1 856 500 856];
RegionsToUse(12, 4, :) = [1 300 1 856];
RegionsToUse(12, 5, :) = [1 300 1 856];


for i = 1:length(Prefixes)
for j = 1:length(AllMaxTIFsV2{i})
    if isempty(AllMaxTIFsV2{i}{j})
        continue
    end
    GFPImage = AllMaxTIFsV2{i}{j};
    if all(~isnan(RegionsToUse(i, j,:)))
        MeanValuesV2(i, j) = mean(GFPImage(RegionsToUse(i,j,1):RegionsToUse(i,j,2),RegionsToUse(i,j,3):RegionsToUse(i,j,4)), 'all');
    end

end
end

%%
