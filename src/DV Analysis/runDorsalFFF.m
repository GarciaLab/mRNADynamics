% prefixes = {'2019-10-22-1Dg11_EfEf3_1','2019-10-23-1Dg11_EfEf3_2','2019-10-23-1Dg11_EfEf3_3',...
%    '2019-10-23-1Dg11_EfEf3_4','2019-10-23-1Dg11_EfEf3_5', '2019-10-24-1Dg11_EfEf3_6'};
% prefixes = {'2019-10-22-1Dg11_EfEf3_1','2019-10-23-1Dg11_EfEf3_2','2019-10-23-1Dg11_EfEf3_3'};

%%%%%
dataType = '1Dg-8D_FFF';
[~, prefixes, ~] = LoadMS2Sets(dataType, 'justPrefixes', 'noCompiledNuclei');
% parpool(18)
for i = 1:length(prefixes)
    if i ~=1
%            ExportDataForLivemRNA(prefixes{i},  'lowbit', 'keepTifs')
            filterMovie(prefixes{i},'highPrecision','customFilter','Difference_of_Gaussian_3D', {2,4}, 'saveAsMat', 'nogpu', 'keepPool', 'nWorkers', 1)
        TrackNuclei(prefixes{i}, 'nWorkers', 1)
    
    segmentSpots(prefixes{i},10025, 'Shadows', 1,'track', 'keepPool', 'keepProcessedData', 'nuclearMask', 'saveAsMat', 'nWorkers', 1);
%             TrackmRNADynamics(prefixes{i}, 'noRetracking');
    %         FindAPAxisFullEmbryo(prefixes{i}, 'CorrectAxis')
    %     AddParticlePosition(prefixes{i},  'yToManualAlignmentPrompt');
    %     CheckDivisionTimes(prefixes{i}, 'lazy')
    end
          nSpots = 1; fit3DGaussiansToAllSpots(prefixes{i}, nSpots, 'nWorkers', 1)
        CompileParticles(prefixes{i}, 'SkipAll', 'ApproveAll', 'minBinSize', .3, 'MinParticles', 0, 'yToManualAlignmentPrompt');
end
%
addDVStuffToSchnitzCells(dataType)
compileAllProjects(dataType)

%%%%%%


% dataTypes = {'1Dg', '1DgW_FFF', '1DgVW_FFF', '0Dg'};
% dataTypes = {'1Dg', '1Dg-5_FFF', '0Dg'};
% dataTypes = {'1Dg','1DgVW_FFF','2Dgc_FFF', '3Dg_Leica'};
% dataTypes = {'1Dg_2xDl', '1DgW_2x_Leica', '1DgW_FFF', '1Dg', '1Dg-5_FFF', '1DgVW_FFF'};
% dataTypes = {'1DgW_2x_Leica', '1DgW_FFF'};
% dataTypes = {'1Dg_2xDl', '1DgW_2x_Leica'};
% dataTypes = {'1Dg_2xDl', '1DgW_FFF', '1Dg'};
% activities = {'fraction', 'maxFluo', 'turnOn'};
% activities = {'fraction'}
dataTypes = {'1DgVW_FFF'};
for i = 1:length(dataTypes)
%     [~, ~, prefixes] = getDorsalPrefixes(dataTypes{i});
    for k = 7:length(prefixes)
%         try
%         AddParticlePosition(prefixes{k},  'yToManualAlignmentPrompt');
%         end
      nSpots = 1; fit3DGaussiansToAllSpots(prefixes{k}, nSpots)
      CompileParticles(prefixes{k}, 'SkipAll', 'ApproveAll', 'minBinSize', .3, 'MinParticles', 0, 'yToManualAlignmentPrompt');
    end
    addDVStuffToSchnitzCells(dataTypes{i})
    compileAllProjects(dataTypes{i})
end
%     addDVStuffToSchnitzCells(dataTypes{i})
%     compileAllProjects(dataTypes{i})
for j = 1:length(activities)
    for i = 1:length(dataTypes)
        if j == 1
            compileAllProjects(dataTypes{i})
        end
        plotFracByDlFluo2(dataTypes{i}, activities{j});
%         xlim([0, 3000])
        if i == 1
            ax1dg =gca;
        else
            ax = gca;
            copyPlot(ax, ax1dg);
        end
    end
%     legend(dataTypes{:});
end


%%%
% dataTypes = {'1Dg', '1DgW_FFF', '1Dg-5_FFF', '0Dg'};
% dataTypes = {'1Dg-5_FFF', '0Dg'};
% dataTypes = {'1DgVW_FFF'};
nSpots = 1;
dataTypes = {'1DgW_2x_Leica'};
for k = 1:length(dataTypes)
    [~, prefixes, ~] = LoadMS2Sets(dataTypes{k}, 'justPrefixes', 'noCompiledNuclei');
    for i= 1:length(prefixes)
        %         TrackNuclei(prefixes{i}, 'nWorkers', 1, 'retrack', 'integrate');
        %        TrackmRNADynamics(prefixes{i}, 'noRetracking');
        %        AddParticlePosition(prefixes{i},  'yToManualAlignmentPrompt');
        fit3DGaussiansToAllSpots(prefixes{i}, nSpots)
        CompileParticles(prefixes{i}, 'SkipAll', 'ApproveAll', 'minBinSize', .3, 'MinParticles', 0, 'yToManualAlignmentPrompt');
        %         CompileNuclearProtein(prefixes{i})
    end
    addDVStuffToSchnitzCells(dataTypes{k})
end


%%
