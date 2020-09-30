function PlotNuclearFluoProfilesAP(Prefix, varargin)
% PlotNuclearFluoProfilesAP.m
% author: Gabriella Martini
% date created: 9/26/20
% date last modified: 9/29/20

nuclear_cycles = [12, 13, 14];
min_nuc_count = 5;

if ~isempty(varargin)
    x = 1;
    while x <= length(varargin)
        switch varargin{x}
            case{'NuclearCycles'}
                 nuclear_cycles = varargin{x+1};
                 x = x+1;
            case{'MinNucBinCount'}
                min_nuc_count = varargin{x+1};
                x = x+1;
            otherwise
                error(['Flag "', varargin{x},'" not valid'])
        end
        x = x +1;
    end
end
%% 

[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);


if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledNucleiTable.mat'], 'file') 
    load([DropboxFolder,filesep,Prefix,filesep,'CompiledNucleiTable.mat']);
else
    CompiledNucleiTable = ConvertCompiledNucleiToTableArray(Prefix);
end
%Get the folders, including the default Dropbox one
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);

CompiledNucleiTable((CompiledNucleiTable.nc >= min(nuclear_cycles)) & (CompiledNucleiTable.nc <= max(nuclear_cycles)), :);
%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];
Prefix_label =  strrep(Prefix, '_5C', '.5C');
APbins = 0:APResolution:1;

colors = [0, 0.4470, 0.7410];
nc_frames = [nc9, nc10, nc11, nc12, nc13, nc14];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];
Prefix_label =  strrep(Prefix, '_5C', '.5C');
profFolder = 'APprofiles';
timeFolder = 'Timeprofiles';
% First, get the name of the folder you're using.
% For example if your folder is 'D:\photos\Info', parentFolder  would = 'D:\photos, and deepestFolder would = 'Info'.
[parentFolder deepestFolder] = fileparts([DropboxFolder,filesep,Prefix]);
% Next, create a name for a subfolder within that.
% For example 'D:\photos\Info\DATA-Info'
newSubFolder = sprintf('%s/%s/%s', parentFolder, deepestFolder, profFolder);
% Finally, create the folder if it doesn't exist already.
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end
newSubFolder = sprintf('%s/%s/%s', parentFolder, deepestFolder, timeFolder);
% Finally, create the folder if it doesn't exist already.
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

%% 

for nc=nuclear_cycles
    load([DropboxFolder,filesep,Prefix,filesep,'MeanAPposNC', num2str(nc),'_FramesNC.mat']);
    load([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAPNC', num2str(nc),'_FramesNC.mat']);
    load([DropboxFolder,filesep,Prefix,filesep,'StdFluoAPNC',num2str(nc),'_FramesNC.mat']);
    load([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC',num2str(nc),'_FramesNC.mat']);
    load([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAPNC',num2str(nc),'_FramesNC.mat']);
    load([DropboxFolder,filesep,Prefix,filesep,'DiffStdErrorAPNC',num2str(nc),'_FramesNC.mat']);
    load([DropboxFolder,filesep,Prefix,filesep,'TimeNC',num2str(nc),'_FramesNC.mat']);
    if NCTimes(1) > 0
        NCTimes = NCTimes - min(NCTimes);
    end
    
    
    %%
    close all

    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);

    eb = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors);
    hold off
    xlabel('Fraction Embryo Length')
    ylabel('Fluo (AU)')
    ylim([0, max(max(max(MeanFluoAP+StdFluoAP)),1)])
    xlim([0, 1])
    title(FrameProfAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(0), ', NC Time: ',num2str(-1), ' min' ]})
    
    for i = 1:size(MeanFluoAP, 1)
        
        use_idx = NumNucAP(i,:) >= min_nuc_count;
        if sum(use_idx) == 0 
           FrameProfAx.Children(1).XData = APbins;
           FrameProfAx.Children(1).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).XData = APbins;
           FrameProfAx.Children(2).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).YPositiveDelta = .1*ones(1, length(APbins));
           FrameProfAx.Children(2).YNegativeDelta = .1*ones(1, length(APbins));
        else
           FrameProfAx.Children(1).YData = MeanFluoAP(i, use_idx);
           FrameProfAx.Children(1).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YData = MeanFluoAP(i, use_idx);
           FrameProfAx.Children(2).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YPositiveDelta = StdFluoAP(i, use_idx)./NumNucAP(i, use_idx);
           FrameProfAx.Children(2).YNegativeDelta  = StdFluoAP(i, use_idx)./NumNucAP(i, use_idx);
            
        end
           
        try
        title(FrameProfAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(i), ', NC Time: ',num2str(NCTimes(i)), ' min' ]})
        end
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameProfFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'FluoProfile_frame', num2str(i),'.png']);

    end
    
    close all

    FrameTraceFig = figure(1);
    FrameTraceAx = axes(FrameTraceFig);
    Frames = 1:size(MeanFluoAP, 1);
    eb = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors);
    hold off
    xlabel('Time (min)')
    ylabel('Fluo (AU)')
    ylim([0, max(max(max(MeanFluoAP+StdFluoAP)),1)])
    xlim([0, max(NCTimes)])
    title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', AP bin: ',num2str(0)]})
    
    for i = 1:size(MeanFluoAP,2)
        if ~sum(~isnan(MeanFluoAP(:, i))) 
            continue 
        end
        
        use_idx = NumNucAP(:,i) >= min_nuc_count;
        if sum(use_idx) == 0 
           FrameTraceAx.Children(1).XData = NCTimes;
           FrameTraceAx.Children(1).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).XData = NCTimes;
           FrameTraceAx.Children(2).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).YPositiveDelta = .1*ones(1, length(NCTimes));
           FrameTraceAx.Children(2).YNegativeDelta = .1*ones(1, length(NCTimes));
        else
           FrameTraceAx.Children(1).YData = MeanFluoAP(use_idx,i);
           FrameTraceAx.Children(1).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YData = MeanFluoAP(use_idx,i);
           FrameTraceAx.Children(2).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YPositiveDelta = StdFluoAP(use_idx,i)./NumNucAP(use_idx, i);
           FrameTraceAx.Children(2).YNegativeDelta  = StdFluoAP(use_idx, i)./NumNucAP(use_idx, i);
            
        end
           

        title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Fraction Embryo Length: ',num2str((i-1)*APResolution)]})
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameTraceFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'FluoTrace_bin', num2str(i),'.png']);

    end
    close all
    
    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);

    eb = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors);
    hold off
    xlabel('Fraction Embryo Length')
    ylabel('Fluo (AU)')
    ylim([0, max(max(max(DiffMeanFluoAP+DiffStdErrorAP)),1)])
    xlim([0, 1])
    title(FrameProfAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(0), ', NC Time: ',num2str(-1), ' min' ]})
    
    for i = 1:size(DiffMeanFluoAP, 1)
        if ~sum(~isnan(DiffMeanFluoAP(i,:))) 
            continue 
        end
        use_idx = NumNucAP(i,:) >= min_nuc_count;
        if sum(use_idx) == 0 
           FrameProfAx.Children(1).XData = APbins;
           FrameProfAx.Children(1).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).XData = APbins;
           FrameProfAx.Children(2).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).YPositiveDelta = .1*ones(1, length(APbins));
           FrameProfAx.Children(2).YNegativeDelta = .1*ones(1, length(APbins));
        else
           FrameProfAx.Children(1).YData = DiffMeanFluoAP(i, use_idx);
           FrameProfAx.Children(1).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YData = DiffMeanFluoAP(i, use_idx);
           FrameProfAx.Children(2).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YPositiveDelta = DiffStdErrorAP(i, use_idx);
           FrameProfAx.Children(2).YNegativeDelta  = DiffStdErrorAP(i, use_idx);
            
        end
           
        %try
            title(FrameProfAx, {Prefix_label,...
                ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(i), ', NC Time: ',num2str(NCTimes(i)), ' min' ]})
        %end
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameProfFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'DiffFluoProfile_frame', num2str(i),'.png']);

    end
    
    close all

    FrameTraceFig = figure(1);
    FrameTraceAx = axes(FrameTraceFig);
    Frames = 1:size(MeanFluoAP, 1);
    eb = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors);
    hold off
    xlabel('Time (min)')
    ylabel('Fluo (AU)')
    ylim([0, max(max(max(DiffMeanFluoAP+DiffStdErrorAP)),1)])
    xlim([0, max(NCTimes)])
    title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', AP bin: ',num2str(0)]})
    
    for i = 1:size(MeanFluoAP,2)
        if ~sum(~isnan(DiffMeanFluoAP(:, i))) 
            continue 
        end
        
        use_idx = NumNucAP(:,i) >= min_nuc_count;
        if sum(use_idx) == 0
           FrameTraceAx.Children(1).XData = NCTimes;
           FrameTraceAx.Children(1).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).XData = NCTimes;
           FrameTraceAx.Children(2).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).YPositiveDelta = .1*ones(1, length(NCTimes));
           FrameTraceAx.Children(2).YNegativeDelta = .1*ones(1, length(NCTimes));
        else
           FrameTraceAx.Children(1).YData = DiffMeanFluoAP(use_idx,i);
           FrameTraceAx.Children(1).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YData = DiffMeanFluoAP(use_idx,i);
           FrameTraceAx.Children(2).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YPositiveDelta = DiffStdErrorAP(use_idx,i);
           FrameTraceAx.Children(2).YNegativeDelta  = DiffStdErrorAP(use_idx, i);
            
        end
           

        title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Fraction Embryo Length: ',num2str((i-1)*APResolution)]})
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameTraceFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'DiffFluoTrace_bin', num2str(i),'.png']);

    end
    
    try
        load([DropboxFolder,filesep,Prefix,filesep,'MeanAPposNC', num2str(nc),'_timesSinceAnaphase.mat']);
        load([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAPNC', num2str(nc),'_timesSinceAnaphase.mat']);
        MeanFluoAPV2 = MeanFluoAP;
        load([DropboxFolder,filesep,Prefix,filesep,'StdFluoAPNC',num2str(nc),'__timesSinceAnaphase.mat']);
        StdFluoAPV2 = StdFluoAP;
        load([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC',num2str(nc),'__timesSinceAnaphase.mat']);
        NumNucAPV2 = NumNucAP;
        load([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAPNC',num2str(nc),'__timesSinceAnaphase.mat']);
        DiffMeanFluoAPV2 = DiffMeanFluoAP;
        load([DropboxFolder,filesep,Prefix,filesep,'DiffStdErrorAPNC',num2str(nc),'__timesSinceAnaphase.mat']);
        DiffStdErrorAPV2 = DiffStdErrorAP;
        load([DropboxFolder,filesep,Prefix,filesep,'BinnedTimesSinceAnaphaseNC',num2str(nc),'.mat']);
        NCTimesV2 = NCTimes;
    catch
        continue
    end
    
    close all

    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);

    eb = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors);
    hold off
    xlabel('Fraction Embryo Length')
    ylabel('Fluo (AU)')
    ylim([0, max(max(max(MeanFluoAP+StdFluoAP)),1)])
    xlim([0, 1])
    title(FrameProfAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(0), ', Time since anaphase: ',num2str(-1), ' min' ]})
    
    for i = 1:size(MeanFluoAP, 1)
        
        use_idx = NumNucAP(i,:) >= min_nuc_count;
        if sum(use_idx) == 0 
           FrameProfAx.Children(1).XData = APbins;
           FrameProfAx.Children(1).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).XData = APbins;
           FrameProfAx.Children(2).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).YPositiveDelta = .1*ones(1, length(APbins));
           FrameProfAx.Children(2).YNegativeDelta = .1*ones(1, length(APbins));
        else
           FrameProfAx.Children(1).YData = MeanFluoAP(i, use_idx);
           FrameProfAx.Children(1).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YData = MeanFluoAP(i, use_idx);
           FrameProfAx.Children(2).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YPositiveDelta = StdFluoAP(i, use_idx)./NumNucAP(i, use_idx);
           FrameProfAx.Children(2).YNegativeDelta  = StdFluoAP(i, use_idx)./NumNucAP(i, use_idx);
            
        end
           
        try
        title(FrameProfAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(i), ', Time since anaphase: ',num2str(NCTimes(i)), ' min' ]})
        end
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameProfFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'FluoProfileV2_frame', num2str(i),'.png']);

    end
    %% 
    
    close all

    FrameTraceFig = figure(1);
    FrameTraceAx = axes(FrameTraceFig);
    Frames = 1:size(MeanFluoAP, 1);
    eb = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors);
    hold off
    xlabel('Time since anaphase (min)')
    ylabel('Fluo (AU)')
    ylim([0, max(max(max(MeanFluoAP+StdFluoAP)),1)])
    xlim([0, max(NCTimes)])
    title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', AP bin: ',num2str(0)]})
    
    for i = 1:size(MeanFluoAP,2)
        if ~sum(~isnan(MeanFluoAP(:, i))) 
            continue 
        end
        
        use_idx = NumNucAP(:,i) >= min_nuc_count;
        if sum(use_idx) == 0 
           FrameTraceAx.Children(1).XData = NCTimes;
           FrameTraceAx.Children(1).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).XData = NCTimes;
           FrameTraceAx.Children(2).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).YPositiveDelta = .1*ones(1, length(NCTimes));
           FrameTraceAx.Children(2).YNegativeDelta = .1*ones(1, length(NCTimes));
        else
           FrameTraceAx.Children(1).YData = MeanFluoAP(use_idx,i);
           FrameTraceAx.Children(1).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YData = MeanFluoAP(use_idx,i);
           FrameTraceAx.Children(2).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YPositiveDelta = StdFluoAP(use_idx,i)./NumNucAP(use_idx, i);
           FrameTraceAx.Children(2).YNegativeDelta  = StdFluoAP(use_idx, i)./NumNucAP(use_idx, i);
            
        end
           

        title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Fraction Embryo Length: ',num2str((i-1)*APResolution)]})
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameTraceFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'FluoTraceV2_bin', num2str(i),'.png']);

    end
    close all
    
    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);

    eb = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors);
    hold off
    xlabel('Fraction Embryo Length')
    ylabel('Fluo (AU)')
    ylim([0, max(max(max(DiffMeanFluoAP+DiffStdErrorAP)),1)])
    xlim([0, 1])
    title(FrameProfAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(0), ', Time since anaphase: ',num2str(-1), ' min' ]})
    
    for i = 1:size(DiffMeanFluoAP, 1)
        
        use_idx = NumNucAP(i,:) >= min_nuc_count;
        if sum(use_idx) == 0 
           FrameProfAx.Children(1).XData = APbins;
           FrameProfAx.Children(1).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).XData = APbins;
           FrameProfAx.Children(2).YData = zeros(1, length(APbins));
           FrameProfAx.Children(2).YPositiveDelta = .1*ones(1, length(APbins));
           FrameProfAx.Children(2).YNegativeDelta = .1*ones(1, length(APbins));
        else
           FrameProfAx.Children(1).YData = DiffMeanFluoAP(i, use_idx);
           FrameProfAx.Children(1).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YData = DiffMeanFluoAP(i, use_idx);
           FrameProfAx.Children(2).XData = MeanAP(i, use_idx);
           FrameProfAx.Children(2).YPositiveDelta = DiffStdErrorAP(i, use_idx);
           FrameProfAx.Children(2).YNegativeDelta  = DiffStdErrorAP(i, use_idx);
            
        end
           

        title(FrameProfAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(i), ', Time Since Anaphase: ',num2str(NCTimes(i)), ' min' ]})
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameProfFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'DiffFluoProfileV2_frame', num2str(i),'.png']);

    end
    
    close all

    FrameTraceFig = figure(1);
    FrameTraceAx = axes(FrameTraceFig);
    Frames = 1:size(MeanFluoAP, 1);
    eb = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    hold on 
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors);
    hold off
    xlabel('Time since anaphase (min)')
    ylabel('Fluo (AU)')
    ylim([0,max(max(max(DiffMeanFluoAP+DiffStdErrorAP)),1)])
    xlim([0, max(NCTimes)])
    title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', AP bin: ',num2str(0)]})
    
    for i = 1:size(MeanFluoAP,2)
        if ~sum(~isnan(DiffMeanFluoAP(:, i))) 
            continue 
        end
        
        use_idx = NumNucAP(:,i) >= min_nuc_count;
        if sum(use_idx) == 0 
           FrameTraceAx.Children(1).XData = NCTimes;
           FrameTraceAx.Children(1).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).XData = NCTimes;
           FrameTraceAx.Children(2).YData = zeros(1, length(NCTimes));
           FrameTraceAx.Children(2).YPositiveDelta = .1*ones(1, length(NCTimes));
           FrameTraceAx.Children(2).YNegativeDelta = .1*ones(1, length(NCTimes));
        else
           FrameTraceAx.Children(1).YData = DiffMeanFluoAP(use_idx,i);
           FrameTraceAx.Children(1).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YData = DiffMeanFluoAP(use_idx,i);
           FrameTraceAx.Children(2).XData = NCTimes(use_idx);
           FrameTraceAx.Children(2).YPositiveDelta = DiffStdErrorAP(use_idx,i);
           FrameTraceAx.Children(2).YNegativeDelta  = DiffStdErrorAP(use_idx, i);
            
        end
           

        title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Fraction Embryo Length: ',num2str((i-1)*APResolution)]})
        if sum(use_idx) == 0  
            continue 
        end
        saveas(FrameTraceFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'DiffFluoTraceV2_bin', num2str(i),'.png']);

    end
    
% 

end



load([DropboxFolder,filesep,Prefix,filesep,'MeanAPpos.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAP.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'StdFluoAP.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'NumNucAP.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAP.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'DiffStdErrorAP.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'AllTimes.mat']);

%%


close all

FrameTraceFig = figure(1);
FrameTraceAx = axes(FrameTraceFig);
Frames = 1:size(MeanFluoAP, 1);
for j=1:length(nuclear_cycles)
    xline(FrameInfo(nc_frames(nuclear_cycles(j)-8)).Time/60, 'r--')
    hold on 
end

eb = errorbar(Times, ones(1, length(Times)), .1*ones(1, length(Times)), 'vertical', 'LineStyle', 'none');
hold on 
set(eb, 'color', colors, 'LineWidth', 1);
set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
plot(Times, ones(1, length(Times)), '.-', 'Color', colors);
hold off
xlabel('Time (min)')
ylabel('Fluo (AU)')
ylim([0, max(max(max(MeanFluoAP+StdFluoAP)),1)])
xlim([0, max(Times)])
title(FrameTraceAx, {Prefix_label,...
        ['Nuclear Cycle ', num2str(nc), ', AP bin: ',num2str(0)]})

for i = 1:size(MeanFluoAP,2)
    if ~sum(~isnan(MeanFluoAP(:, i))) 
        continue 
    end

    use_idx = NumNucAP(:,i) >= min_nuc_count;
    if sum(use_idx) == 0 
       FrameTraceAx.Children(1).XData = Times;
       FrameTraceAx.Children(1).YData = zeros(1, length(Times));
       FrameTraceAx.Children(2).XData = Times;
       FrameTraceAx.Children(2).YData = zeros(1, length(Times));
       FrameTraceAx.Children(2).YPositiveDelta = .1*ones(1, length(Times));
       FrameTraceAx.Children(2).YNegativeDelta = .1*ones(1, length(Times));
    else
       FrameTraceAx.Children(1).YData = MeanFluoAP(use_idx,i);
       FrameTraceAx.Children(1).XData = Times(use_idx);
       FrameTraceAx.Children(2).YData = MeanFluoAP(use_idx,i);
       FrameTraceAx.Children(2).XData = Times(use_idx);
       FrameTraceAx.Children(2).YPositiveDelta = StdFluoAP(use_idx,i)./NumNucAP(use_idx, i);
       FrameTraceAx.Children(2).YNegativeDelta  = StdFluoAP(use_idx, i)./NumNucAP(use_idx, i);

    end


    title(FrameTraceAx, {Prefix_label,...
        ['Nuclear Cycle ', num2str(nc), ', Fraction Embryo Length: ',num2str((i-1)*APResolution)]})
    if sum(use_idx) == 0  
        continue 
    end
    saveas(FrameTraceFig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep,'FullFluoTrace_bin', num2str(i),'.png']);

end