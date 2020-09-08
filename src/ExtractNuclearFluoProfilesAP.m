function ExtractNuclearFluoProfilesAP(Prefix, varargin)
% ExtractNuclearFluoProfilesAP.m
% author: Gabriella Martini
% date created: 8/27/20
% date last modified: 8/27/20

nuclear_cycles = [9, 10, 11, 12, 13, 14];
quantile_cutoff = .05;
min_nuc_count = 5;

if ~isempty(varargin)
    x = 1;
    while x <= length(varargin)
        switch varargin{x}
            case{'NuclearCycles'}
                 nuclear_cycles = varargin{x+1};
                 x = x+1;
            case{'QuantileCutoff'}
                quantile_cutoff = varargin{x+1};
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
nc_frames = [nc9, nc10, nc11, nc12, nc13, nc14];
nc_frames = nc_frames(nuclear_cycles - 8);
colors = [0, 0.4470, 0.7410];
for nc=nuclear_cycles
    cnt  = CompiledNucleiTable(CompiledNucleiTable.nc == nc, :);
    cnt  = cnt(cnt.Fluo > 0,:);
    cnt  = cnt(cnt.Flag1 == false,:);
    cnt  = cnt(cnt.Flag2 == false,:);
    cnt  = cnt(cnt.Flag3 == false,:);
    cnt  = cnt(cnt.Flag4 == false,:);
    if isempty(cnt)
        continue
    end
    %%
    close all
    numFrames = max(cnt.FrameNC)+1;
    numBins = 1/APResolution+1;
    MeanFluoAP = zeros(numFrames, numBins);
    StdFluoAP = zeros(numFrames, numBins);
    NumNucAP = zeros(numFrames, numBins);
    DiffMeanFluoAP = zeros(numFrames, numBins);
    DiffStdFluoAP = zeros(numFrames, numBins);
    APbins = 0:APResolution:1;
    fc = 1;
    NCFrames = unique(cnt.FrameNC);
    NCTimes = [];
    for i=1:length(NCFrames)
        NCTimes(i) =  mean(cnt.TimeNC(cnt.FrameNC == NCFrames(i)));
    end
    %NCTimes = 
    for i = 1:length(NCFrames)
        f = NCFrames(i);
        t = NCTimes(i);
        close all
        for b = 0:(numBins-1)
            APbin = b*APResolution;
            sub_table = cnt((cnt.FrameNC == f) & (abs(cnt.APbin- APbin) < 0.0001), :);
            if size(sub_table, 5) >= 10
                cutoff = quantile(sub_table.Fluo2, quantile_cutoff);
                fluoVector = sub_table.Fluo2(sub_table.Fluo2 >= cutoff);
            else
                fluoVector = sub_table.Fluo2;
            end
            MeanFluo = mean(fluoVector);
            StdFluo = std(fluoVector);
            MeanFluoAP(i, b+1) = MeanFluo;
            StdFluoAP(i, b+1) = StdFluo;
            NumNucAP(i, b+1) = size(fluoVector, 1);
        end
        fig = figure(fc);
        ax = axes('Parent', fig);
        use_idx = NumNucAP(i,:) >= min_nuc_count;
        if all(isnan(MeanFluoAP(i,use_idx)))
            continue
        end
        minFluoBin= find(MeanFluoAP(i,:) == min(MeanFluoAP(i,use_idx)));
        
        DiffMeanFluoAP(i,:) = MeanFluoAP(i,:)- MeanFluoAP(i,minFluoBin(1));
        DiffStdFluoAP(i,:) = sqrt(StdFluoAP(i,:).^2+StdFluoAP(i,minFluoBin( 1)).^2);
        eb = errorbar(APbins(use_idx), MeanFluoAP(i, use_idx), StdFluoAP(i, use_idx)./NumNucAP(i, use_idx), 'vertical', 'LineStyle', 'none');
        hold on
        set(eb, 'color', colors, 'LineWidth', 1);
        set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        plot(APbins(use_idx), MeanFluoAP(i, use_idx), '.-', 'Color', colors);
        
        %histogram(sub_table.Fluo)
        xlabel('Fraction Embryo Length')
        ylabel('Fluo (AU)')
        ylim([0, (max(cnt.Fluo2)-min(cnt.Fluo2))*1.1])
        xlim([0, 0.95])
        %CurrentFrame = sub_table.Frame(1);

    %     ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
    %             FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
    %     imshow(ImageHis,'DisplayRange',[]);
        title(ax, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(f), ', NC Time: ',num2str(t), ' min' ]})
    %     hold on 


    %     scatter(sub14.xPos, sub14.yPos,40, sub14.Fluo, 'filled')
    %     plot(sub14.MedianDV, sub14.Fluo,'o')


        %ylabel('Fluo (AU)')
        %xlabel('Median DV (pixels)')
        saveas(fig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'NC', num2str(nc), 'FluoProfile_frame', num2str(f),'.png']);
        fc = fc+1;
        hold off
        close all
    end
% 

    close all

    save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAPNC', num2str(nc),'_Frames.mat'],...
            'MeanFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAPNC',num2str(nc),'_Frames.mat'],...
            'StdFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC',num2str(nc),'_Frames.mat'],...
            'NumNucAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAPNC',num2str(nc),'_Frames.mat'],...
            'DiffMeanFluoAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffStdFluoAPNC',num2str(nc),'_Frames.mat'],...
            'DiffStdFluoAP');

end

for nc=nuclear_cycles
    cnt  = CompiledNucleiTable(CompiledNucleiTable.nc == nc, :);
    cnt  = cnt(cnt.Fluo > 0,:);
    cnt  = cnt(cnt.Flag1 == false,:);
    cnt  = cnt(cnt.Flag2 == false,:);
    cnt  = cnt(cnt.Flag3 == false,:);
    cnt  = cnt(cnt.Flag4 == false,:);
    if isempty(cnt)
        continue
    end
    %%
    close all
    NCTimes = unique(round(cnt.TimeNC, 1));
    cnt.RoundedTimeNC = round(cnt.TimeNC, 1);
    numTimes = length(NCTimes);
    numBins = 1/APResolution+1;
    MeanFluoAP = zeros(numTimes, numBins);
    StdFluoAP = zeros(numTimes, numBins);
    NumNucAP = zeros(numTimes, numBins);
    DiffMeanFluoAP = zeros(numTimes, numBins);
    DiffStdFluoAP = zeros(numTimes, numBins);
    
    APbins = 0:APResolution:1;
    fc = 1;
    
    %NCTimes = 
    for i = 1:numTimes
        t = NCTimes(i);
        close all
        for b = 0:(numBins-1)
            APbin = b*APResolution;
            sub_table = cnt((cnt.RoundedTimeNC == t) & (abs(cnt.APbin- APbin) < 0.0001), :);
            if size(sub_table, 5) >= 10
                cutoff = quantile(sub_table.Fluo2, quantile_cutoff);
                fluoVector = sub_table.Fluo2(sub_table.Fluo2 >= cutoff);
            else
                fluoVector = sub_table.Fluo2;
            end
            MeanFluo = mean(fluoVector);
            StdFluo = std(fluoVector);
            MeanFluoAP(i, b+1) = MeanFluo;
            StdFluoAP(i, b+1) = StdFluo;
            NumNucAP(i, b+1) = size(fluoVector, 1);
        end
        use_idx = NumNucAP(i,:) >= min_nuc_count;
        if all(isnan(MeanFluoAP(i,use_idx)))
            continue
        end
        minFluoBin= find(MeanFluoAP(i,:) == min(MeanFluoAP(i,use_idx)));
        
        DiffMeanFluoAP(i,:) = MeanFluoAP(i,:)- MeanFluoAP(i,minFluoBin(1));
        DiffStdFluoAP(i,:) = sqrt(StdFluoAP(i,:).^2+StdFluoAP(i,minFluoBin( 1)).^2);

    end


    close all

    save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAPNC', num2str(nc),'_Times.mat'],...
            'MeanFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAPNC',num2str(nc),'_Times.mat'],...
            'StdFluoAP');

    save([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC',num2str(nc),'_Times.mat'],...
            'NumNucAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAPNC',num2str(nc),'_Times.mat'],...
            'DiffMeanFluoAP');
    save([DropboxFolder,filesep,Prefix,filesep,'DiffStdFluoAPNC',num2str(nc),'_Times.mat'],...
            'DiffStdFluoAP');
    save([DropboxFolder,filesep,Prefix,filesep,'TimesNC',num2str(nc),'.mat'],...
        'NCTimes');

end


%%
close all
cnt  = CompiledNucleiTable(CompiledNucleiTable.Fluo > 0,:);
cnt  = cnt(cnt.Flag1 == false,:);
cnt  = cnt(cnt.Flag2 == false,:);
cnt  = cnt(cnt.Flag3 == false,:);
cnt  = cnt(cnt.Flag4 == false,:);
numFrames =  length(unique(cnt.Frame));
numBins = 1/APResolution+1;
MeanFluoAP = zeros(numFrames, numBins);
StdFluoAP = zeros(numFrames, numBins);
NumNucAP = zeros(numFrames, numBins);
DiffMeanFluoAP = zeros(numFrames, numBins);
DiffStdFluoAP = zeros(numFrames, numBins);
profFolder = 'APprofiles';
timeFolder = 'Timeprofiles';
% First, get the name of the folder you're using.

% Finally, create the folder if it doesn't exist already.

APbins = 0:APResolution:1;
fc = 1;
AllFrames = unique(cnt.Frame);
AllTimes = [];
for i=1:length(AllFrames)
    AllTimes(i) =  mean(cnt.Time(cnt.Frame == AllFrames(i)));
end
for i = 1:length(AllFrames)
    f = AllFrames(i);
    t = AllTimes(i);
    close all
    for b = 0:(numBins-1)
        APbin = b*APResolution;
        sub_table = cnt((cnt.Frame == f) & (abs(cnt.APbin- APbin) < 0.0001), :);
        cutoff = quantile(sub_table.Fluo2, quantile_cutoff);
        if size(sub_table, 5) >= 10
            cutoff = quantile(sub_table.Fluo2, quantile_cutoff);
            fluoVector = sub_table.Fluo2(sub_table.Fluo2 >= cutoff);
        else
            fluoVector = sub_table.Fluo2;
        end
        MeanFluo = mean(fluoVector);
        StdFluo = std(fluoVector);
        MeanFluoAP(i, b+1) = MeanFluo;
        StdFluoAP(i, b+1) = StdFluo;
        NumNucAP(i, b+1) = size(fluoVector, 1);
    end
%     fig = figure(fc);
%     ax = axes('Parent', fig);
    use_idx = NumNucAP(i,:) >= min_nuc_count;
    if all(isnan(MeanFluoAP(i,use_idx)))
        continue
    end
    minFluoBin= find(MeanFluoAP(i,:) == min(MeanFluoAP(i,use_idx)));

    DiffMeanFluoAP(i,:) = MeanFluoAP(i,:)- MeanFluoAP(i,minFluoBin(1));
    DiffStdFluoAP(i,:) = sqrt(StdFluoAP(i,:).^2+StdFluoAP(i,minFluoBin( 1)).^2);
%     eb = errorbar(APbins(use_idx), DiffMeanFluoAP(i, use_idx), DiffStdFluoAP(i, use_idx)./NumNucAP(i, use_idx), 'vertical', 'LineStyle', 'none');
%     hold on
%     set(eb, 'color', colors, 'LineWidth', 1);
%     set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
%     plot(APbins(use_idx), DiffMeanFluoAP(i,use_idx), '.-');
% 
%     %histogram(sub_table.Fluo)
%     xlabel('Fraction Embryo Length')
%     ylabel('Fluo (AU)')
%     ylim([0, (max(cnt.Fluo)-min(cnt.Fluo))*1.1])
%     xlim([0, 0.95])
%     %CurrentFrame = sub_table.Frame(1);
% 
% %     ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
% %             FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
% %     imshow(ImageHis,'DisplayRange',[]);
%     title(ax, {Prefix_label,...
%          ['Frame: ',num2str(f), ', Time: ',num2str(t)]})
%     hold on 


%     scatter(sub14.xPos, sub14.yPos,40, sub14.Fluo, 'filled')
%     plot(sub14.MedianDV, sub14.Fluo,'o')


    %ylabel('Fluo (AU)')
    %xlabel('Median DV (pixels)')
%     saveas(fig,[DropboxFolder,filesep,Prefix,filesep, profFolder, filesep, 'FluoProfile_frame', num2str(f),'.png']);
%     fc = fc+1;
%     close all
end

close all
fc = 1;

for b = 0:(numBins-1)
    use_idx = NumNucAP(:,b+1) >= min_nuc_count;
    if all(isnan(DiffMeanFluoAP(:, b+1)))
        continue 
    end
    close all
    fig = figure(fc);
    ax = axes('Parent', fig);
    for j=1:length(nc_frames)
        nc_start_frame = nc_frames(j);
        xline(nc_start_frame, 'r')
        hold on
    end
    
    eb = errorbar(AllFrames(use_idx), MeanFluoAP(use_idx,b+1), StdFluoAP(use_idx,b+1)./NumNucAP(use_idx,b+1), 'vertical', 'LineStyle', 'none');
    hold on
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(AllFrames(use_idx), MeanFluoAP(use_idx,b+1), '.-', 'color', colors);

    %histogram(sub14.Fluo)
    xlabel('Frame')
    ylabel('Fluo (AU)')
    ylim([0, (max(cnt.Fluo))*1.1])
    %xlim([0.05, 0.95])
    %CurrentFrame = sub14.Frame(1);

%     ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
%             FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
%     imshow(ImageHis,'DisplayRange',[]);
    title(ax, {Prefix_label,...
        ['Fraction Embryo Length: ',num2str(b*APResolution)]})
%     hold on 
    if ~all(isnan(MeanFluoAP(use_idx, b+1)))
        saveas(fig,[DropboxFolder,filesep,Prefix,filesep, timeFolder, filesep,'FrameProfile_AP', num2str(b),'.png']);
        fc = fc+1;
    end
    
   
end

close all

save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAP_Frames.mat'],...
        'MeanFluoAP');

save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAP_Frames.mat'],...
        'StdFluoAP');

save([DropboxFolder,filesep,Prefix,filesep,'NumNucAP_Frames.mat'],...
        'NumNucAP');
    
save([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAP_Frames.mat'],...
    'DiffMeanFluoAP');
save([DropboxFolder,filesep,Prefix,filesep,'DiffStdFluoAP_Frames.mat'],...
    'DiffStdFluoAP');





%%
close all
AllTimes = unique(round(cnt.Time, 1));
cnt.RoundedTime = round(cnt.Time, 1);
numTimes = length(AllTimes);
numBins = 1/APResolution+1;
MeanFluoAP = zeros(numTimes, numBins);
StdFluoAP = zeros(numTimes, numBins);
NumNucAP = zeros(numTimes, numBins);
DiffMeanFluoAP = zeros(numTimes, numBins);
DiffStdFluoAP = zeros(numTimes, numBins);

APbins = 0:APResolution:1;
fc = 1;

for i = 1:length(AllTimes)
    t = AllTimes(i);
    close all
    for b = 0:(numBins-1)
        APbin = b*APResolution;
        sub_table = cnt((cnt.RoundedTime == t) & (abs(cnt.APbin- APbin) < 0.0001), :);
        cutoff = quantile(sub_table.Fluo2, quantile_cutoff);
        if size(sub_table, 5) >= 10
            cutoff = quantile(sub_table.Fluo2, quantile_cutoff);
            fluoVector = sub_table.Fluo2(sub_table.Fluo2 >= cutoff);
        else
            fluoVector = sub_table.Fluo2;
        end
        MeanFluo = mean(fluoVector);
        StdFluo = std(fluoVector);
        MeanFluoAP(i, b+1) = MeanFluo;
        StdFluoAP(i, b+1) = StdFluo;
        NumNucAP(i, b+1) = size(fluoVector, 1);
    end

end


%NCTimes = 
for i = 1:numTimes
    t = AllTimes(i);
    close all
    for b = 0:(numBins-1)
        APbin = b*APResolution;
        sub_table = cnt((cnt.RoundedTime == t) & (abs(cnt.APbin- APbin) < 0.0001), :);
        if size(sub_table, 5) >= 10
            cutoff = quantile(sub_table.Fluo2, quantile_cutoff);
            fluoVector = sub_table.Fluo2(sub_table.Fluo2 >= cutoff);
        else
            fluoVector = sub_table.Fluo2;
        end
        MeanFluo = mean(fluoVector);
        StdFluo = std(fluoVector);
        MeanFluoAP(i, b+1) = MeanFluo;
        StdFluoAP(i, b+1) = StdFluo;
        NumNucAP(i, b+1) = size(fluoVector, 1);
    end
    use_idx = NumNucAP(i,:) >= min_nuc_count;
    if all(isnan(MeanFluoAP(i,use_idx)))
        continue
    end
    minFluoBin= find(MeanFluoAP(i,:) == min(MeanFluoAP(i,use_idx)));

    DiffMeanFluoAP(i,:) = MeanFluoAP(i,:)- MeanFluoAP(i,minFluoBin(1));
    DiffStdFluoAP(i,:) = sqrt(StdFluoAP(i,:).^2+StdFluoAP(i,minFluoBin( 1)).^2);

end

fc = 1;
for b = 0:(numBins-1)
    use_idx = NumNucAP(:,b+1) >= min_nuc_count;
    if all(isnan(DiffMeanFluoAP(:, b+1)))
        continue 
    end
    close all
    fig = figure(fc);
    ax = axes('Parent', fig);
    for j=1:length(nc_frames)
        nc_start_frame = APDivision(uint16(nuclear_cycles(j)),b+1)+1;
        nc_start_time = FrameInfo(nc_start_frame).Time/60;
        xline(nc_start_time, 'r')
        hold on
    end
    eb = errorbar(AllTimes(use_idx), MeanFluoAP(use_idx,b+1), StdFluoAP(use_idx,b+1)./NumNucAP(use_idx,b+1), 'vertical', 'LineStyle', 'none');
    hold on
    set(eb, 'color', colors, 'LineWidth', 1);
    set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(AllTimes(use_idx), MeanFluoAP(use_idx,b+1), '.-', 'color', colors);

    %histogram(sub14.Fluo)
    xlabel('Time (min)')
    ylabel('Fluo (AU)')
    ylim([0, (max(cnt.Fluo))*1.1])
    xlim([0.0, max(AllTimes)+1])
    %CurrentFrame = sub14.Frame(1);

%     ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
%             FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
%     imshow(ImageHis,'DisplayRange',[]);
    title(ax, {Prefix_label,...
        ['Fraction Embryo Length: ',num2str(b*APResolution)]})
%     hold on 
    if ~all(isnan(MeanFluoAP(use_idx, b+1)))
        saveas(fig,[DropboxFolder,filesep,Prefix,filesep, timeFolder, filesep,'TimeProfile_AP', num2str(t),'min.png']);
        fc = fc+1;
    end
    
   
end

close all

save([DropboxFolder,filesep,Prefix,filesep,'MeanFluoAP_Times.mat'],...
        'MeanFluoAP');

save([DropboxFolder,filesep,Prefix,filesep,'StdFluoAP_Times.mat'],...
        'StdFluoAP');

save([DropboxFolder,filesep,Prefix,filesep,'NumNucAPNC_Times.mat'],...
        'NumNucAP');
save([DropboxFolder,filesep,Prefix,filesep,'DiffMeanFluoAP_Frames.mat'],...
    'DiffMeanFluoAP');
save([DropboxFolder,filesep,Prefix,filesep,'DiffStdFluoAP_Frames.mat'],...
    'DiffStdFluoAP');
save([DropboxFolder,filesep,Prefix,filesep,'Times.mat'],...
    'AllTimes');

