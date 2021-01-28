function PlotMS2FluoProfilesAP(Prefix, varargin)
% PlotNuclearFluoProfilesAP.m
% author: Gabriella Martini
% date created: 12/30/20
% date last modified: 1/6/20

nuclear_cycles = [10, 11, 12, 13, 14];
min_nuc_count = 5;
SkipProfiles = false;
SkipTraces = false;

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
            case{'PlotTitle'}
                PlotTitle = varargin{x+1};
                x = x+1;
            case{'PlottingColors'}
                PlottingColors = varargin{x+1};
                x = x+1;
            case{'SkipProfiles'}
                SkipProfiles = true;
            case{'SkipTraces'}
                SkipTraces = true;    
            otherwise
                error(['Flag "', varargin{x},'" not valid'])
        end
        x = x +1;
    end
end


if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif  ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default", "pboc", or "gradient".') % change to error
end


if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
end
%%

liveExperiment = LiveExperiment(Prefix);
APbins = 0:liveExperiment.APResolution:1;
FrameInfo = getFrameInfo(liveExperiment);
FrameTimes = [FrameInfo(:).Time]; % in seconds
numFrames = length(FrameTimes);
nc_info = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, length(FrameInfo)];
load([liveExperiment.resultsFolder, 'MeanProfiles.mat'])
DataFolder=liveExperiment.resultsFolder;
Prefix_label =  strrep(Prefix, '_5C', '.5C');

profSubFolder = 'APprofiles';
timeSubFolder = 'Timeprofiles';
% First, get the name of the folder you're using.
% For example if your folder is 'D:\photos\Info', parentFolder  would = 'D:\photos, and deepestFolder would = 'Info'.
[parentFolder deepestFolder] = fileparts(DataFolder);
% Next, create a name for a subfolder within that.
% For example 'D:\photos\Info\DATA-Info'
profdir = sprintf('%s/%s/%s', parentFolder, deepestFolder, profSubFolder);
% Finally, create the folder if it doesn't exist already.
if ~exist(profdir, 'dir')
    mkdir(profdir);
end
timedir = sprintf('%s/%s/%s', parentFolder, deepestFolder, timeSubFolder);
% Finally, create the folder if it doesn't exist already.
if ~exist(timedir, 'dir')
    mkdir(timedir);
end

%%
if ~SkipProfiles
for nc=nuclear_cycles
    NCNonNaNIndices = find(sum(~isnan(squeeze(UnalignedCycleMeanTraces(:,:,nc-8))), 2)>0);

    if isempty(NCNonNaNIndices)
        continue
    end
    profdir2 = [profdir, filesep, 'nc', num2str(nc)];
    if ~exist(profdir2, 'dir')
        mkdir(profdir2)
    end
    profdir3 = [profdir2, filesep, datestr(now, 'yyyymmdd')];
    if ~exist(profdir3, 'dir')
        mkdir(profdir3)
    end
    
    
    
    NCTimes = UnalignedCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(NCTimes);
    NCcycleMeanTraces = UnalignedCycleMeanTraces(IncludedRows,:,:);
    NCcycleStdErrors = UnalignedCycleTraceStdErrors(IncludedRows,:,:);
    NCcycleNumNuclei = UnalignedCycleNumOnNuclei(IncludedRows,:,:);
    UnalignedNCNonNaNIndices = find(sum(~isnan(squeeze(NCcycleMeanTraces(:,:,nc-8))), 2)>0);
    NCTimes = NCTimes(UnalignedNCNonNaNIndices);
    NCcycleMeanTraces = squeeze(NCcycleMeanTraces(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycleStdErrors = squeeze(NCcycleStdErrors(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycleNumNuclei = squeeze(NCcycleNumNuclei(UnalignedNCNonNaNIndices,:,nc-8));

    
    NCTimes3D = Unaligned3DCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(NCTimes3D);
    NCcycle3DMeanTraces = Unaligned3DCycleMeanTraces(IncludedRows,:,:);
    NCcycle3DStdErrors = Unaligned3DCycleTraceStdErrors(IncludedRows,:,:);
    NCcycle3DNumNuclei = Unaligned3DCycleNumOnNuclei(IncludedRows,:,:);
    UnalignedNCNonNaNIndices = find(sum(~isnan(squeeze(NCcycle3DMeanTraces(:,:,nc-8))), 2)>0);
    NCTimes3D = NCTimes3D(UnalignedNCNonNaNIndices);
    NCcycle3DMeanTraces = squeeze(NCcycle3DMeanTraces(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycle3DStdErrors = squeeze(NCcycle3DStdErrors(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycle3DNumNuclei = squeeze(NCcycle3DNumNuclei(UnalignedNCNonNaNIndices,:,nc-8));
    
    
    AlignedNCTimes = AnaphaseAlignedCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(AlignedNCTimes);
    AlignedNCcycleMeanTraces = AnaphaseAlignedCycleMeanTraces(IncludedRows,:,:);
    AlignedNCcycleStdErrors = AnaphaseAlignedCycleTraceStdErrors(IncludedRows,:,:);
    AlignedNCcycleNumNuclei = AnaphaseAlignedCycleNumOnNuclei(IncludedRows,:,:);
    AlignedNCNonNaNIndices = find(sum(~isnan(squeeze(AlignedNCcycleMeanTraces(:,:,nc-8))), 2)>0);
    AlignedNCTimes = AlignedNCTimes(AlignedNCNonNaNIndices);
    AlignedNCcycleMeanTraces = squeeze(AlignedNCcycleMeanTraces(AlignedNCNonNaNIndices,:,nc-8));
    AlignedNCcycleStdErrors = squeeze(AlignedNCcycleStdErrors(AlignedNCNonNaNIndices,:,nc-8));
    AlignedNCcycleNumNuclei = squeeze(AlignedNCcycleNumNuclei(AlignedNCNonNaNIndices,:,nc-8));
    
    
    AlignedNC3DTimes = AnaphaseAligned3DCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(AlignedNC3DTimes);
    AlignedNCcycle3DMeanTraces = AnaphaseAligned3DCycleMeanTraces(IncludedRows,:,:);
    AlignedNCcycle3DStdErrors = AnaphaseAligned3DCycleTraceStdErrors(IncludedRows,:,:);
    AlignedNCcycle3DNumNuclei = AnaphaseAligned3DCycleNumOnNuclei(IncludedRows,:,:);
    AlignedNC3DNonNaNIndices = find(sum(~isnan(squeeze(AlignedNCcycle3DMeanTraces(:,:,nc-8))), 2)>0);
    AlignedNC3DTimes = AlignedNC3DTimes(AlignedNC3DNonNaNIndices);
    AlignedNCcycle3DMeanTraces = squeeze(AlignedNCcycle3DMeanTraces(AlignedNC3DNonNaNIndices,:,nc-8));
    AlignedNCcycle3DStdErrors = squeeze(AlignedNCcycle3DStdErrors(AlignedNC3DNonNaNIndices,:,nc-8));
    AlignedNCcycle3DNumNuclei = squeeze(AlignedNCcycle3DNumNuclei(AlignedNC3DNonNaNIndices,:,nc-8));
    
    
    TbinnedNCTimes = TbinnedCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(TbinnedNCTimes);
    TbinnedNCcycleMeanTraces = TbinnedCycleMeanTraces(IncludedRows,:,:);
    TbinnedNCcycleStdErrors = TbinnedCycleTraceStdErrors(IncludedRows,:,:);
    TbinnedNCcycleNumNuclei = TbinnedCycleNumOnNuclei(IncludedRows,:,:);
    TbinnedNCNonNaNIndices = find(sum(~isnan(squeeze(TbinnedNCcycleMeanTraces(:,:,nc-8))), 2)>0);
    TbinnedNCTimes = TbinnedNCTimes(TbinnedNCNonNaNIndices);
    TbinnedNCcycleMeanTraces = squeeze(TbinnedNCcycleMeanTraces(TbinnedNCNonNaNIndices,:,nc-8));
    TbinnedNCcycleStdErrors = squeeze(TbinnedNCcycleStdErrors(TbinnedNCNonNaNIndices,:,nc-8));
    TbinnedNCcycleNumNuclei = squeeze(TbinnedNCcycleNumNuclei(TbinnedNCNonNaNIndices,:,nc-8));
    
    TbinnedNC3DTimes = Tbinned3DCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(TbinnedNC3DTimes);
    TbinnedNCcycle3DMeanTraces = Tbinned3DCycleMeanTraces(IncludedRows,:,:);
    TbinnedNCcycle3DStdErrors = Tbinned3DCycleTraceStdErrors(IncludedRows,:,:);
    TbinnedNCcycle3DNumNuclei = Tbinned3DCycleNumOnNuclei(IncludedRows,:,:);
    TbinnedNC3DNonNaNIndices = find(sum(~isnan(squeeze(TbinnedNCcycle3DMeanTraces(:,:,nc-8))), 2)>0);
    TbinnedNC3DTimes = TbinnedNC3DTimes(TbinnedNC3DNonNaNIndices);
    TbinnedNCcycle3DMeanTraces = squeeze(TbinnedNCcycle3DMeanTraces(TbinnedNC3DNonNaNIndices,:,nc-8));
    TbinnedNCcycle3DStdErrors = squeeze(TbinnedNCcycle3DStdErrors(TbinnedNC3DNonNaNIndices,:,nc-8));
    TbinnedNCcycle3DNumNuclei = squeeze(TbinnedNCcycle3DNumNuclei(TbinnedNC3DNonNaNIndices,:,nc-8));
    %%
    close all
    
    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);
    
    eb{1} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on
    set(eb{1}, 'color', colors(1,:), 'LineWidth', 1);
    set(get(get(eb{1}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(1,:));
    eb{2} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    set(eb{2}, 'color', colors(2,:), 'LineWidth', 1);
    set(get(get(eb{2}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(2,:));
    hold off
    xlabel('Fraction Embryo Length')
    ylabel('Fluo (AU)')
    ylim([0, max([max(max(NCcycleMeanTraces+NCcycleStdErrors)),max(max(NCcycle3DMeanTraces+NCcycle3DStdErrors)),1])*1.2])
    xlim([0, 1])
    title(FrameProfAx, {Prefix_label,...
        ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(0), ', NC Time: ',num2str(-1), ' min' ]})
    
    for i = 1:size(NCcycleMeanTraces, 1)
        
        use_idx = NCcycleNumNuclei(i,:) >= min_nuc_count;
        if sum(use_idx) == 0
            FrameProfAx.Children(3).XData = APbins;
            FrameProfAx.Children(3).YData = zeros(1, length(APbins));
            FrameProfAx.Children(4).XData = APbins;
            FrameProfAx.Children(4).YData = zeros(1, length(APbins));
            FrameProfAx.Children(4).YPositiveDelta = .1*ones(1, length(APbins));
            FrameProfAx.Children(4).YNegativeDelta = .1*ones(1, length(APbins));
            set(FrameProfAx.Children(3), 'Visible', 'off')
            set(FrameProfAx.Children(4), 'Visible', 'off')
        else
            FrameProfAx.Children(3).XData = APbins(use_idx);
            FrameProfAx.Children(3).YData = NCcycleMeanTraces(i, use_idx);
            FrameProfAx.Children(4).YData = NCcycleMeanTraces(i, use_idx);
            FrameProfAx.Children(4).XData = APbins(use_idx);
            FrameProfAx.Children(4).YPositiveDelta = NCcycleStdErrors(i, use_idx);
            FrameProfAx.Children(4).YNegativeDelta  = NCcycleStdErrors(i, use_idx);
            set(FrameProfAx.Children(3), 'Visible', 'on')
            set(FrameProfAx.Children(4), 'Visible', 'on')
            
        end
        
        use_idx2 = NCcycle3DNumNuclei(i,:) >= min_nuc_count;
        if sum(use_idx2) == 0
            FrameProfAx.Children(1).XData = APbins;
            FrameProfAx.Children(1).YData = zeros(1, length(APbins));
            FrameProfAx.Children(2).XData = APbins;
            FrameProfAx.Children(2).YData = zeros(1, length(APbins));
            FrameProfAx.Children(2).YPositiveDelta = .1*ones(1, length(APbins));
            FrameProfAx.Children(2).YNegativeDelta = .1*ones(1, length(APbins));
            set(FrameProfAx.Children(1), 'Visible', 'off')
            set(FrameProfAx.Children(2), 'Visible', 'off')
        else
            FrameProfAx.Children(1).XData = APbins(use_idx2);
            FrameProfAx.Children(1).YData = NCcycle3DMeanTraces(i, use_idx2);
            FrameProfAx.Children(2).YData = NCcycle3DMeanTraces(i, use_idx2);
            FrameProfAx.Children(2).XData = APbins(use_idx2);
            FrameProfAx.Children(2).YPositiveDelta = NCcycle3DStdErrors(i, use_idx2);
            FrameProfAx.Children(2).YNegativeDelta  = NCcycle3DStdErrors(i, use_idx2);
            set(FrameProfAx.Children(1), 'Visible', 'on')
            set(FrameProfAx.Children(2), 'Visible', 'on')
            
        end
        if ~sum(use_idx2) & ~sum(use_idx)
            continue
        end
        
        
        try
            title(FrameProfAx, {Prefix_label,...
                ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(i), ', NC Time: ',num2str(round(NCTimes(i), 2)), ' min' ]})
        end
        hlegend = legend({'Fluo', 'Fluo3DGauss'}, 'Location', 'northeast',...
            'FontSize', 10);
        %         if (sum(use_idx) == 0) & (sum(use_idx2) == 0)
        %             continue
        %         end
        saveas(FrameProfFig,[profdir3,filesep,  'FluoProfile_frame', num2str(i),'.png']);
        
    end
    
    close all
    
    FrameProfFig = figure(1);
    FrameProfAx = axes(FrameProfFig);
    
    eb{1} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    hold on
    set(eb{1}, 'color', colors(1,:), 'LineWidth', 1);
    set(get(get(eb{1}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(1,:));
    eb{2} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    set(eb{2}, 'color', colors(2,:), 'LineWidth', 1);
    set(get(get(eb{2}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(2,:));
    eb{3} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    set(eb{3}, 'color', colors(3,:), 'LineWidth', 1);
    set(get(get(eb{3}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(3,:));
    eb{4} = errorbar(APbins, ones(1, length(APbins)), .1*ones(1, length(APbins)), 'vertical', 'LineStyle', 'none');
    set(eb{4}, 'color', colors(4,:), 'LineWidth', 1);
    set(get(get(eb{4}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(APbins, ones(1, length(APbins)), '.-', 'Color', colors(4,:));
    hold off
    xlabel('Fraction Embryo Length')
    ylabel('Fluo (AU)')
    ylim([0, max([max(max(AlignedNCcycleMeanTraces+AlignedNCcycleStdErrors)),max(max(AlignedNCcycle3DMeanTraces+AlignedNCcycle3DStdErrors)),...
        max(max(TbinnedNCcycleMeanTraces+TbinnedNCcycleStdErrors)),max(max(TbinnedNCcycle3DMeanTraces+TbinnedNCcycle3DStdErrors)), 1])*1.2])
    xlim([0, 1])
    title(FrameProfAx, {Prefix_label,...
        ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(0), ', NC Time: ',num2str(-1), ' min' ]})
    deltaT = diff(AlignedNCTimes);
    deltaT = deltaT(1);
    SharedMaxLength = max([length(AlignedNCTimes), length(AlignedNC3DTimes), length(TbinnedNCTimes), length(TbinnedNC3DTimes)]);
    if length(AlignedNCTimes) == SharedMaxLength
        SharedTimes = AlignedNCTimes;
    elseif length(AlignedNC3DTimes) == SharedMaxLength
        SharedTimes = AlignedNC3DTimes;
    elseif length(TbinnedNCTimes) == SharedMaxLength
        SharedTimes = TbinnedNCTimes;
    elseif length(TbinnedNC3DTimes) == SharedMaxLength
        SharedTimes = TbinnedNC3DTimes;
    end
    
    for i = 1:SharedMaxLength
        saveplot = 0;
        if length(AlignedNCTimes) < i
            FrameProfAx.Children(3).XData = APbins;
            FrameProfAx.Children(3).YData = zeros(1, length(APbins));
            FrameProfAx.Children(4).XData = APbins;
            FrameProfAx.Children(4).YData = zeros(1, length(APbins));
            FrameProfAx.Children(4).YPositiveDelta = .1*ones(1, length(APbins));
            FrameProfAx.Children(4).YNegativeDelta = .1*ones(1, length(APbins));
            set(FrameProfAx.Children(3), 'Visible', 'off')
            set(FrameProfAx.Children(4), 'Visible', 'off')
        else
            use_idx = AlignedNCcycleNumNuclei(i,:) >= min_nuc_count;
            if sum(use_idx) == 0
                FrameProfAx.Children(3).XData = APbins;
                FrameProfAx.Children(3).YData = zeros(1, length(APbins));
                FrameProfAx.Children(4).XData = APbins;
                FrameProfAx.Children(4).YData = zeros(1, length(APbins));
                FrameProfAx.Children(4).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(4).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(3), 'Visible', 'off')
                set(FrameProfAx.Children(4), 'Visible', 'off')
            else
                saveplot = 1;
                FrameProfAx.Children(3).XData = APbins(use_idx);
                FrameProfAx.Children(3).YData = AlignedNCcycleMeanTraces(i, use_idx);
                FrameProfAx.Children(4).YData = AlignedNCcycleMeanTraces(i, use_idx);
                FrameProfAx.Children(4).XData = APbins(use_idx);
                FrameProfAx.Children(4).YPositiveDelta = AlignedNCcycleStdErrors(i, use_idx);
                FrameProfAx.Children(4).YNegativeDelta  = AlignedNCcycleStdErrors(i, use_idx);
                set(FrameProfAx.Children(3), 'Visible', 'on')
                set(FrameProfAx.Children(4), 'Visible', 'on')
            end
        end
        if length(AlignedNC3DTimes) < i
            FrameProfAx.Children(1).XData = APbins;
            FrameProfAx.Children(1).YData = zeros(1, length(APbins));
            FrameProfAx.Children(2).XData = APbins;
            FrameProfAx.Children(2).YData = zeros(1, length(APbins));
            FrameProfAx.Children(2).YPositiveDelta = .1*ones(1, length(APbins));
            FrameProfAx.Children(2).YNegativeDelta = .1*ones(1, length(APbins));
            set(FrameProfAx.Children(1), 'Visible', 'off')
            set(FrameProfAx.Children(2), 'Visible', 'off')
        else
            use_idx2 = AlignedNCcycle3DNumNuclei(i,:) >= min_nuc_count;
            if sum(use_idx2) == 0
                FrameProfAx.Children(1).XData = APbins;
                FrameProfAx.Children(1).YData = zeros(1, length(APbins));
                FrameProfAx.Children(2).XData = APbins;
                FrameProfAx.Children(2).YData = zeros(1, length(APbins));
                FrameProfAx.Children(2).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(2).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(1), 'Visible', 'off')
                set(FrameProfAx.Children(2), 'Visible', 'off')
            else
                saveplot = 1;
                FrameProfAx.Children(1).XData = APbins(use_idx2);
                FrameProfAx.Children(1).YData = AlignedNCcycle3DMeanTraces(i, use_idx2);
                FrameProfAx.Children(2).YData = AlignedNCcycle3DMeanTraces(i, use_idx2);
                FrameProfAx.Children(2).XData = APbins(use_idx2);
                FrameProfAx.Children(2).YPositiveDelta = AlignedNCcycle3DStdErrors(i, use_idx2);
                FrameProfAx.Children(2).YNegativeDelta  = AlignedNCcycle3DStdErrors(i, use_idx2);
                set(FrameProfAx.Children(1), 'Visible', 'on')
                set(FrameProfAx.Children(2), 'Visible', 'on')
            end
        end
        
        if length(TbinnedNCTimes) < i
            FrameProfAx.Children(7).XData = APbins;
            FrameProfAx.Children(7).YData = zeros(1, length(APbins));
            FrameProfAx.Children(8).XData = APbins;
            FrameProfAx.Children(8).YData = zeros(1, length(APbins));
            FrameProfAx.Children(8).YPositiveDelta = .1*ones(1, length(APbins));
            FrameProfAx.Children(8).YNegativeDelta = .1*ones(1, length(APbins));
            set(FrameProfAx.Children(7), 'Visible', 'off')
            set(FrameProfAx.Children(8), 'Visible', 'off')
        else
            use_idx3 = TbinnedNCcycleNumNuclei(i,:) >= min_nuc_count;
            if sum(use_idx3) == 0
                FrameProfAx.Children(7).XData = APbins;
                FrameProfAx.Children(7).YData = zeros(1, length(APbins));
                FrameProfAx.Children(8).XData = APbins;
                FrameProfAx.Children(8).YData = zeros(1, length(APbins));
                FrameProfAx.Children(8).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(8).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(7), 'Visible', 'off')
                set(FrameProfAx.Children(8), 'Visible', 'off')
            else
                saveplot = 1;
                FrameProfAx.Children(7).XData = APbins(use_idx3);
                FrameProfAx.Children(7).YData = TbinnedNCcycleMeanTraces(i, use_idx3);
                FrameProfAx.Children(8).YData = TbinnedNCcycleMeanTraces(i, use_idx3);
                FrameProfAx.Children(8).XData = APbins(use_idx3);
                FrameProfAx.Children(8).YPositiveDelta = TbinnedNCcycleStdErrors(i, use_idx3);
                FrameProfAx.Children(8).YNegativeDelta  = TbinnedNCcycleStdErrors(i, use_idx3);
                set(FrameProfAx.Children(7), 'Visible', 'on')
                set(FrameProfAx.Children(8), 'Visible', 'on')
            end
        end
        if length(TbinnedNC3DTimes) < i
            FrameProfAx.Children(5).XData = APbins;
            FrameProfAx.Children(5).YData = zeros(1, length(APbins));
            FrameProfAx.Children(6).XData = APbins;
            FrameProfAx.Children(6).YData = zeros(1, length(APbins));
            FrameProfAx.Children(6).YPositiveDelta = .1*ones(1, length(APbins));
            FrameProfAx.Children(6).YNegativeDelta = .1*ones(1, length(APbins));
            set(FrameProfAx.Children(5), 'Visible', 'off')
            set(FrameProfAx.Children(6), 'Visible', 'off')
        else
            use_idx4 = TbinnedNCcycle3DNumNuclei(i,:) >= min_nuc_count;
            if sum(use_idx4) == 0
                FrameProfAx.Children(5).XData = APbins;
                FrameProfAx.Children(5).YData = zeros(1, length(APbins));
                FrameProfAx.Children(6).XData = APbins;
                FrameProfAx.Children(6).YData = zeros(1, length(APbins));
                FrameProfAx.Children(6).YPositiveDelta = .1*ones(1, length(APbins));
                FrameProfAx.Children(6).YNegativeDelta = .1*ones(1, length(APbins));
                set(FrameProfAx.Children(5), 'Visible', 'off')
                set(FrameProfAx.Children(6), 'Visible', 'off')
            else
                saveplot = 1;
                FrameProfAx.Children(5).XData = APbins(use_idx4);
                FrameProfAx.Children(5).YData = TbinnedNCcycle3DMeanTraces(i, use_idx4);
                FrameProfAx.Children(6).YData = TbinnedNCcycle3DMeanTraces(i, use_idx4);
                FrameProfAx.Children(6).XData = APbins(use_idx4);
                FrameProfAx.Children(6).YPositiveDelta = TbinnedNCcycle3DStdErrors(i, use_idx4);
                FrameProfAx.Children(6).YNegativeDelta  = TbinnedNCcycle3DStdErrors(i, use_idx4);
                set(FrameProfAx.Children(5), 'Visible', 'on')
                set(FrameProfAx.Children(6), 'Visible', 'on')
            end
        end
        if ~saveplot
            continue
        end
        
        
        
        try
            title(FrameProfAx, {Prefix_label,...
                ['Nuclear Cycle ', num2str(nc), ', Frame: ',num2str(i), ', NC Time: ',num2str(round(SharedTimes(i), 2)), ' min' ]})
        end
        hlegend = legend({'Time binned Fluo', 'Time binned Fluo3DGauss', 'Anaphase Aligned Fluo', 'Anaphase Aligned Fluo3DGauss'}, 'Location', 'northeast',...
            'FontSize', 10);
        %         if (sum(use_idx) == 0) & (sum(use_idx2) == 0)
        %             continue
        %         end
        saveas(FrameProfFig,[profdir3,filesep,  'TBinnedFluoProfile_frame', num2str(i),'.png']);
        
    end
end
end
%%
if ~SkipTraces
for nc=nuclear_cycles
    NCNonNaNIndices = find(sum(~isnan(squeeze(UnalignedCycleMeanTraces(:,:,nc-8))), 2)>0);
    if isempty(NCNonNaNIndices)
        continue
    end
    timedir2 = [timedir, filesep, 'nc', num2str(nc)];
    if ~exist(timedir2, 'dir')
        mkdir(timedir2)
    end
    timedir3 = [timedir2, filesep, datestr(now, 'yyyymmdd')];
    if ~exist(timedir3, 'dir')
        mkdir(timedir3)
    end
     
   
    NCTimes = UnalignedCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(NCTimes);
    NCcycleMeanTraces = UnalignedCycleMeanTraces(IncludedRows,:,:);
    NCcycleStdErrors = UnalignedCycleTraceStdErrors(IncludedRows,:,:);
    NCcycleNumNuclei = UnalignedCycleNumOnNuclei(IncludedRows,:,:);
    UnalignedNCNonNaNIndices = find(sum(~isnan(squeeze(NCcycleMeanTraces(:,:,nc-8))), 2)>0);
    NCTimes = NCTimes(UnalignedNCNonNaNIndices);
    NCcycleMeanTraces = squeeze(NCcycleMeanTraces(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycleStdErrors = squeeze(NCcycleStdErrors(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycleNumNuclei = squeeze(NCcycleNumNuclei(UnalignedNCNonNaNIndices,:,nc-8));

    
    NCTimes3D = Unaligned3DCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(NCTimes3D);
    NCcycle3DMeanTraces = Unaligned3DCycleMeanTraces(IncludedRows,:,:);
    NCcycle3DStdErrors = Unaligned3DCycleTraceStdErrors(IncludedRows,:,:);
    NCcycle3DNumNuclei = Unaligned3DCycleNumOnNuclei(IncludedRows,:,:);
    UnalignedNCNonNaNIndices = find(sum(~isnan(squeeze(NCcycle3DMeanTraces(:,:,nc-8))), 2)>0);
    NCTimes3D = NCTimes3D(UnalignedNCNonNaNIndices);
    NCcycle3DMeanTraces = squeeze(NCcycle3DMeanTraces(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycle3DStdErrors = squeeze(NCcycle3DStdErrors(UnalignedNCNonNaNIndices,:,nc-8));
    NCcycle3DNumNuclei = squeeze(NCcycle3DNumNuclei(UnalignedNCNonNaNIndices,:,nc-8));
    
    
    AlignedNCTimes = AnaphaseAlignedCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(AlignedNCTimes);
    AlignedNCcycleMeanTraces = AnaphaseAlignedCycleMeanTraces(IncludedRows,:,:);
    AlignedNCcycleStdErrors = AnaphaseAlignedCycleTraceStdErrors(IncludedRows,:,:);
    AlignedNCcycleNumNuclei = AnaphaseAlignedCycleNumOnNuclei(IncludedRows,:,:);
    AlignedNCNonNaNIndices = find(sum(~isnan(squeeze(AlignedNCcycleMeanTraces(:,:,nc-8))), 2)>0);
    AlignedNCTimes = AlignedNCTimes(AlignedNCNonNaNIndices);
    AlignedNCcycleMeanTraces = squeeze(AlignedNCcycleMeanTraces(AlignedNCNonNaNIndices,:,nc-8));
    AlignedNCcycleStdErrors = squeeze(AlignedNCcycleStdErrors(AlignedNCNonNaNIndices,:,nc-8));
    AlignedNCcycleNumNuclei = squeeze(AlignedNCcycleNumNuclei(AlignedNCNonNaNIndices,:,nc-8));
    
    
    AlignedNC3DTimes = AnaphaseAligned3DCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(AlignedNC3DTimes);
    AlignedNCcycle3DMeanTraces = AnaphaseAligned3DCycleMeanTraces(IncludedRows,:,:);
    AlignedNCcycle3DStdErrors = AnaphaseAligned3DCycleTraceStdErrors(IncludedRows,:,:);
    AlignedNCcycle3DNumNuclei = AnaphaseAligned3DCycleNumOnNuclei(IncludedRows,:,:);
    AlignedNC3DNonNaNIndices = find(sum(~isnan(squeeze(AlignedNCcycle3DMeanTraces(:,:,nc-8))), 2)>0);
    AlignedNC3DTimes = AlignedNC3DTimes(AlignedNC3DNonNaNIndices);
    AlignedNCcycle3DMeanTraces = squeeze(AlignedNCcycle3DMeanTraces(AlignedNC3DNonNaNIndices,:,nc-8));
    AlignedNCcycle3DStdErrors = squeeze(AlignedNCcycle3DStdErrors(AlignedNC3DNonNaNIndices,:,nc-8));
    AlignedNCcycle3DNumNuclei = squeeze(AlignedNCcycle3DNumNuclei(AlignedNC3DNonNaNIndices,:,nc-8));
    
    
    TbinnedNCTimes = TbinnedCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(TbinnedNCTimes);
    TbinnedNCcycleMeanTraces = TbinnedCycleMeanTraces(IncludedRows,:,:);
    TbinnedNCcycleStdErrors = TbinnedCycleTraceStdErrors(IncludedRows,:,:);
    TbinnedNCcycleNumNuclei = TbinnedCycleNumOnNuclei(IncludedRows,:,:);
    TbinnedNCNonNaNIndices = find(sum(~isnan(squeeze(TbinnedNCcycleMeanTraces(:,:,nc-8))), 2)>0);
    TbinnedNCTimes = TbinnedNCTimes(TbinnedNCNonNaNIndices);
    TbinnedNCcycleMeanTraces = squeeze(TbinnedNCcycleMeanTraces(TbinnedNCNonNaNIndices,:,nc-8));
    TbinnedNCcycleStdErrors = squeeze(TbinnedNCcycleStdErrors(TbinnedNCNonNaNIndices,:,nc-8));
    TbinnedNCcycleNumNuclei = squeeze(TbinnedNCcycleNumNuclei(TbinnedNCNonNaNIndices,:,nc-8));
    
    TbinnedNC3DTimes = Tbinned3DCycleFrameTimes{nc-8}/60;
    IncludedRows = 1:length(TbinnedNC3DTimes);
    TbinnedNCcycle3DMeanTraces = Tbinned3DCycleMeanTraces(IncludedRows,:,:);
    TbinnedNCcycle3DStdErrors = Tbinned3DCycleTraceStdErrors(IncludedRows,:,:);
    TbinnedNCcycle3DNumNuclei = Tbinned3DCycleNumOnNuclei(IncludedRows,:,:);
    TbinnedNC3DNonNaNIndices = find(sum(~isnan(squeeze(TbinnedNCcycle3DMeanTraces(:,:,nc-8))), 2)>0);
    TbinnedNC3DTimes = TbinnedNC3DTimes(TbinnedNC3DNonNaNIndices);
    TbinnedNCcycle3DMeanTraces = squeeze(TbinnedNCcycle3DMeanTraces(TbinnedNC3DNonNaNIndices,:,nc-8));
    TbinnedNCcycle3DStdErrors = squeeze(TbinnedNCcycle3DStdErrors(TbinnedNC3DNonNaNIndices,:,nc-8));
    TbinnedNCcycle3DNumNuclei = squeeze(TbinnedNCcycle3DNumNuclei(TbinnedNC3DNonNaNIndices,:,nc-8));
    close all
    
    FrameTraceFig = figure(1);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.6]);
    FrameTraceAx = axes(FrameTraceFig);
    Frames = 1:size(NCcycleMeanTraces, 1);
    eb{1} = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    hold on
    set(eb{1}, 'color', colors(5,:), 'LineWidth', 1);
    set(get(get(eb{1}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors(5,:));
    eb{2} = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    set(eb{2}, 'color', colors(6,:), 'LineWidth', 1);
    set(get(get(eb{2}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors(6,:));
    eb{3} = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    set(eb{3}, 'color', colors(1,:), 'LineWidth', 1);
    set(get(get(eb{3}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors(1,:));
    eb{4} = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    set(eb{4}, 'color', colors(2,:), 'LineWidth', 1);
    set(get(get(eb{4}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors(2,:));
    eb{5} = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    set(eb{5}, 'color', colors(3,:), 'LineWidth', 1);
    set(get(get(eb{5}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors(3,:));
    eb{6} = errorbar(NCTimes, ones(1, length(NCTimes)), .1*ones(1, length(NCTimes)), 'vertical', 'LineStyle', 'none');
    set(eb{6}, 'color', colors(4,:), 'LineWidth', 1);
    set(get(get(eb{6}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    plot(NCTimes, ones(1, length(NCTimes)), '.-', 'Color', colors(4,:));
    hold off
    xlabel('Time (min)')
    ylabel('Fluo (AU)')
    ylim([0, max([max(max(NCcycleMeanTraces+NCcycleStdErrors)),max(max(NCcycle3DMeanTraces+NCcycle3DStdErrors)),...
        max(max(AlignedNCcycleMeanTraces+AlignedNCcycleStdErrors)),max(max(AlignedNCcycle3DMeanTraces+AlignedNCcycle3DStdErrors)),...
        max(max(TbinnedNCcycleMeanTraces+TbinnedNCcycleStdErrors)),max(max(TbinnedNCcycle3DMeanTraces+TbinnedNCcycle3DStdErrors)),1])*1.2])
    xlim([0, max([max(NCTimes), max(AlignedNCTimes), max(AlignedNC3DTimes),max(TbinnedNCTimes),  max(TbinnedNC3DTimes)])])
    title(FrameTraceAx, {Prefix_label,...
        ['Nuclear Cycle ', num2str(nc), ', AP bin: ',num2str(0)]})
    
    for i = 1:size(NCcycleMeanTraces,2)
        if ~sum(~isnan(NCcycleMeanTraces(:, i))) & ~sum(~isnan(AlignedNCcycleMeanTraces(:, i))) & ~sum(~isnan(TbinnedNCcycleMeanTraces(:, i)))
            continue
        end
        if sum(~isnan(NCcycleMeanTraces(:, i)))
            use_idx = NCcycleNumNuclei(:,i) >= min_nuc_count;
            if sum(use_idx) == 0
                FrameTraceAx.Children(3).XData = NCTimes;
                FrameTraceAx.Children(3).YData = zeros(1, length(NCTimes));
                FrameTraceAx.Children(4).XData = NCTimes;
                FrameTraceAx.Children(4).YData = zeros(1, length(NCTimes));
                FrameTraceAx.Children(4).YPositiveDelta = .1*ones(1, length(NCTimes));
                FrameTraceAx.Children(4).YNegativeDelta = .1*ones(1, length(NCTimes));
                set(FrameTraceAx.Children(3), 'Visible', 'off')
                set(FrameTraceAx.Children(4), 'Visible', 'off')
            else
                FrameTraceAx.Children(3).YData = NCcycleMeanTraces(use_idx,i);
                FrameTraceAx.Children(3).XData = NCTimes(use_idx);
                FrameTraceAx.Children(4).YData = NCcycleMeanTraces(use_idx,i);
                FrameTraceAx.Children(4).XData = NCTimes(use_idx);
                FrameTraceAx.Children(4).YPositiveDelta = NCcycleStdErrors(use_idx,i);
                FrameTraceAx.Children(4).YNegativeDelta  = NCcycleStdErrors(use_idx, i);
                set(FrameTraceAx.Children(3), 'Visible', 'on')
                set(FrameTraceAx.Children(4), 'Visible', 'on')
            end
            
            use_idx2 = NCcycle3DNumNuclei(:,i) >= min_nuc_count;
            if sum(use_idx2) == 0
                FrameTraceAx.Children(1).XData = NCTimes;
                FrameTraceAx.Children(1).YData = zeros(1, length(NCTimes));
                FrameTraceAx.Children(2).XData = NCTimes;
                FrameTraceAx.Children(2).YData = zeros(1, length(NCTimes));
                FrameTraceAx.Children(2).YPositiveDelta = .1*ones(1, length(NCTimes));
                FrameTraceAx.Children(2).YNegativeDelta = .1*ones(1, length(NCTimes));
                set(FrameTraceAx.Children(1), 'Visible', 'off')
                set(FrameTraceAx.Children(2), 'Visible', 'off')
            else
                FrameTraceAx.Children(1).YData = NCcycle3DMeanTraces(use_idx2,i);
                FrameTraceAx.Children(1).XData = NCTimes(use_idx2);
                FrameTraceAx.Children(2).YData = NCcycle3DMeanTraces(use_idx2,i);
                FrameTraceAx.Children(2).XData = NCTimes(use_idx2);
                FrameTraceAx.Children(2).YPositiveDelta = NCcycle3DStdErrors(use_idx2,i);
                FrameTraceAx.Children(2).YNegativeDelta  = NCcycle3DStdErrors(use_idx2, i);
                set(FrameTraceAx.Children(1), 'Visible', 'on')
                set(FrameTraceAx.Children(2), 'Visible', 'on')
            end
        end
        if i <= size(AlignedNCcycleMeanTraces, 2)
            if sum(~isnan(AlignedNCcycleMeanTraces(:, i)))
                use_idx = AlignedNCcycleNumNuclei(:,i) >= min_nuc_count;
                if sum(use_idx) == 0
                    FrameTraceAx.Children(7).XData = NCTimes;
                    FrameTraceAx.Children(7).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(8).XData = NCTimes;
                    FrameTraceAx.Children(8).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(8).YPositiveDelta = .1*ones(1, length(NCTimes));
                    FrameTraceAx.Children(8).YNegativeDelta = .1*ones(1, length(NCTimes));
                    set(FrameTraceAx.Children(7), 'Visible', 'off')
                    set(FrameTraceAx.Children(8), 'Visible', 'off')
                else
                    FrameTraceAx.Children(7).YData = AlignedNCcycleMeanTraces(use_idx,i);
                    FrameTraceAx.Children(7).XData = AlignedNCTimes(use_idx);
                    FrameTraceAx.Children(8).YData = AlignedNCcycleMeanTraces(use_idx,i);
                    FrameTraceAx.Children(8).XData = AlignedNCTimes(use_idx);
                    FrameTraceAx.Children(8).YPositiveDelta = AlignedNCcycleStdErrors(use_idx,i);
                    FrameTraceAx.Children(8).YNegativeDelta  = AlignedNCcycleStdErrors(use_idx, i);
                    set(FrameTraceAx.Children(7), 'Visible', 'on')
                    set(FrameTraceAx.Children(8), 'Visible', 'on')
                end
                
                use_idx2 = AlignedNCcycle3DNumNuclei(:,i) >= min_nuc_count;
                if sum(use_idx2) == 0
                    FrameTraceAx.Children(5).XData = NCTimes;
                    FrameTraceAx.Children(5).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(6).XData = NCTimes;
                    FrameTraceAx.Children(6).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(6).YPositiveDelta = .1*ones(1, length(NCTimes));
                    FrameTraceAx.Children(6).YNegativeDelta = .1*ones(1, length(NCTimes));
                    set(FrameTraceAx.Children(5), 'Visible', 'off')
                    set(FrameTraceAx.Children(6), 'Visible', 'off')
                else
                    FrameTraceAx.Children(5).YData = AlignedNCcycle3DMeanTraces(use_idx2,i);
                    FrameTraceAx.Children(5).XData = AlignedNC3DTimes(use_idx2);
                    FrameTraceAx.Children(6).YData = AlignedNCcycle3DMeanTraces(use_idx2,i);
                    FrameTraceAx.Children(6).XData = AlignedNC3DTimes(use_idx2);
                    FrameTraceAx.Children(6).YPositiveDelta = AlignedNCcycle3DStdErrors(use_idx2,i);
                    FrameTraceAx.Children(6).YNegativeDelta  = AlignedNCcycle3DStdErrors(use_idx2, i);
                    set(FrameTraceAx.Children(5), 'Visible', 'on')
                    set(FrameTraceAx.Children(6), 'Visible', 'on')
                end
            end
        end
        if i <= size(TbinnedNCcycleMeanTraces, 2)
            if sum(~isnan(TbinnedNCcycleMeanTraces(:, i)))
                use_idx = TbinnedNCcycleNumNuclei(:,i) >= min_nuc_count;
                if sum(use_idx) == 0
                    FrameTraceAx.Children(11).XData = NCTimes;
                    FrameTraceAx.Children(11).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(12).XData = NCTimes;
                    FrameTraceAx.Children(12).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(12).YPositiveDelta = .1*ones(1, length(NCTimes));
                    FrameTraceAx.Children(12).YNegativeDelta = .1*ones(1, length(NCTimes));
                    set(FrameTraceAx.Children(11), 'Visible', 'off')
                    set(FrameTraceAx.Children(12), 'Visible', 'off')
                else
                    FrameTraceAx.Children(11).YData = TbinnedNCcycleMeanTraces(use_idx,i);
                    FrameTraceAx.Children(11).XData = TbinnedNCTimes(use_idx);
                    FrameTraceAx.Children(12).YData = TbinnedNCcycleMeanTraces(use_idx,i);
                    FrameTraceAx.Children(12).XData = TbinnedNCTimes(use_idx);
                    FrameTraceAx.Children(12).YPositiveDelta = TbinnedNCcycleStdErrors(use_idx,i);
                    FrameTraceAx.Children(12).YNegativeDelta  = TbinnedNCcycleStdErrors(use_idx, i);
                    set(FrameTraceAx.Children(11), 'Visible', 'on')
                    set(FrameTraceAx.Children(12), 'Visible', 'on')
                end
                
                use_idx2 = TbinnedNCcycle3DNumNuclei(:,i) >= min_nuc_count;
                if sum(use_idx2) == 0
                    FrameTraceAx.Children(9).XData = NCTimes;
                    FrameTraceAx.Children(9).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(10).XData = NCTimes;
                    FrameTraceAx.Children(10).YData = zeros(1, length(NCTimes));
                    FrameTraceAx.Children(10).YPositiveDelta = .1*ones(1, length(NCTimes));
                    FrameTraceAx.Children(10).YNegativeDelta = .1*ones(1, length(NCTimes));
                    set(FrameTraceAx.Children(9), 'Visible', 'off')
                    set(FrameTraceAx.Children(10), 'Visible', 'off')
                else
                    FrameTraceAx.Children(9).YData = TbinnedNCcycle3DMeanTraces(use_idx2,i);
                    FrameTraceAx.Children(9).XData = TbinnedNC3DTimes(use_idx2);
                    FrameTraceAx.Children(10).YData =  TbinnedNCcycle3DMeanTraces(use_idx2,i);
                    FrameTraceAx.Children(10).XData = TbinnedNC3DTimes(use_idx2);
                    FrameTraceAx.Children(10).YPositiveDelta =  TbinnedNCcycle3DStdErrors(use_idx2,i);
                    FrameTraceAx.Children(10).YNegativeDelta  =  TbinnedNCcycle3DStdErrors(use_idx2, i);
                    set(FrameTraceAx.Children(9), 'Visible', 'on')
                    set(FrameTraceAx.Children(10), 'Visible', 'on')
                end
            end
        end
        
        
        
        title(FrameTraceAx, {Prefix_label,...
            ['Nuclear Cycle ', num2str(nc), ', Fraction Embryo Length: ',num2str((i-1)*liveExperiment.APResolution)]})
        if sum(use_idx) == 0
            continue
        end
        hlegend = legend({'Time binned Fluo', 'Time binned Fluo3DGauss', 'Anaphase Aligned Fluo',...
            'Anaphase Aligned Fluo3DGauss','Fluo', 'Fluo3DGauss'}, 'Location', 'eastoutside',...
            'FontSize', 10);
        %hlegend.NumColumns = 2;
        saveas(FrameTraceFig,[timedir3,filesep, 'FluoTrace_bin', num2str(i),'.png']);
        
    end
end

%%
FrameTimesMin = FrameTimes/60;
Ntimepoints = 0;
NCStartFrames = NaN(1,6);
NCEndFrames= NaN(1,6);
MeanFluoAP = NaN(size(UnalignedCycleMeanTraces, 1)+50, size(UnalignedCycleMeanTraces, 2),...
    size(UnalignedCycleMeanTraces, 3));
StdFluoAP = NaN(size(UnalignedCycleMeanTraces, 1)+50, size(UnalignedCycleMeanTraces, 2),...
    size(UnalignedCycleMeanTraces, 3));
NumNucAP = NaN(size(UnalignedCycleMeanTraces, 1)+50, size(UnalignedCycleMeanTraces, 2),...
    size(UnalignedCycleMeanTraces, 3));
Mean3DFluoAP = NaN(size(UnalignedCycleMeanTraces, 1)+50, size(UnalignedCycleMeanTraces, 2),...
    size(UnalignedCycleMeanTraces, 3));
Std3DFluoAP = NaN(size(UnalignedCycleMeanTraces, 1)+50, size(UnalignedCycleMeanTraces, 2),...
    size(UnalignedCycleMeanTraces, 3));
NumNuc3DAP = NaN(size(UnalignedCycleMeanTraces, 1)+50, size(UnalignedCycleMeanTraces, 2),...
    size(UnalignedCycleMeanTraces, 3));


AllFrameTimes = NaN(1, size(UnalignedCycleMeanTraces, 1)+50);
for NC = 9:14
    nc_idx = NC-8;
    if ~isempty(UnalignedCycleFrameTimes{nc_idx}) 
        
        if all(isnan(NCStartFrames))
            NCStartFrames(nc_idx) = 1;
            NCEndFrames(nc_idx) = length(UnalignedCycleFrameTimes{nc_idx});
            AllFrameTimes(NCStartFrames(nc_idx):NCEndFrames(nc_idx)) = UnalignedCycleFrameTimes{nc_idx};
            Ntimepoints = Ntimepoints+ length(UnalignedCycleFrameTimes{nc_idx});
        else
            NCStartFrames(nc_idx) = NCEndFrames(nc_idx-1);
            NCEndFrames(nc_idx) = length(UnalignedCycleFrameTimes{nc_idx})+NCStartFrames(nc_idx)-1;
            Ntimepoints = Ntimepoints+ length(UnalignedCycleFrameTimes{nc_idx})-1;
            AllFrameTimes(NCStartFrames(nc_idx):NCEndFrames(nc_idx)) = UnalignedCycleFrameTimes{nc_idx}+...
                AllFrameTimes(NCStartFrames(nc_idx));
        end
        MeanFluoAP(NCStartFrames(nc_idx):NCEndFrames(nc_idx),:,nc_idx) =...
            UnalignedCycleMeanTraces(1:length(UnalignedCycleFrameTimes{nc_idx}),:,nc_idx);
        StdFluoAP(NCStartFrames(nc_idx):NCEndFrames(nc_idx),:,nc_idx) =...
            UnalignedCycleTraceStdErrors(1:length(UnalignedCycleFrameTimes{nc_idx}),:,nc_idx);
        NumNucAP(NCStartFrames(nc_idx):NCEndFrames(nc_idx),:,nc_idx) =...
            UnalignedCycleNumNuclei(1:length(UnalignedCycleFrameTimes{nc_idx}),:,nc_idx);
        Mean3DFluoAP(NCStartFrames(nc_idx):NCEndFrames(nc_idx),:,nc_idx) =...
            Unaligned3DCycleMeanTraces(1:length(UnalignedCycleFrameTimes{nc_idx}),:,nc_idx);
        Std3DFluoAP(NCStartFrames(nc_idx):NCEndFrames(nc_idx),:,nc_idx) =...
            Unaligned3DCycleTraceStdErrors(1:length(UnalignedCycleFrameTimes{nc_idx}),:,nc_idx);
        NumNuc3DAP(NCStartFrames(nc_idx):NCEndFrames(nc_idx),:,nc_idx) =...
            Unaligned3DCycleNumNuclei(1:length(UnalignedCycleFrameTimes{nc_idx}),:,nc_idx);
        
    end
    
end



MeanFluoAP = nansum(MeanFluoAP(~isnan(AllFrameTimes),:,:), 3);
StdFluoAP = nansum(StdFluoAP(~isnan(AllFrameTimes),:,:), 3);
NumNucAP = nansum(NumNucAP(~isnan(AllFrameTimes),:,:), 3);
Mean3DFluoAP = nansum(Mean3DFluoAP(~isnan(AllFrameTimes),:,:), 3);
Std3DFluoAP = nansum(Std3DFluoAP(~isnan(AllFrameTimes),:,:), 3);
NumNuc3DAP = nansum(NumNuc3DAP(~isnan(AllFrameTimes),:,:), 3);
AllFrameTimes = AllFrameTimes(~isnan(AllFrameTimes))/60;
nc_info = NCStartFrames;

close all

FrameTraceFig = figure(1);
FrameTraceAx = axes(FrameTraceFig);
NFrames = 1:size(MeanFluoAP, 1);
linehandles = {};
for j=1:length(nuclear_cycles)
    if nc_info(nuclear_cycles(j)-8) > 0
        linehandles{j} = xline(AllFrameTimes(nc_info(nuclear_cycles(j)-8)), 'r--');
        set(get(get(linehandles{j} , 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        hold on
    end
end

eb = errorbar(AllFrameTimes, ones(1, length(AllFrameTimes)), .1*ones(1, length(AllFrameTimes)), 'vertical', 'LineStyle', 'none');
hold on
set(eb, 'color', colors(4,:), 'LineWidth', 1);
set(get(get(eb, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
plot(AllFrameTimes, ones(1, length(AllFrameTimes)), '.-', 'Color', colors(4,:));
eb2 = errorbar(AllFrameTimes, ones(1, length(AllFrameTimes)), .1*ones(1, length(AllFrameTimes)), 'vertical', 'LineStyle', 'none');
set(eb2, 'color', colors(5,:), 'LineWidth', 1);
set(get(get(eb2, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
plot(AllFrameTimes, ones(1, length(AllFrameTimes)), '.-', 'Color', colors(5,:));
hold off
xlabel('Time (min)')
ylabel('Fluo (AU)')
ylim([0, max([max(max(MeanFluoAP+StdFluoAP)),max(max(Mean3DFluoAP+Std3DFluoAP)),1])])
xlim([0, max(AllFrameTimes)])
title(FrameTraceAx, {Prefix_label,...
    ['Nuclear Cycle ', num2str(nc), ', AP bin: ',num2str(0)]})

for i = 1:size(MeanFluoAP,2)
    if ~sum(~isnan(MeanFluoAP(:, i)))
        continue
    end
    
    use_idx = NumNucAP(:,i) >= min_nuc_count;
    if sum(use_idx) == 0
        FrameTraceAx.Children(3).XData = FrameTimes/60;
        FrameTraceAx.Children(3).YData = zeros(1, length(FrameTimes));
        FrameTraceAx.Children(4).XData = FrameTimes/60;
        FrameTraceAx.Children(4).YData = zeros(1, length(FrameTimes));
        FrameTraceAx.Children(4).YPositiveDelta = .1*ones(1, length(FrameTimes));
        FrameTraceAx.Children(4).YNegativeDelta = .1*ones(1, length(FrameTimes));
    else
        FrameTraceAx.Children(3).YData = MeanFluoAP(use_idx,i);
        FrameTraceAx.Children(3).XData = AllFrameTimes(use_idx);
        FrameTraceAx.Children(4).YData = MeanFluoAP(use_idx,i);
        FrameTraceAx.Children(4).XData = AllFrameTimes(use_idx);
        FrameTraceAx.Children(4).YPositiveDelta = StdFluoAP(use_idx,i);
        FrameTraceAx.Children(4).YNegativeDelta  = StdFluoAP(use_idx, i);
    end
    
    use_idx2 = NumNuc3DAP(:,i) >= min_nuc_count;
    if sum(use_idx2) == 0
        FrameTraceAx.Children(1).XData = FrameTimes/60;
        FrameTraceAx.Children(1).YData = zeros(1, length(FrameTimes));
        FrameTraceAx.Children(2).XData = FrameTimes/60;
        FrameTraceAx.Children(2).YData = zeros(1, length(FrameTimes));
        FrameTraceAx.Children(2).YPositiveDelta = .1*ones(1, length(FrameTimes));
        FrameTraceAx.Children(2).YNegativeDelta = .1*ones(1, length(FrameTimes));
    else
        FrameTraceAx.Children(1).YData = Mean3DFluoAP(use_idx2,i);
        FrameTraceAx.Children(1).XData = AllFrameTimes(use_idx2);
        FrameTraceAx.Children(2).YData = Mean3DFluoAP(use_idx2,i);
        FrameTraceAx.Children(2).XData = AllFrameTimes(use_idx2);
        FrameTraceAx.Children(2).YPositiveDelta = Std3DFluoAP(use_idx2,i);
        FrameTraceAx.Children(2).YNegativeDelta  = Std3DFluoAP(use_idx2, i);
    end
    
    
    title(FrameTraceAx, {Prefix_label,...
        ['Nuclear Cycle ', num2str(nc), ', Fraction Embryo Length: ',num2str((i-1)*liveExperiment.APResolution)]})
    if sum(use_idx) == 0
        continue
    end
    hlegend = legend({'Fluo', 'Fluo3DGauss'}, 'Location', 'northeast',...
        'FontSize', 10);
    saveas(FrameTraceFig,[timedir3, filesep,'FullFluoTrace_bin', num2str(i),'.png']);
    
end
end
close all