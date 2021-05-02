function GenerateTwoProjectComps(ResultsPaths, outdir)
close all




if ~exist(outdir, 'dir')
    mkdir(outdir)
end

R = 8.314*10^(-3);
UseSharedYAxis = false;
UseSharedYAxisFrequency = false;
UseSharedYAxisInit = false;

APResolution = 0.025;
APbins = (0:APResolution:1)*100;

CompiledParameters = CompileParameters(ResultsPaths);
NumSets = size(CompiledParameters.InitiationRates, 1);
IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 1)) > 0);
MinAPbin = min(IncludedAPbins);
MaxAPbin = max(IncludedAPbins);

NumTimeBins  =  size(CompiledParameters.InitiationRates, 2);
IncludedTimeBins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 2)) > 0).';

InitYmin = nanmin(nanmin(nanmin(CompiledParameters.InitiationRates+CompiledParameters.InitiationRatesStdErr)));
InitYmax = nanmax(nanmax(nanmax(CompiledParameters.InitiationRates+CompiledParameters.InitiationRatesStdErr)));

InitYmin = max(floor(InitYmin/50)*50, 0);
InitYmax = ceil(InitYmax/50)*50;

DurYmin = nanmin(nanmin(nanmin(CompiledParameters.Durations+CompiledParameters.DurationsStdErr)));
DurYmax = nanmax(nanmax(nanmax(CompiledParameters.Durations+CompiledParameters.DurationsStdErr)));

DurYmin = max(floor(DurYmin/2)*2, 0);
DurYmax = min(ceil(DurYmax/2)*2, 20);


FreqYmin = nanmin(nanmin(nanmin(CompiledParameters.Frequencies+CompiledParameters.FrequenciesStdErr)));
FreqYmax = nanmax(nanmax(nanmax(CompiledParameters.Frequencies+CompiledParameters.FrequenciesStdErr)));

FreqYmin = max(floor(FreqYmin/0.1)*0.1, 0);
FreqYmax = min(ceil(FreqYmax/0.1)*0.1, 4);


SubplotDims = numSubplots(length(ResultsPaths));
SubplotDims(2) = SubplotDims(2)+1;

SubFigDims = [0.85, 0.8];% 0.9*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
SubFigDims(2) = min([SubFigDims(2), 0.95]);
SubFigDims = round(SubFigDims, 2);

LegendSubplots = (1:SubplotDims(1))*SubplotDims(2);


LegendWidth = 0.2;
LegendXPosition = 0.78;
SubplotXBuffer = 0.05;
SubplotYBuffer = 0.1;

SubplotWidth = .23; %(LegendXPosition-(0.05-SubplotXBuffer))/(SubplotDims(2)-1);
SubplotHeight = .4;  %(1-(0.15-SubplotYBuffer))/SubplotDims(1);
SubplotXPositions = 0.025+(0:(SubplotDims(2)-2))*(SubplotWidth);
SubplotYPositions = 0.0+(0:(SubplotDims(1)-1))*(SubplotHeight);

SubplotYPositions = [.050, .5];


%%
FigAx = cell(1, NumSets);

dur_fig = figure(1);

set(dur_fig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
set(gcf,'color','w');
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);
hold on
SubplotIndex = 0;
SubplotIndexList = zeros(1, NumSets);
for SetIndex = 1:NumSets
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        FigAx{SubplotIndex} = subplot(2,4,(SetIndex-1)*4+1, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        FigAx{SubplotIndex} = subplot(2,4, (SetIndex-1)*4+1);
    end
    SubplotIndexList(SubplotIndex) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Durations(SetIndex,:,:));
    y_val_errs = squeeze(CompiledParameters.DurationsStdErr(SetIndex,:,:));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(t,:)))
            continue
        end
        set_ids = ~isnan(y_vals(t,:));
        
        x = APbins(set_ids);
        y = y_vals(t, set_ids);
        y_err = y_val_errs(t,set_ids);
        
        errorbar(FigAx{SubplotIndex}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(FigAx{SubplotIndex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on

    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FigAx{SubplotIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FigAx{SubplotIndex}, 'position', pos);
    
    xlabel('AP position')
    ylabel('burst duration (min)')
    set(FigAx{SubplotIndex} ,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
     title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C w', num2str(CompiledParameters.nSteps(SetIndex))]);
    newylim = get(FigAx{SubplotIndex}, 'ylim');
    newylim(1) = DurYmin;
    if newylim(2) > 15
        newylim(2)  = 15;
    end
%     set(FigAx{SetIndex}, 'ylim', newylim)
    
   SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        FigAx{SubplotIndex} = subplot( 2,4,(SetIndex-1)*4+2, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        FigAx{SubplotIndex} = subplot(2,4, (SetIndex-1)*4+2);
    end
    SubplotIndexList(SubplotIndex) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Frequencies(SetIndex,:,:));
    y_val_errs = squeeze(CompiledParameters.FrequenciesStdErr(SetIndex,:,:));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(t,:)))
            continue
        end
        set_ids = ~isnan(y_vals(t,:));
        
        x = APbins(set_ids);
        y = y_vals(t, set_ids);
        y_err = y_val_errs(t,set_ids);
        
        errorbar(FigAx{SubplotIndex}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(FigAx{SubplotIndex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FigAx{SubplotIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FigAx{SubplotIndex}, 'position', pos);
    
    xlabel('AP position')
    ylabel('burst frequency (1/min)')
    set(FigAx{SubplotIndex} ,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
   title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C w', num2str(CompiledParameters.nSteps(SetIndex))]);
    newylim = get(FigAx{SubplotIndex}, 'ylim');
    newylim(1) = FreqYmin;
    if newylim(2) > 2
        newylim(2)  = 2;
    end
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        FigAx{SubplotIndex} = subplot(2,4, (SetIndex-1)*4+3, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        FigAx{SubplotIndex} = subplot(2,4, (SetIndex-1)*4+3);
    end
    SubplotIndexList(SubplotIndex) = SubplotIndex;
    
    y_vals = squeeze(CompiledParameters.InitiationRates(SetIndex,:,:));
    y_val_errs = squeeze(CompiledParameters.InitiationRatesStdErr(SetIndex,:,:));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(t,:)))
            continue
        end
        set_ids = ~isnan(y_vals(t,:));
        
        x = APbins(set_ids);
        y = y_vals(t, set_ids);
        y_err = y_val_errs(t,set_ids);
        
        errorbar(FigAx{SubplotIndex}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(FigAx{SubplotIndex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    

    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FigAx{SubplotIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FigAx{SubplotIndex}, 'position', pos);
    
    xlabel('AP position')
    ylabel('initiation rate (au/min)')
    set(FigAx{SubplotIndex} ,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C']);
    newylim = get(FigAx{SubplotIndex}, 'ylim');
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);
set(LegendAx, 'Clim', [1, length(IncludedTimeBins)])
% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
h = colorbar;
h.Ticks = 1:length(IncludedTimeBins);
h.TickLabels = CompiledParameters.TimeVector(IncludedTimeBins);
ylabel(h,'time cohort (minutes into nc14)')
hold off
set(LegendAx,'Fontsize',16)
axis off

pos = get(LegendAx, 'Position');
pos(1) = .6;
pos(3) = .2;
set(LegendAx, 'Position', pos);


for SubplotIndex = 1:7
    if isempty(FigAx{SubplotIndex})
        continue
    end
    if SubplotIndex < 5
        ColumnIndex = SubplotIndex;
        RowIndex = 1;
    else
        ColumnIndex = SubplotIndex-4;
        RowIndex = 2;
    end

    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FigAx{SubplotIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FigAx{SubplotIndex}, 'position', pos);
end

saveas(dur_fig,[outdir, 'CompProjectsubplots_allparams.png'])

%%
