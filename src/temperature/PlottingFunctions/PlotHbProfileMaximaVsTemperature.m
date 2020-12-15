function PlotHbProfileMaximaVsTemperature(this, AllMaxFluos, HalfMaxVectors, SetsContainMax, outdir, varargin)
x = 1;
while x < length(varargin)
    if strcmp(lower(varargin{x}), 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    end
    x = x+ 1;
end
if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    disp('Invalid choice of plotting colors. Can use either "default" or "pboc".') % change to error
end
if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
end
outdir2 = [outdir, filesep, 'NC14AnteriorMaxima'];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end
outdir3 = [outdir2, filesep, datestr(now, 'yyyymmdd')];
if ~exist(outdir3, 'dir')
    mkdir(outdir3)
end
Temp_obs = this.Temp_obs(this.ProcessedExperiments);
R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
IsAnterior = zeros(1, length(this.ProcessedExperiments));

for i = 1:length(this.ProcessedExperiments)
    if strcmp(lower(this.Regions{this.ProcessedExperiments(i)}), 'anterior')
        IsAnterior(i) = 1;
    end
end
%for nc_idx=1:length(this.IncludedNCs)
nc_idx = find(this.IncludedNCs == 14);

MaxFluos = AllMaxFluos{nc_idx};
HalfMaxVector = HalfMaxVectors{nc_idx};



close all
clear eb prof
MaximaFig = figure(1);
set(MaximaFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
set(gcf,'color','w');
% MaximaAx = axes(MaximaFig);
MaximaAx = subplot(1, 2, 1);
p1 =  scatter(MaximaAx, Temp_obs(SetsContainMax{nc_idx} == 1 & IsAnterior), MaxFluos(SetsContainMax{nc_idx} == 1 & IsAnterior),50,...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));
xlabel('Temperature (°C)', 'FontSize', 16)
ylabel('Maximum Nuclear Fluorescence (AU)')
xlim([15, 29])
MaximaAx.XAxis.FontSize = 16;
MaximaAx.YAxis.FontSize = 16;
MaximaAx2 = subplot(1, 2, 2);
p2 =  scatter(MaximaAx2, Temp_obs(SetsContainMax{nc_idx} == 1 & IsAnterior), HalfMaxVector(SetsContainMax{nc_idx} == 1 & IsAnterior),50,...
    'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', colors(2,:));
xlabel('Temperature (°C)', 'FontSize', 16)
ylabel('Fraction of Embryo Length at Hb Boundary ')
MaximaAx2.XAxis.FontSize = 16;
MaximaAx2.YAxis.FontSize = 16;
xlim([15, 29])
%ylim([13.5, 16.5])

saveas(MaximaFig,[outdir3, filesep, 'NC14Maxima.png']);



end