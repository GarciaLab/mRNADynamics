
function PlotHbNC14TimeToMaxProfile(this, AllNCMaxSetTimes, SetsContainMax, outdir, varargin)
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
outdir2 = [outdir, filesep, 'NC14TimeToAnteriorMaximum'];
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
nc_idx = find(this.IncludedNCs == 14);
% if isempty(nc_idx)
%     error('NC 14 is not included in these data sets.')
% end
MaxSetTimes = AllNCMaxSetTimes{nc_idx};
TemperatureVector = -1./(Temp_obs(SetsContainMax{nc_idx} == 1 & IsAnterior)+273);
FitParams = polyfit(TemperatureVector,log(MaxSetTimes(SetsContainMax{nc_idx} == 1 & IsAnterior)),1);
alpha = FitParams(1);
A = FitParams(2);
Ea = alpha*R;
Subplot1DataYmin = min(MaxSetTimes(SetsContainMax{nc_idx} == 1& IsAnterior));
Subplot1DataYmax = max(MaxSetTimes(SetsContainMax{nc_idx} == 1& IsAnterior));
if Subplot1DataYmin/10 == floor(Subplot1DataYmin/10)
    Plot1Ymin = Subplot1DataYmin-10;
else
    Plot1Ymin =  floor(Subplot1DataYmin/10)*10;
end
if Subplot1DataYmax/10 == floor(Subplot1DataYmax/10)
    Plot1Ymax = Subplot1DataYmax+10;
else
    Plot1Ymax =  ceil(Subplot1DataYmax/10)*10;
end


Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;

Subplot2DataYmin = min(min(log(MaxSetTimes(SetsContainMax{nc_idx} == 1 & IsAnterior))), min(A+alpha*[Subplot2DataXmin,Subplot2DataXmax]));
Subplot2DataYmax = max(max(log(MaxSetTimes(SetsContainMax{nc_idx} == 1 & IsAnterior))), max(A+alpha*[Subplot2DataXmin,Subplot2DataXmax]));
Subplot2Yspan = Subplot2DataYmax-Subplot2DataYmin;
Plot2Ymin = Subplot2DataYmin - Subplot2Yspan*.05;
Plot2Ymax = Subplot2DataYmax + Subplot2Yspan*.05;

%%

%PlotXMin
%PlotXMax


TimeMaximaFig = figure(1);
set(TimeMaximaFig,'units', 'normalized', 'position',[0.01, 0.05, .8, .6]);
set(gcf,'color','w');
TimeMaximaAx1 = subplot(1, 2, 1);
p1 =  scatter(TimeMaximaAx1, Temp_obs(SetsContainMax{nc_idx} == 1 & IsAnterior), MaxSetTimes(SetsContainMax{nc_idx} == 1& IsAnterior),50,...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));
xlabel('Temperature (°C)', 'FontSize', 16)
ylabel('time (min)')
xlim([15, 29])
ylim([Plot1Ymin, Plot1Ymax])
TimeMaximaAx1.XAxis.FontSize = 16;
TimeMaximaAx1.YAxis.FontSize = 16;

TimeMaximaAx2 = subplot(1, 2, 2);
p2 =  scatter(TimeMaximaAx2, TemperatureVector, log(MaxSetTimes(SetsContainMax{nc_idx} == 1 & IsAnterior)),50,...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));
hold on
p3 = plot(TimeMaximaAx2, [Plot2Xmin, Plot2Xmax], A+alpha*[Plot2Xmin, Plot2Xmax],'-', 'Color', colors(1,:), 'MarkerEdgeColor', colors(1,:));
xlabel('-\alpha/(T + 273) (1/K)', 'FontSize', 16)
ylabel('log(time)', 'FontSize', 16)
TimeMaximaAx2.XAxis.FontSize = 16;
TimeMaximaAx2.YAxis.FontSize = 16;
ylim([Plot2Ymin, Plot2Ymax])
xlim([Plot2Xmin, Plot2Xmax])
%sgtitle( 'Measured time maximum Hunchback concentration', 'FontSize', 20)
legend([p2 p3], {'Data', ['Arrehnius fit, E_a = ', num2str(round(Ea, 2)),' kJ/mol']}, 'Location', 'northeast',...
    'FontSize', 14)

%end
saveas(TimeMaximaFig,[outdir3,  filesep, 'NC14TimeToMaxima.png']);
end