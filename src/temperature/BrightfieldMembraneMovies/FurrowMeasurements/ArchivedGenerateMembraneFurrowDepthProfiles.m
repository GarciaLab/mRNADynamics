clear all
% T = 17.5C: 33-37 m 
% T = 20C: (27) 29-31.5 (33) m
% T = 22.5C: 19-22 m
% T = 25C: 16-18 m
% T = 27.5C: (13) 15-16 m

Prefixes = {'2021-08-25-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo1',...
    '2021-08-25-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...
    '2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo3',...
    '2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo1',...
    '2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo2',...
    '2021-08-26-BrightfieldMembraneFurrow-T25C-Embryo1',...
    };%,...
Prefixes2 = {'2022-04-11-BrightfieldMembraneFurrow-T25C-Embryo1',...
    '2022-04-11-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo1',... 
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...% Good Length NC13
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo3',...% Good Length NC13
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo4',...% Good Length NC13
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo5',...% Good Length NC13
    '2022-04-13-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo6',...% Good Length NC13
    '2022-04-18-BrightfieldMembraneFurrow-HisRFP-T27_5C-Embryo1',...% Good Length NC13
    '2022-04-18-BrightfieldMembraneFurrow-HisRFP-T27_5C-Embryo2',...% Good Length NC13
    '2022-04-18-BrightfieldMembraneFurrow-HisRFP-T27_5C-Embryo3',...% Good Length NC13
    '2022-04-19-BrightfieldMembraneFurrow-HisRFP-T22_5C-Embryo1',...% Good Length NC13
    '2022-04-19-BrightfieldMembraneFurrow-HisRFP-T22_5C-Embryo2',...% Good Length NC13
    '2022-04-17-BrightfieldMembraneFurrow-HisRFP-T20C-Embryo1',...% Good Length NC13
    '2022-04-17-BrightfieldMembraneFurrow-HisRFP-T20C-Embryo2',...
    '2022-04-17-BrightfieldMembraneFurrow-HisRFP-T20C-Embryo3',...% Good Length NC13
    '2022-04-15-BrightfieldMembraneFurrow-HisRFP-T17_5C-Embryo2'... (ok)
    '2022-04-15-BrightfieldMembraneFurrow-HisRFP-T17_5C-Embryo3',...% Good Length NC13
    '2022-04-16-BrightfieldMembraneFurrow-HisRFP-T17_5C-Embryo4'...% Good Length NC13
    };
Prefixes = [Prefixes Prefixes2];
% Prefixes = {'2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...
%     '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo3',...
%     '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo4',...
%     '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo5',...
%     '2022-04-13-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo6'};

T_sets = [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,25,  27.5, 27.5, 27.5, 22.5, 22.5, 20, 20, 20, 17.5, 17.5, 17.5];
T_obs = [25, 25, 25, 25, 25, 25, 25, 25.1, 25.1, 25.1, 25, 24.9, 24.9, 27.6, 27.5, 27.5, 22.5, 22.5, 20.3, 20.2, 20,17.7,  17.8, 17.8];
colors = brewermap(6,'Spectral');
markers = {'o', '^', '>', 's', 'd', '<', 'p', 'h'};
unique_temps = unique(T_sets);
color_index = [2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 4, 4, 5, 5, 5,  6, 6];
marker_index = zeros(1, length(T_sets));
for i = 1:6
    matches = find(color_index == i);
    marker_index(matches) = 1:length(matches);
end


liveExperiments = {};
APLengths =[];
DVLengths = [];
FrameInfos = {};
anaphaseFrames= {};
DeltaFCs = {};
StdDeltaFCs = {};
CountFCs = {};
NC14FrameTimes = {};


SEDeltaFCs = {};


dDelta_dts ={};
sigma_ts = {};


nc14s = [];
nc13s = [];        

for i = 1:length(Prefixes)
    disp([num2str(i),' of ', num2str(length(Prefixes))])
    liveExperiments{i} = LiveExperiment(Prefixes{i});
    if i == 5
        APLengths(i) = NaN;
        DVLengths(i) = NaN;
    else
        if isfile([liveExperiments{i}.resultsFolder,filesep,'APDetection.mat'])
            try
                load([liveExperiments{i}.resultsFolder,filesep,'APDetection.mat'], 'APLength', 'DVLength')
                APLengths(i)=  APLength;
                DVLengths(i)=  DVLength;
            catch
                [APLengths(i), DVLengths(i)] = GetAPAxisLength(Prefixes{i});
            end
            
        else
            [APLengths(i), DVLengths(i)] = GetAPAxisLength(Prefixes{i});
        end
    end
    FrameInfos{i} = getFrameInfo(liveExperiments{i});
    FrameTimes = [FrameInfos{i}(:).Time];
    anaphaseFrames{i} = liveExperiments{i}.anaphaseFrames;
    nc14s(i) =anaphaseFrames{i}(6);
    nc13s(i) =anaphaseFrames{i}(5);
    try
    NC13lengths(i) = (FrameTimes(nc14s(i))-FrameTimes(nc13s(i)))/60;
    catch
       NC13lengths(i) = NaN;
    end
    load([liveExperiments{i}.resultsFolder,filesep,'FurrowCanalDepthMeasurements.mat'])
    DeltaFCs{i} = DeltaFC_um(nc14s(i):end,1).';
    StdDeltaFCs{i} = DeltaFC_um(nc14s(i):end,2).';
    CountFCs{i} = NumFurrowMeasurements(nc14s(i):end);
    SEDeltaFCs{i} = StdDeltaFCs{i}./sqrt(CountFCs{i});
    NC14FrameTimes{i} = [FrameInfos{i}(nc14s(i):end).Time]/60;
    NC14FrameTimes{i} = NC14FrameTimes{i}-min(NC14FrameTimes{i});
    deriv = diff(DeltaFCs{i})./diff(NC14FrameTimes{i});
    dDelta_dts{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    sigma_ts{i} = SEDeltaFCs{i}./abs(dDelta_dts{i});
    sigma_ts{i}(sigma_ts{i} > 10) = NaN;
end
mkdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements');
plotdir = 'S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/Figures';
mkdir(plotdir);
%%

close all
RawDeltaFCFig = figure(1);
set(RawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(RawDeltaFCFig);
for i = 1:length(Prefixes)
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{i}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
hold on 
end
grid on 
hold off
hlegend = legend(Prefixes, 'Location', 'eastoutside',...
            'FontSize', 14);

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
outpath = [plotdir, filesep, 'AllTemperaturesRawDeltaFCs_noErrorBar.png'];
saveas(RawDeltaFCFig,outpath);
%%
close all
T25CAllRawDeltaFCFig = figure(1);
set(T25CAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T25CAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 25
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 25ºC')
outpath = [plotdir, filesep, 'T25CRawDeltaFCs_noErrorBar.png'];
saveas(T25CAllRawDeltaFCFig,outpath);

%%
close all
T25CHisRFPAllRawDeltaFCFig = figure(1);
set(T25CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T25CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 25 & color_index(i) == 2 
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 25ºC')
outpath = [plotdir, filesep, 'T25CHisRFPRawDeltaFCs_noErrorBar.png'];
saveas(T25CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T25CHisRFPAllRawDeltaFCFig = figure(1);
set(T25CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T25CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 25 & color_index(i) == 2 & i > 6
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 25ºC')
outpath = [plotdir, filesep, 'T25CUnsquishedHisRFPRawDeltaFCs_noErrorBar.png'];
saveas(T25CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T25CHisRFPAllRawDeltaFCFig = figure(1);
set(T25CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T25CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 25 & color_index(i) == 2 & i <= 6
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 25ºC')
outpath = [plotdir, filesep, 'T25CSquishedHisRFPRawDeltaFCs_noErrorBar.png'];
saveas(T25CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T25CHisRFPAllRawDeltaFCFig = figure(1);
set(T25CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T25CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 25 & color_index(i) == 3
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 25ºC')
outpath = [plotdir, filesep, 'T25CNoHisRawDeltaFCs_noErrorBar.png'];
saveas(T25CHisRFPAllRawDeltaFCFig,outpath);


%%
close all
T25CHisRFPAllRawDeltaFCFig = figure(1);
set(T25CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T25CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 25 & color_index(i) == 3 & i > 6
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 25ºC')
outpath = [plotdir, filesep, 'T25CNoHisUnsquishedRawDeltaFCs_noErrorBar.png'];
saveas(T25CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T25CHisRFPAllRawDeltaFCFig = figure(1);
set(T25CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T25CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 25 & color_index(i) == 3 & i <= 6
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 25ºC')
outpath = [plotdir, filesep, 'T25CNoHisSquishedRawDeltaFCs_noErrorBar.png'];
saveas(T25CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T27_5CHisRFPAllRawDeltaFCFig = figure(1);
set(T27_5CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T27_5CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 27.5 
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 27.5ºC')
outpath = [plotdir, filesep, 'T27_5CHisRFPUnsquishedRawDeltaFCs_noErrorBar.png'];
saveas(T27_5CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T22_5CHisRFPAllRawDeltaFCFig = figure(1);
set(T22_5CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T22_5CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 22.5 
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 22.5ºC')
outpath = [plotdir, filesep, 'T22_5CHisRFPUnsquishedRawDeltaFCs_noErrorBar.png'];
saveas(T22_5CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T20CHisRFPAllRawDeltaFCFig = figure(1);
set(T20CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T20CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 20 
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 20ºC')
outpath = [plotdir, filesep, 'T20CHisRFPUnsquishedRawDeltaFCs_noErrorBar.png'];
saveas(T20CHisRFPAllRawDeltaFCFig,outpath);

%%
close all
T17_5CHisRFPAllRawDeltaFCFig = figure(1);
set(T17_5CHisRFPAllRawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(T17_5CHisRFPAllRawDeltaFCFig);
for i = 1:length(Prefixes)
    if T_sets(i) == 17.5
    scatter(NC14FrameTimes{i}, DeltaFCs{i}, 20, 'filled', markers{marker_index(i)}, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(color_index(i),:));
hold on
    end
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
title('T = 17.5ºC')
outpath = [plotdir, filesep, 'T17_5CHisRFPUnsquishedRawDeltaFCs_noErrorBar.png'];
saveas(T17_5CHisRFPAllRawDeltaFCFig,outpath);



%%
close all
RawDeltaFCFig2 = figure(2);
set(RawDeltaFCFig2,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax2 = axes(RawDeltaFCFig2);
for i = 1:length(Prefixes)
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(NC14FrameTimes{i}, DeltaFCs{i},SEDeltaFCs{i},[markers{marker_index(i)}, '-'],...
      'MarkerSize', 10,...
      'MarkerEdgeColor', 'k','MarkerFaceColor',colors(color_index(i),:),...
      'Color', colors(color_index(i),:));%,'LineStyle', 'none')
 
 
 hold on 
end
grid on 
hold off

xlabel('Time into cycle 14 (minutes)', 'FontSize', 16)
ylabel('Membrane Furrow Depth (microns)', 'FontSize', 16)
%xlim([0, 65])
set(RawDeltaFCs_ax2,'Fontsize',18)
outpath = [plotdir, filesep, 'RawDeltaFCs.png'];
saveas(RawDeltaFCFig2,outpath);
%%
MaxDeltaShift =1;
MaxTShift = 10;

l = max([length(DeltaFCs{1}), length(DeltaFCs{2}), length(DeltaFCs{3}),...
    length(DeltaFCs{4}),  length(DeltaFCs{5}), length(DeltaFCs{6})]);

MaxT_wt = floor(max([NC14FrameTimes{3},NC14FrameTimes{4},NC14FrameTimes{5},NC14FrameTimes{6},NC14FrameTimes{1},NC14FrameTimes{2}]));
InterpolatedDeltaFCs = cell(1,6);
InterpolatedSEs = cell(1,6);
Interpolated_dDelta_dts =cell(1,6);
Interpolated_sigma_ts = cell(1,6);
InterpolatedTimes = 0:MaxT_wt;
for i = 1:6
    InterpolatedDeltaFCs{i} = interp1(NC14FrameTimes{i},DeltaFCs{i}, 0:MaxT_wt);
    InterpolatedSEs{i} = NaN(1,length(InterpolatedDeltaFCs{i}));
    for j = 1:length(InterpolatedTimes)
        t = InterpolatedTimes(j);
        Tindex = find(round(NC14FrameTimes{i},5) == t);
        if isempty(Tindex)
            Tlow = find(NC14FrameTimes{i} < t, 1, 'last');
            Thigh = find(NC14FrameTimes{i} > t, 1);
            if ~isempty(Tlow) & ~isempty(Thigh)
                t1 = NC14FrameTimes{i}(Tlow);
                s1  = SEDeltaFCs{i}(Tlow);
                t2 = NC14FrameTimes{i}(Thigh);
                s2  = SEDeltaFCs{i}(Thigh);
                dfdt1 = (t2-t)/(t2-t1);
                dfdt2 = (t-t1)/(t2-t1);
                InterpolatedSEs{i}(j) = sqrt(dfdt1^2*s1^2+dfdt2^2*s2^2);
            end
        else 
           InterpolatedSEs{i}(j) =  SEDeltaFCs{i}(Tindex);
        end
    end
    
    deriv = diff(InterpolatedDeltaFCs{i})./diff(InterpolatedTimes);
    Interpolated_dDelta_dts{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    Interpolated_sigma_ts{i} = InterpolatedSEs{i}./abs(Interpolated_dDelta_dts{i});
    Interpolated_sigma_ts{i}(Interpolated_sigma_ts{i} > 10) = NaN;
    
end
%%

RawDeltaFCFig3 = figure(3);
set(RawDeltaFCFig3,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax3 = axes(RawDeltaFCFig3);

%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(InterpolatedTimes, InterpolatedDeltaFCs{3},InterpolatedSEs{3},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:),...
      'Color', colors(1,:));%,'LineStyle', 'none')
 
 
 hold on 
  errorbar(InterpolatedTimes, InterpolatedDeltaFCs{4},InterpolatedSEs{4},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(2,:),'MarkerFaceColor',colors(2,:),...
      'Color', colors(2,:));%,'LineStyle', 'none')
   errorbar(InterpolatedTimes, InterpolatedDeltaFCs{5},InterpolatedSEs{5},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(3,:),'MarkerFaceColor',colors(3,:),...
      'Color', colors(3,:));%,'LineStyle', 'none')
   errorbar(InterpolatedTimes, InterpolatedDeltaFCs{6},InterpolatedSEs{6},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(4,:),'MarkerFaceColor',colors(4,:),...
      'Color', colors(4,:));%,'LineStyle', 'none')
 
    errorbar(InterpolatedTimes, InterpolatedDeltaFCs{1},InterpolatedSEs{1},'o-','MarkerSize', 5,...
      'MarkerEdgeColor', 'k','MarkerFaceColor',colors(5,:),...
      'Color', colors(5,:));%,'LineStyle', 'none')
   errorbar(InterpolatedTimes, InterpolatedDeltaFCs{2},InterpolatedSEs{2},'o-','MarkerSize', 5,...
      'MarkerEdgeColor', 'k','MarkerFaceColor',colors(6,:),...
      'Color', colors(6,:));%,'LineStyle', 'none')
 

grid on
hold off
xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 65])
set(RawDeltaFCs_ax3,'Fontsize',18)
outpath = [plotdir, filesep, 'InterpolatedDeltaFCs.png'];
saveas(RawDeltaFCFig3,outpath);

%%

y_temp = NaN(4, l);
for i =3:6
    y_temp(i-2,15:length(InterpolatedDeltaFCs{i})) = InterpolatedDeltaFCs{i}(15:end);
end
%%
MatSize = [4,3*l-2];
LengthY = l;
fun = @(x)DeltaFC_Chi2(x,y_temp, MatSize, LengthY);
x0 = [0, 0, 0, 0, 0, 0];
lb = [-10, -10, -10, -1, -1, -1];
ub = [10,10,10,1,1,1];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = [1,2,3];
x_wt= ga(fun,6, A,b,Aeq,beq,lb,ub, nonlcon, intcon)


saveVars = [saveVars, 'x_wt'];
WildTypeProfile = {};
MinTime = 0;
wt_Tshifts = ([x_wt(1), x_wt(2), x_wt(3) 0]);
MaxTime =  MaxT_wt+(max(wt_Tshifts)-min(wt_Tshifts)); 
WildTypeProfile.Times = 0:MaxTime;
RefTimes = InterpolatedTimes-min([wt_Tshifts 0]);
ShiftedProfiles = NaN(4,length(WildTypeProfile.Times));
ShiftedProfileSEs = NaN(4,length(WildTypeProfile.Times));
ShiftedProfiles(1,RefTimes+1) = InterpolatedDeltaFCs{3};
ShiftedProfileSEs(1,RefTimes+1) = InterpolatedSEs{3};
for i = 4:6
    ShiftedProfiles(i-2,RefTimes+1+x_wt(i-3)) = InterpolatedDeltaFCs{i}+x_wt(i);
    ShiftedProfileSEs(i-2,RefTimes+1+x_wt(i-3)) = InterpolatedSEs{i};
end

WTProfile = mean(ShiftedProfiles,1,'omitnan');
WildTypeProfile.Profile = WTProfile;
WTProfileSE = std(ShiftedProfiles,1,'omitnan')./sqrt(sum(~isnan(ShiftedProfiles),1));
WildTypeProfile.SE =WTProfileSE;
%WildTypeProfile.Times 

deriv = diff(WTProfile)./diff(WildTypeProfile.Times);
    WTDiff  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    WTProfileSigmaT = WTProfileSE./abs(WTDiff);
    WTProfileSigmaT(WTProfileSigmaT > 10) = NaN;
   WildTypeProfile.SigmaT = WTProfileSigmaT; 
   saveVars = [saveVars, 'WildTypeProfile'];
% 
%  scatter(InterpolatedTimes, InterpolatedDeltaFCs{3},'k.')
%  hold on 
%   scatter(InterpolatedTimes+x_wt(1), InterpolatedDeltaFCs{4}+x_wt(4),'r.')
%  scatter(InterpolatedTimes+x_wt(2), InterpolatedDeltaFCs{5}+x_wt(5),'g.') 
%   scatter(InterpolatedTimes+x_wt(3), InterpolatedDeltaFCs{6}+x_wt(6),'b.') 
%  
% 
% %    scatter(NC14FrameTimes{3}, DeltaFCs{3}(nc14s(3):end),'g.')
% %     scatter(NC14FrameTimes{4}, DeltaFCs{4}(nc14s(4):end),'b.')
% % 
% % scatter(NC14FrameTimes{5}, DeltaFCs{5}(nc14s(5):end),'c.')
% % scatter(NC14FrameTimes{6}, DeltaFCs{6}(nc14s(6):end),'m.')
% hold off
% xlabel('Time into cycle 14 (minutes)')
% ylabel('Membrane Furrow Depth (microns)')


RawDeltaFCFig4 = figure;
set(RawDeltaFCFig4,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax4 = axes(RawDeltaFCFig4);

%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(InterpolatedTimes-min(wt_Tshifts), InterpolatedDeltaFCs{3},InterpolatedSEs{3},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:),...
      'Color', colors(1,:));%,'LineStyle', 'none')
 
 
 hold on 
  errorbar(InterpolatedTimes+x_wt(1)-min(wt_Tshifts), InterpolatedDeltaFCs{4}+x_wt(4),InterpolatedSEs{4},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(2,:),'MarkerFaceColor',colors(2,:),...
      'Color', colors(2,:));%,'LineStyle', 'none')
   errorbar(InterpolatedTimes+x_wt(2)-min(wt_Tshifts), InterpolatedDeltaFCs{5}+x_wt(5),InterpolatedSEs{5},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(3,:),'MarkerFaceColor',colors(3,:),...
      'Color', colors(3,:));%,'LineStyle', 'none')
   errorbar(InterpolatedTimes+x_wt(3)-min(wt_Tshifts), InterpolatedDeltaFCs{6}+x_wt(6),InterpolatedSEs{6},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(4,:),'MarkerFaceColor',colors(4,:),...
      'Color', colors(4,:));%,'LineStyle', 'none')

grid on



   errorbar(WildTypeProfile.Times, WTProfile,WTProfileSE,'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax4,'Fontsize',18)
outpath = [plotdir, filesep, 'ShiftedDeltaFCs_wt.png'];
saveas(RawDeltaFCFig4,outpath);
%%



%%

y_temp2 = NaN(2, l);
for i =1:2
    y_temp2(i,15:length(InterpolatedDeltaFCs{i})) = InterpolatedDeltaFCs{i}(15:end);
end
%%
MatSize = [2,3*l-2];
LengthY = l;
fun2 = @(x)DeltaFC_Chi2(x,y_temp2, MatSize, LengthY);
lb = [-2, -1];
ub = [2,1];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1;
x_his= ga(fun2,2, A,b,Aeq,beq,lb,ub, nonlcon, intcon)

saveVars = [saveVars, 'x_his'];
HisRFPProfile = {};
MinTime = 0;
his_Tshifts = ([x_his(1) 0 ]);
MaxTime =  MaxT_wt+(max(his_Tshifts)-min(his_Tshifts)); 
HisRFPProfile.Times = 0:MaxTime;
RefTimes = InterpolatedTimes-min([his_Tshifts 0]);
ShiftedProfiles = NaN(2,length(HisRFPProfile.Times));
ShiftedProfileSEs = NaN(2,length(HisRFPProfile.Times));
ShiftedProfiles(1,RefTimes+1) = InterpolatedDeltaFCs{1};
ShiftedProfileSEs(1,RefTimes+1) = InterpolatedSEs{1};

ShiftedProfiles(2,RefTimes+1+x_his(1)) = InterpolatedDeltaFCs{2}+x_his(2);
ShiftedProfileSEs(2,RefTimes+1+x_his(1)) = InterpolatedSEs{2};


HisProfile = mean(ShiftedProfiles,1,'omitnan');
HisRFPProfile.Profile = HisProfile;
HisProfileSE = std(ShiftedProfiles,1,'omitnan')./sqrt(sum(~isnan(ShiftedProfiles),1));
HisRFPProfile.SE =HisProfileSE;
%WildTypeProfile.Times 

deriv = diff(HisProfile)./diff(HisRFPProfile.Times);
    HisDiff  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    HisProfileSigmaT = HisProfileSE./abs(HisDiff);
    HisProfileSigmaT(HisProfileSigmaT > 10) = NaN;
   HisRFPProfile.SigmaT = HisProfileSigmaT; 
     saveVars = [saveVars, 'HisRFPProfile'];
% 
%  scatter(InterpolatedTimes, InterpolatedDeltaFCs{3},'k.')
%  hold on 
%   scatter(InterpolatedTimes+x_wt(1), InterpolatedDeltaFCs{4}+x_wt(4),'r.')
%  scatter(InterpolatedTimes+x_wt(2), InterpolatedDeltaFCs{5}+x_wt(5),'g.') 
%   scatter(InterpolatedTimes+x_wt(3), InterpolatedDeltaFCs{6}+x_wt(6),'b.') 
%  
% 
% %    scatter(NC14FrameTimes{3}, DeltaFCs{3}(nc14s(3):end),'g.')
% %     scatter(NC14FrameTimes{4}, DeltaFCs{4}(nc14s(4):end),'b.')
% % 
% % scatter(NC14FrameTimes{5}, DeltaFCs{5}(nc14s(5):end),'c.')
% % scatter(NC14FrameTimes{6}, DeltaFCs{6}(nc14s(6):end),'m.')
% hold off
% xlabel('Time into cycle 14 (minutes)')
% ylabel('Membrane Furrow Depth (microns)')


RawDeltaFCFig5 = figure;
set(RawDeltaFCFig5,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax5 = axes(RawDeltaFCFig5);

%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(InterpolatedTimes-min(his_Tshifts), InterpolatedDeltaFCs{1},InterpolatedSEs{1},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(5,:),'MarkerFaceColor',colors(5,:),...
      'Color', colors(5,:));%,'LineStyle', 'none')
 
 
 hold on 
  errorbar(InterpolatedTimes+x_his(1)-min(his_Tshifts), InterpolatedDeltaFCs{2}+x_his(2),InterpolatedSEs{2},'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(6,:),'MarkerFaceColor',colors(6,:),...
      'Color', colors(6,:));%,'LineStyle', 'none')


grid on



   errorbar(HisRFPProfile.Times, HisProfile,HisProfileSE,'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax5,'Fontsize',18)
outpath = [plotdir, filesep, 'ShiftedDeltaFCs_histone.png'];
saveas(RawDeltaFCFig5,outpath);


%%
l  = max([length(WildTypeProfile.Times), length(HisRFPProfile.Times)]);
MaxT =  max([max(WildTypeProfile.Times), max(HisRFPProfile.Times)]);
y_temp3 = NaN(2, l);
y_temp3(1,1:length(HisRFPProfile.Times)) = HisRFPProfile.Profile;
y_temp3(2,1:length(WildTypeProfile.Times)) = WildTypeProfile.Profile;

%%
MatSize = [2,3*l-2];
LengthY = l;
fun2 = @(x)DeltaFC_Chi2(x,y_temp3, MatSize, LengthY);
lb = [-20, -1];
ub = [20,1];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1;
x_final= ga(fun2,2, A,b,Aeq,beq,lb,ub, nonlcon, intcon)
saveVars = [saveVars, 'x_final'];
%%
WildTypeProfile.InitProfile = WildTypeProfile.Profile ;
WildTypeProfile.InitProfileSE = WildTypeProfile.SE ;
WildTypeProfile.InitProfileSigmaT =WildTypeProfile.SigmaT;
WildTypeProfile.InitTimes = WildTypeProfile.Times;

WildTypeProfile.Times = 0:(max(WildTypeProfile.InitTimes)+x_final(1));
WildTypeProfile.Profile = NaN(1,length(WildTypeProfile.Times));
WildTypeProfile.SE = NaN(1,length(WildTypeProfile.Times));
WildTypeProfile.SigmaT = NaN(1,length(WildTypeProfile.Times));
if x_final(1) >= 0
    WildTypeProfile.Profile((1+x_final(1)):(length(WildTypeProfile.InitProfile)+x_final(1))) =...
        WildTypeProfile.InitProfile+x_final(2);
    WildTypeProfile.SE((1+x_final(1)):(length(WildTypeProfile.InitProfile)+x_final(1))) =...
        WildTypeProfile.InitProfileSE;
    WildTypeProfile.SigmaT((1+x_final(1)):(length(WildTypeProfile.InitProfile)+x_final(1))) =...
        WildTypeProfile.InitProfileSigmaT;

else
   WildTypeProfile.Profile(1: (length(WildTypeProfile.InitProfile)+x_final(1)))  = ...
       WildTypeProfile.InitProfile(1-x_final(1):end)+x_final(2);
   WildTypeProfile.SE(1: (length(WildTypeProfile.InitProfile)+x_final(1)))  = ...
       WildTypeProfile.InitProfileSE(1-x_final(1):end);
   WildTypeProfile.SigmaT(1: (length(WildTypeProfile.InitProfile)+x_final(1)))  = ...
       WildTypeProfile.InitProfileSigmaT(1-x_final(1):end);
end
while isnan(WildTypeProfile.Profile(end))
    WildTypeProfile.SE = WildTypeProfile.SE(1:end-1);
    WildTypeProfile.SigmaT = WildTypeProfile.SigmaT(1:end-1);
    WildTypeProfile.Times = WildTypeProfile.Times(1:end-1);
    WildTypeProfile.Profile = WildTypeProfile.Profile(1:end-1);
end

%%




RawDeltaFCFig6 = figure;
set(RawDeltaFCFig6,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax6 = axes(RawDeltaFCFig6);

%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...

   errorbar(HisRFPProfile.Times, HisProfile,HisProfileSE,'.-','MarkerSize', 20,...
      'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:),...
      'Color', colors(1,:));%,'LineStyle', 'none')
 
 hold on 


grid on



   errorbar(WildTypeProfile.Times, WildTypeProfile.Profile,WildTypeProfile.SE,'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  hold off

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax6,'Fontsize',18)
outpath = [plotdir, filesep, 'ShiftedDeltaFCs_FinalProfiles.png'];
saveas(RawDeltaFCFig6,outpath);


save([outdir,filesep,'MembraneFurrowProfileInfo.mat'],saveVars{:});