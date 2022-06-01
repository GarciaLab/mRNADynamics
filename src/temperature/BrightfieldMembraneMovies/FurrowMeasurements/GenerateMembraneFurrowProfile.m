% function FurrowProfileInfo = GenerateMembraneFurrowProfile(Prefixes,Outfile)
% author: Gabriella Martini
% date created: 5/31/22

% T = 17.5C: 33-37 m 
% T = 20C: (27) 29-31.5 (33) m
% T = 22.5C: 19-22 m
% T = 25C: 16-18 m
% T = 27.5C: (13) 15-16 m

Prefixes = {'2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo3',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo4',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo5',...
    '2022-04-13-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo6'};

SetLabel = 'T25C_UnsquishedHisRFPNoPV';


%%
NSets = length(Prefixes);
colors = brewermap(min(length(Prefixes), 8),'Spectral');
saveVars = {};


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
NC13durations = [];

for i = 1:NSets
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
    if ~isnan(nc14s(i)) & ~isnan(nc13s(i))
        NC13durations(i) = (FrameTimes(nc14s(i))-FrameTimes(nc13s(i)))/60;
    else
       NC13durations(i) = NaN;
    end
    load([liveExperiments{i}.resultsFolder,filesep,'FurrowCanalDepthMeasurements.mat'])
    DeltaFCs{i} = DeltaFCnoPV_um(nc14s(i):end,1).';
    StdDeltaFCs{i} = DeltaFCnoPV_um(nc14s(i):end,2).';
    CountFCs{i} = NumFurrowMeasurements(nc14s(i):end);
    SEDeltaFCs{i} = StdDeltaFCs{i}./sqrt(CountFCs{i});
    NC14FrameTimes{i} = [FrameInfos{i}(nc14s(i):end).Time]/60;
    NC14FrameTimes{i} = NC14FrameTimes{i}-min(NC14FrameTimes{i});
    deriv = diff(DeltaFCs{i})./diff(NC14FrameTimes{i});
    dDelta_dts{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    sigma_ts{i} = SEDeltaFCs{i}./abs(dDelta_dts{i});
    sigma_ts{i}(sigma_ts{i} > 10) = NaN;
end
if ~isdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements')
mkdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements');
end
plotdir = ['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/Figures',filesep, SetLabel];
if ~isdir(plotdir)
mkdir(plotdir);
end
savedir = ['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/',filesep, SetLabel];
if ~isdir(savedir)
mkdir(savedir);
end
%%

close all
RawDeltaFCFig = figure(1);
set(RawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(RawDeltaFCFig);
for i = 1:NSets
    plot(NC14FrameTimes{i}, DeltaFCs{i}, '-o', 'Color', colors(i,:), 'MarkerSize', 5,...
         'MarkerFaceColor', colors(i,:));
hold on 
end
grid on 
hold off
hlegend = legend(Prefixes, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
%xlim([0, 65])
set(RawDeltaFCs_ax,'Fontsize',18)
outpath = [plotdir, filesep, 'RawDeltaFCs_noErrorBar.png'];
saveas(RawDeltaFCFig,outpath);

%%
close all
RawDeltaFCFig2 = figure(2);
set(RawDeltaFCFig2,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax2 = axes(RawDeltaFCFig2);
for i = 1:NSets
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(NC14FrameTimes{i}, DeltaFCs{i},SEDeltaFCs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
 
 
 hold on 
end
grid on 
hold off
hlegend = legend(Prefixes, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)', 'FontSize', 16)
ylabel('Membrane Furrow Depth (microns)', 'FontSize', 16)
%xlim([0, 65])
set(RawDeltaFCs_ax2,'Fontsize',18)
outpath = [plotdir, filesep, 'RawDeltaFCs.png'];
saveas(RawDeltaFCFig2,outpath);
%%
MaxDeltaShift =1;
MaxTShift = 10;
LengthNC14s = zeros(1, NSets);
LengthNC14FrameTimes = zeros(1,NSets);
MaxDeltaMeasureFrames = zeros(1,NSets);
MaxDeltaMeasureTimes = zeros(1,NSets);
for i = 1:NSets
    GoodMeasurements = ~isnan(DeltaFCs{i});
    DeltaFCs{i} = DeltaFCs{i}(GoodMeasurements);
    NC14FrameTimes{i} = NC14FrameTimes{i}(GoodMeasurements);
    LengthNC14s(i) = length(DeltaFCs{i});
    LengthNC14FrameTimes(i) = max(NC14FrameTimes{i});

   
end
   
l = max(LengthNC14s);


MaxT_wt = floor(max(LengthNC14FrameTimes));
HalfT = MaxT_wt/2;


InterpolatedDeltaFCs = cell(1,NSets);
InterpolatedSEs = cell(1,NSets);
Interpolated_dDelta_dts =cell(1,NSets);
Interpolated_sigma_ts = cell(1,NSets);
InterpolatedTimes = 0:MaxT_wt;
for i = 1:NSets
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
close all
RawDeltaFCFig3 = figure(3);
set(RawDeltaFCFig3,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax3 = axes(RawDeltaFCFig3);
for i = 1:NSets
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(InterpolatedTimes, InterpolatedDeltaFCs{i},InterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
 
 
 hold on 
end
grid on 
hold off
hlegend = legend(Prefixes, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
x_limit = 5*ceil(max(InterpolatedTimes)/5);
xlim([0, x_limit])
set(RawDeltaFCs_ax3,'Fontsize',18)
outpath = [plotdir, filesep, 'InterpolatedDeltaFCs.png'];
saveas(RawDeltaFCFig3,outpath);

%%
lInt = length(InterpolatedTimes);
MinTimeToUse = find(InterpolatedTimes >= HalfT, 1);
y_temp = NaN(NSets, lInt);
for i = 1:NSets
    y_temp(i,MinTimeToUse:length(InterpolatedDeltaFCs{i})) = InterpolatedDeltaFCs{i}(MinTimeToUse:end);
end
%%
MatSize = [NSets,(NSets-1)*lInt-2];
LengthY = lInt;
fun = @(x)DeltaFC_Chi2(x,y_temp, MatSize, LengthY);
x0 = zeros(1, 2*(NSets-1));
lb = [-10*ones(1, NSets-1) -2*ones(1, NSets-1)];
ub = -1*lb;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1:(NSets-1);
x_wt= ga(fun,2*(NSets-1), A,b,Aeq,beq,lb,ub, nonlcon, intcon);


saveVars = [saveVars, 'x_wt'];
WildTypeProfile = {};
MinTime = 0;
%%
wt_Tshifts = zeros(1, NSets);
for i = 1:NSets-1
    wt_Tshifts(i) = x_wt(i);
end
AverageDeltaT = round(mean(wt_Tshifts));
StartMinDeltaT = min([wt_Tshifts 0]);
DiffDeltaT = AverageDeltaT - StartMinDeltaT;
MaxTime =  MaxT_wt+(max(wt_Tshifts)-min(wt_Tshifts)); 
WildTypeProfile.Times = 0:MaxTime;
RefTimes = InterpolatedTimes-min([wt_Tshifts 0]);
ShiftedProfiles = NaN(NSets,length(WildTypeProfile.Times));
ShiftedProfileSEs = NaN(NSets,length(WildTypeProfile.Times));
ShiftedProfiles(1,RefTimes+1) = InterpolatedDeltaFCs{1};
ShiftedProfileSEs(1,RefTimes+1) = InterpolatedSEs{1};
for i = 2:NSets
    ShiftedProfiles(i,RefTimes+1+x_wt(i-1)) = InterpolatedDeltaFCs{i}+x_wt(i-1+NSets-1);
    ShiftedProfileSEs(i,RefTimes+1+x_wt(i-1)) = InterpolatedSEs{i};
end


WildTypeProfile.Times = WildTypeProfile.Times-DiffDeltaT;

WTProfile = mean(ShiftedProfiles,1,'omitnan');
WildTypeProfile.Profile = WTProfile;
WTProfileSE = std(ShiftedProfiles,1,'omitnan')./sqrt(sum(~isnan(ShiftedProfiles),1));
WildTypeProfile.SE =WTProfileSE;
WTProfileSD = std(ShiftedProfiles,1,'omitnan');
WildTypeProfile.SD =WTProfileSD;
%WildTypeProfile.Times 

deriv = diff(WTProfile)./diff(WildTypeProfile.Times);
    WTDiff  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    WTProfileSigmaT = WTProfileSE./abs(WTDiff);
    WTProfileSigmaT(WTProfileSigmaT > 10) = NaN;
   WildTypeProfile.SigmaT = WTProfileSigmaT; 
   saveVars = [saveVars, 'WildTypeProfile'];
%%
close all
RawDeltaFCFig4 = figure;
set(RawDeltaFCFig4,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax4 = axes(RawDeltaFCFig4);
for i = 1:NSets
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
    if i == 1
  errorbar(InterpolatedTimes-min(wt_Tshifts)-DiffDeltaT, InterpolatedDeltaFCs{i},InterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    else
        errorbar(InterpolatedTimes-min(wt_Tshifts)-DiffDeltaT+x_wt(i-1), InterpolatedDeltaFCs{i}+x_wt(i+NSets-2),InterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    end
 
 
 hold on 
end





   errorbar(WildTypeProfile.Times, WTProfile,WTProfileSE,'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  
  grid on
  hold off
  
  Plot_Labels = [Prefixes 'Average Profile'];
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 80])
set(RawDeltaFCs_ax4,'Fontsize',18)
outpath = [plotdir, filesep, 'ShiftedDeltaFCs.png'];
saveas(RawDeltaFCFig4,outpath);

%%
close all
RawDeltaFCFig5 = figure;
set(RawDeltaFCFig5,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax5 = axes(RawDeltaFCFig5);

   scatter(WildTypeProfile.Times, WildTypeProfile.SE, 20,'filled',...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k');%,'LineStyle', 'none')
  






  grid on
  hold off
  
  Plot_Labels = [Prefixes 'Average Profile'];
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('$\sigma_{\delta_{FC}}\,\, (\mu m)$','interpreter','latex')
xlim([0, 80])
set(RawDeltaFCs_ax5,'Fontsize',18)
outpath = [plotdir, filesep, 'StdErrDeltaFCs.png'];
saveas(RawDeltaFCFig5,outpath);

%%
close all
RawDeltaFCFig5 = figure;
set(RawDeltaFCFig5,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax5 = axes(RawDeltaFCFig5);

   scatter(WildTypeProfile.Times, WildTypeProfile.SigmaT, 20,'filled',...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k');%,'LineStyle', 'none')
  






  grid on
  hold off
  
  Plot_Labels = [Prefixes 'Average Profile'];
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('$\sigma_{t}\,\, (m)$','interpreter','latex')
xlim([0, 80])
set(RawDeltaFCs_ax5,'Fontsize',18)
outpath = [plotdir, filesep, 'SigmaTDeltaFCs.png'];
saveas(RawDeltaFCFig5,outpath);


%%
save([savedir, filesep, 'MembraneFurrowProfile.mat'], 'saveVars');


%%
