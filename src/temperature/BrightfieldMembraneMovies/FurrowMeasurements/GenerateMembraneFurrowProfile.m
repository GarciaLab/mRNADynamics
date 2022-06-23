% function FurrowProfileInfo = GenerateMembraneFurrowProfile(Prefixes,Outfile)
% author: Gabriella Martini
% date created: 5/31/22

% T = 17.5C: 33-37 m 
% T = 20C: (27) 29-31.5 (33) m
% T = 22.5C: 19-22 m
% T = 25C: 16-18 m
% T = 27.5C: (13) 15-16 m
clear all
DubuisDataPath = 'S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/DubuisFigure2BMembraneData.csv';
DubuisFurrowData = csvread(DubuisDataPath,1, 0); 
DubuisTimes = DubuisFurrowData(:,1).';
DubuisMeanProfile = DubuisFurrowData(:,10).';
DubuisFixedMeanProfile = DubuisFurrowData(:,11).';

SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
load(SizeDataPath,'MeanAPLength','MeanDVLength');

APScale = MeanAPLength;
DVScale = MeanDVLength;
ComboScale = sqrt(MeanAPLength^2+MeanDVLength^2);


Prefixes = {'2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo3',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo4',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo5',...
    '2022-04-13-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo6'};

SetLabel = 'T25C_UnsquishedHisRFPNoPV';

Prefixes = {'2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo1',...
    '2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo2',...
    '2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo3',...
    '2021-08-26-BrightfieldMembraneFurrow-T25C-Embryo1'};

SetLabel = 'T25C_SquishedywNoPV';

Prefixes = {'2021-08-25-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo1',...
    '2021-08-25-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2'};

SetLabel = 'T25C_SquishedHisRFPNoPV';

Prefixes = {'2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo1',... % O
    '2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo2',...
    '2021-08-11-BrightfieldMembraneFurrow-T25C-Embryo3',...
    '2021-08-26-BrightfieldMembraneFurrow-T25C-Embryo1',...
    '2021-08-25-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo1',...
    '2021-08-25-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo2',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo3',...%'2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo4',...
    '2022-04-12-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo5',...
    '2022-04-13-BrightfieldMembraneFurrow-HisRFP-T25C-Embryo6'... 

    };
GroupStartIndices = [1, 5, 5,7];
GroupEndIndices = [4,10,6,10];
DirLabels = {'T25C_SquishedywNoPV',...
    'T25C_Squishedyw',...
    'T25C_AllHisRFPNoPV',...
    'T25C_AllHisRFP',...
        'T25C_SquishedHisRFPNoPV',...
    'T25C_SquishedHisRFP',...
    'T25C_UnsquishedHisRFPNoPV',...
    'T25C_UnsquishedHisRFP'...
    };

SetLabels = {'T25C Squished yw without PV',...
    'T25C Squished yw',...
    'T25C All HisRFP without PV',...
    'T25C All HisRFP',...
        'T25C Squished HisRFP without PV',...
    'T25C Squished HisRFP',...
    'T25C Unsquished HisRFP without PV',...
    'T25C Unsquished HisRFP'...
    };
 
NumGroups = 4;
    



%%
liveExperiments = cell(1, NumGroups);
APLengths = cell(1, NumGroups);
DVLengths = cell(1, NumGroups);
FrameInfos  = cell(1, NumGroups);
anaphaseFrames = cell(1, NumGroups); 
DeltaFCs = cell(2, NumGroups); 
StdDeltaFCs = cell(2, NumGroups); 
CountFCs = cell(2, NumGroups); 
NC14FrameTimes = cell(2, NumGroups); 


SEDeltaFCs = cell(2, NumGroups); 


dDelta_dts = cell(2, NumGroups); 
sigma_ts = cell(2, NumGroups); 


nc14s = cell(1, NumGroups); 
nc13s = cell(1, NumGroups); 
NC13durations = cell(1, NumGroups); 
InterpolatedDeltaFCs= cell(2, NumGroups); 
InterpolatedSEs = cell(2, NumGroups); 
Interpolated_dDelta_dts =cell(2, NumGroups); 
Interpolated_sigma_ts =cell(2, NumGroups); 
InterpolatedTimes =cell(2, NumGroups); 
WildTypeProfiles = cell(2, NumGroups); 
AverageDeltaTs = NaN(2, NumGroups);
AverageDeltaProfs = NaN(2, NumGroups);
x_wts =cell(2, NumGroups);
%%
for j = 1:NumGroups
    %%
NSets = GroupEndIndices(j)-GroupStartIndices(j)+1;
colors = brewermap(min(NSets, 8),'Spectral');
saveVars = {};


liveExperiments{j} = {};
APLengths{j} =[];
DVLengths{j} = [];
FrameInfos{j} = {};
anaphaseFrames{j}= {};


nc14s{j} = [];
nc13s{j} = []; 
NC13durations{j} = [];



NC14FrameTimes{1,j} = {};
DeltaFCs{1,j} = {};
StdDeltaFCs{1,j} = {};
CountFCs{1,j} = {};
SEDeltaFCs{1,j} = {};
dDelta_dts{1,j} ={};
sigma_ts{1,j} = {};

NC14FrameTimes{2,j} = {};
DeltaFCs{2,j} = {};
StdDeltaFCs{2,j} = {};
CountFCs{2,j} = {};
SEDeltaFCs{2,j} = {};
dDelta_dts{2,j} ={};
sigma_ts{2,j} = {};



for i = 1:NSets
    i_sub = i + GroupStartIndices(j)-1;
    disp(['Set ', num2str(j), ', ', num2str(i),' of ', num2str(NSets)])
    liveExperiments{j}{i} = LiveExperiment(Prefixes{i_sub});
    
    if isfile([liveExperiments{j}{i}.resultsFolder,filesep,'APDetection.mat'])
        try
            load([liveExperiments{j}{i}.resultsFolder,filesep,'APDetection.mat'], 'APLength', 'DVLength')
            APLengths{j}(i)=  APLength;
            DVLengths{j}(i)=  DVLength;
        catch
            [APLengths{j}(i), DVLengths{j}(i)] = GetAPAxisLength(Prefixes{i_sub});
        end
        
    else
        [APLengths{j}(i), DVLengths{j}(i)] = GetAPAxisLength(Prefixes{i_sub});
    end
    


    FrameInfos{j}{i} = getFrameInfo(liveExperiments{j}{i});
    FrameTimes = [FrameInfos{j}{i}(:).Time];
    anaphaseFrames{j}{i} = liveExperiments{j}{i}.anaphaseFrames;
    nc14s{j}(i) =anaphaseFrames{j}{i}(6);
    nc13s{j}(i) =anaphaseFrames{j}{i}(5);
    if ~isnan(nc14s{j}(i)) & ~isnan(nc13s{j}(i)) & (nc13s{j}(i) ~= 0)
        NC13durations{j}(i) = (FrameTimes(nc14s{j}(i))-FrameTimes(nc13s{j}(i)))/60;
    else
       NC13durations{j}(i) = NaN;
    end
    load([liveExperiments{j}{i}.resultsFolder,filesep,'FurrowCanalDepthMeasurements.mat'])
    DeltaFCs{1, j}{i} = ComboScale/sqrt(APLengths{j}(i)^2+DVLengths{j}(i)^2)*DeltaFCnoPV_um(nc14s{j}(i):end,1).';
    StdDeltaFCs{1, j}{i} = ComboScale/sqrt(APLengths{j}(i)^2+DVLengths{j}(i)^2)*DeltaFCnoPV_um(nc14s{j}(i):end,2).';
    CountFCs{1, j}{i} = NumFurrowMeasurements(nc14s{j}(i):end);
    SEDeltaFCs{1, j}{i} = StdDeltaFCs{1, j}{i}./sqrt(CountFCs{1, j}{i});
    NC14FrameTimes{1,j}{i} = [FrameInfos{j}{i}(nc14s{j}(i):end).Time]/60;
    NC14FrameTimes{1,j}{i} = NC14FrameTimes{1,j}{i}-min(NC14FrameTimes{1,j}{i});
    deriv = diff(DeltaFCs{1, j}{i})./diff(NC14FrameTimes{1,j}{i});
    dDelta_dts{1, j}{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    sigma_ts{1, j}{i} = SEDeltaFCs{1, j}{i}./abs(dDelta_dts{1, j}{i});
    sigma_ts{1, j}{i}(sigma_ts{1, j}{i} > 10) = NaN;
    
    DeltaFCs{2, j}{i} = ComboScale/sqrt(APLengths{j}(i)^2+DVLengths{j}(i)^2)*DeltaFC_um(nc14s{j}(i):end,1).';
    StdDeltaFCs{2, j}{i} = ComboScale/sqrt(APLengths{j}(i)^2+DVLengths{j}(i)^2)*DeltaFC_um(nc14s{j}(i):end,2).';
    CountFCs{2, j}{i} = NumFurrowMeasurements(nc14s{j}(i):end);
    NC14FrameTimes{2,j}{i} = [FrameInfos{j}{i}(nc14s{j}(i):end).Time]/60;
    NC14FrameTimes{2,j}{i} = NC14FrameTimes{2,j}{i}-min(NC14FrameTimes{2,j}{i});
    SEDeltaFCs{2, j}{i} = StdDeltaFCs{2, j}{i}./sqrt(CountFCs{2, j}{i});
    deriv = diff(DeltaFCs{2, j}{i})./diff(NC14FrameTimes{2,j}{i});
    dDelta_dts{2, j}{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    sigma_ts{2, j}{i} = SEDeltaFCs{2, j}{i}./abs(dDelta_dts{2, j}{i});
    sigma_ts{2, j}{i}(sigma_ts{2, j}{i} > 10) = NaN;
end
%%
if ~isdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements')
mkdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/');
end
plotdir_nopv = ['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/Figures',filesep, DirLabels{(j-1)*2+1}];
if ~isdir(plotdir_nopv)
mkdir(plotdir_nopv);
end
savedir_nopv = ['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements',filesep, DirLabels{(j-1)*2+1}];
if ~isdir(savedir_nopv)
mkdir(savedir_nopv);
end
%%

close all
RawDeltaFCFig = figure(1);
set(RawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(RawDeltaFCFig);
RawDeltaFCs_ax.FontSize = 14;
for i = 1:NSets
    plot(NC14FrameTimes{1,j}{i}, DeltaFCs{1, j}{i}, '-o', 'Color', colors(i,:), 'MarkerSize', 5,...
         'MarkerFaceColor', colors(i,:));
hold on 
end
grid on 
hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j))];
hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)', 'FontSize', 14)
ylabel('Membrane Furrow Depth (microns)', 'FontSize', 14)
%xlim([0, 65])
ylim([0 45])
title([SetLabels{2*(j-1)+1}, ' without pv membrane'], 'FontSize', 14)
set(RawDeltaFCs_ax,'Fontsize',14)
outpath = [plotdir_nopv, filesep, 'RawDeltaFCs_noErrorBar.png'];
saveas(RawDeltaFCFig,outpath);

%%
close all
RawDeltaFCFig2 = figure(2);
set(RawDeltaFCFig2,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax2 = axes(RawDeltaFCFig2);
RawDeltaFCs_ax2.FontSize = 14;
for i = 1:NSets
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(NC14FrameTimes{1,j}{i}, DeltaFCs{1,j}{i},SEDeltaFCs{1,j}{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
 
 
 hold on 
end
grid on 
hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j))];
hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)', 'FontSize', 16)
ylabel('Membrane Furrow Depth (microns)', 'FontSize', 16)
%xlim([0, 65])
ylim([0 45])
title([SetLabels{2*(j-1)+1}, ' without pv membrane'], 'FontSize', 14)
set(RawDeltaFCs_ax2,'Fontsize',14)
outpath = [plotdir_nopv, filesep, 'RawDeltaFCs.png'];
saveas(RawDeltaFCFig2,outpath);
%%
MaxDeltaShift =1;
MaxTShift = 10;
LengthNC14s = zeros(1, NSets);
LengthNC14FrameTimes = zeros(1,NSets);
MaxDeltaMeasureFrames = zeros(1,NSets);
MaxDeltaMeasureTimes = zeros(1,NSets);
for i = 1:NSets
    GoodMeasurements = ~isnan(DeltaFCs{1,j}{i});
    DeltaFCs{1,j}{i} = DeltaFCs{1,j}{i}(GoodMeasurements);
    NC14FrameTimes{1,j}{i} = NC14FrameTimes{1,j}{i}(GoodMeasurements);
    LengthNC14s(i) = length(DeltaFCs{1,j}{i});
    LengthNC14FrameTimes(i) = max(NC14FrameTimes{1,j}{i});

   
end
   
l = max(LengthNC14s);


MaxT_wt = floor(max(LengthNC14FrameTimes));
HalfT = MaxT_wt/2;
InterpolatedTimes{1,j} = 0:0.25:MaxT_wt;


InterpolatedDeltaFCs{1,j} = cell(1,NSets);
InterpolatedSEs{1,j} = cell(1,NSets);
Interpolated_dDelta_dts{1,j} =cell(1,NSets);
Interpolated_sigma_ts{1,j} = cell(1,NSets);

for i = 1:NSets
    InterpolatedDeltaFCs{1,j}{i} = interp1(NC14FrameTimes{1,j}{i},DeltaFCs{1,j}{i}, InterpolatedTimes{1,j});
    InterpolatedSEs{1,j}{i} = NaN(1,length(InterpolatedDeltaFCs{1,j}{i}));
    for k = 1:length(InterpolatedTimes{1,j})
        t = InterpolatedTimes{1,j}(k);
        Tindex = find(round(NC14FrameTimes{1,j}{i},5) == t);
        if isempty(Tindex)
            Tlow = find(NC14FrameTimes{1,j}{i} < t, 1, 'last');
            Thigh = find(NC14FrameTimes{1,j}{i} > t, 1);
            if ~isempty(Tlow) & ~isempty(Thigh)
                t1 = NC14FrameTimes{1,j}{i}(Tlow);
                s1  = SEDeltaFCs{1,j}{i}(Tlow);
                t2 = NC14FrameTimes{1,j}{i}(Thigh);
                s2  = SEDeltaFCs{1,j}{i}(Thigh);
                dfdt1 = (t2-t)/(t2-t1);
                dfdt2 = (t-t1)/(t2-t1);
                InterpolatedSEs{1,j}{i}(k) = sqrt(dfdt1^2*s1^2+dfdt2^2*s2^2);
            end
        else 
           InterpolatedSEs{1,j}{i}(k) =  SEDeltaFCs{1,j}{i}(Tindex);
        end
    end
    
    deriv = diff(InterpolatedDeltaFCs{1,j}{i})./diff(InterpolatedTimes{1,j});
    Interpolated_dDelta_dts{1,j}{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    Interpolated_sigma_ts{1,j}{i} = InterpolatedSEs{1,j}{i}./abs(Interpolated_dDelta_dts{1,j}{i});
    Interpolated_sigma_ts{1,j}{i}(Interpolated_sigma_ts{1,j}{i} > 10) = NaN;
    
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
  errorbar(InterpolatedTimes{1,j}, InterpolatedDeltaFCs{1,j}{i},InterpolatedSEs{1,j}{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
 
 
 hold on 
end
grid on 
hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j))];
hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
x_limit = 5*ceil(max(InterpolatedTimes{1,j})/5);
xlim([0, x_limit])
ylim([0 50])
set(RawDeltaFCs_ax3,'Fontsize',14)
ylim([0 45])
title([SetLabels{2*(j-1)+1}, ' without pv membrane'], 'FontSize', 14)
outpath = [plotdir_nopv, filesep, 'InterpolatedDeltaFCs.png'];
saveas(RawDeltaFCFig3,outpath);

%%
lInt = length(InterpolatedTimes{1,j});
MinTimeToUse = find(InterpolatedTimes{1,j} >= HalfT/2, 1);
MinTimeToUse = 1;
y_temp = NaN(NSets, lInt);
for i = 1:NSets
    y_temp(i,MinTimeToUse:length(InterpolatedDeltaFCs{1,j}{i})) = InterpolatedDeltaFCs{1,j}{i}(MinTimeToUse:end);
end
%%
MatSize = [NSets,(NSets-1)*lInt-2];
LengthY = lInt;
fun = @(x)DeltaFC_Chi2(x,y_temp, MatSize, LengthY);
x0 = zeros(1, 2*(NSets-1));
lb = [-40*ones(1, NSets-1) -2*ones(1, NSets-1)];
ub = -1*lb;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1:(NSets-1);
x_wt= ga(fun,2*(NSets-1), A,b,Aeq,beq,lb,ub, nonlcon, intcon);
x_wts{1,j} = x_wt;

saveVars = [saveVars, 'x_wt'];
WildTypeProfiles{1,j} = {};
MinTime = 0;
%%
wt_Tshifts = zeros(1, NSets);
wt_ProfShifts = zeros(1, NSets);
for i = 1:NSets-1
    wt_Tshifts(i) = x_wt(i);
    wt_ProfShifts(i) = x_wt((NSets-1)+i);
end
AverageDeltaTs(1,j) = mean(wt_Tshifts);
DiffDeltaT = (AverageDeltaTs(1,j))*0.25;
if j == 2
    DiffDeltaT = 0;
end
MinTime = (min(wt_Tshifts))*0.25 ;
MaxTime =  MaxT_wt+(max(wt_Tshifts))*0.25; 
RefTimes = MinTime:0.25:MaxTime;

AverageDeltaProfs(1,j) = mean(wt_ProfShifts);
ShiftedProfiles = NaN(NSets,length(RefTimes));
ShiftedProfileSEs = NaN(NSets,length(RefTimes));

zeroTimeLocation = find(round(RefTimes, 6) == 0);

ShiftedProfiles(1,zeroTimeLocation:zeroTimeLocation + length( InterpolatedDeltaFCs{1,j}{1})-1) =...
    InterpolatedDeltaFCs{1,j}{1}-AverageDeltaProfs(1,j);
ShiftedProfileSEs(1,zeroTimeLocation:zeroTimeLocation + length( InterpolatedDeltaFCs{1,j}{1})-1) =...
    InterpolatedSEs{1,j}{1};
for i = 2:NSets
    firstSetTimeLocation = find(round(RefTimes, 6) == wt_Tshifts(i-1)*0.25);
    ShiftedProfiles(i,firstSetTimeLocation:firstSetTimeLocation + length( InterpolatedDeltaFCs{1,j}{i})-1) =...
        InterpolatedDeltaFCs{1,j}{i}+x_wt(i-1+NSets-1)-AverageDeltaProfs(1,j);
ShiftedProfileSEs(i,firstSetTimeLocation:firstSetTimeLocation + length( InterpolatedDeltaFCs{1,j}{i})-1) =...
   InterpolatedSEs{1,j}{i};

end


WildTypeProfiles{1,j}.Times =RefTimes-DiffDeltaT;

WTProfile = mean(ShiftedProfiles,1,'omitnan');
WildTypeProfiles{1,j}.Profile = WTProfile;
WTProfileSE = std(ShiftedProfiles,1,'omitnan')./sqrt(sum(~isnan(ShiftedProfiles),1));
WildTypeProfiles{1,j}.SE =WTProfileSE;
WTProfileSD = std(ShiftedProfiles,1,'omitnan');
WildTypeProfiles{1,j}.SD =WTProfileSD;
WildTypeProfiles{1,j}.Counts = sum(~isnan(ShiftedProfiles),1);
%WildTypeProfile.Times 

deriv = diff(WTProfile)./diff(WildTypeProfiles{1,j}.Times);
    WTDiff  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    WTProfileSigmaT = WTProfileSE./abs(WTDiff);
    WTProfileSigmaT(WTProfileSigmaT > 10) = NaN;
   WildTypeProfiles{1,j}.SigmaT = WTProfileSigmaT; 
   
   BadTimePoints = WildTypeProfiles{1,j}.Counts  < 2 | WildTypeProfiles{1,j}.Times < 0 | isnan(WildTypeProfiles{1,j}.Profile < 0) ;
   WildTypeProfiles{1,j}.Times = WildTypeProfiles{1,j}.Times(~BadTimePoints);
   WildTypeProfiles{1,j}.Profile = WildTypeProfiles{1,j}.Profile(~BadTimePoints);
   WildTypeProfiles{1,j}.SE = WildTypeProfiles{1,j}.SE(~BadTimePoints);
   WildTypeProfiles{1,j}.SD = WildTypeProfiles{1,j}.SD(~BadTimePoints);
   WildTypeProfiles{1,j}.Counts = WildTypeProfiles{1,j}.Counts(~BadTimePoints);
   WildTypeProfiles{1,j}.SigmaT = WildTypeProfiles{1,j}.SigmaT(~BadTimePoints);
   
   WildTypeProfile = WildTypeProfiles{1,j};
   AverageDeltaT = AverageDeltaTs(1,j)*0.25;
   AverageProfShift = AverageDeltaProfs(1,j);

   saveVars = [saveVars, 'WildTypeProfile','AverageDeltaT', 'AverageProfShift'];
   save([savedir_nopv, filesep, 'MembraneFurrowProfile.mat'], saveVars{:});
%%
%close all
RawDeltaFCFig4 = figure;
set(RawDeltaFCFig4,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax4 = axes(RawDeltaFCFig4);
for i = 1:NSets
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
    if i == 1
  errorbar(WildTypeProfiles{1,j}.Times, ShiftedProfiles(i,~BadTimePoints),ShiftedProfileSEs(i,~BadTimePoints),'.-',...
      'MarkerSize', 10,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    else
        errorbar(WildTypeProfiles{1,j}.Times, ShiftedProfiles(i,~BadTimePoints),...
            ShiftedProfileSEs(i,~BadTimePoints),'.-',...
      'MarkerSize', 10,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    end
 
 
 hold on 
end





   errorbar(WildTypeProfiles{1,j}.Times, WildTypeProfiles{1,j}.Profile,WildTypeProfiles{1,j}.SE,'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  
  grid on
  hold off
  
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j)) 'Average Profile'];
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 80])
set(RawDeltaFCs_ax4,'Fontsize',14)
ylim([0 45])
title([SetLabels{2*(j-1)+1}, ' without pv membrane'], 'FontSize', 14)
outpath = [plotdir_nopv, filesep, 'ShiftedDeltaFCs.png'];
saveas(RawDeltaFCFig4,outpath);

%%
close all
RawDeltaFCFig5 = figure;
set(RawDeltaFCFig5,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax5 = axes(RawDeltaFCFig5);

   scatter(WildTypeProfiles{1,j}.Times, WildTypeProfiles{1,j}.SE, 20,'filled',...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k');%,'LineStyle', 'none')
  






  grid on
  hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j)) 'Average Profile'];

  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('$\sigma_{\delta_{FC}}\,\, (\mu m)$','interpreter','latex')
xlim([0, 80])
set(RawDeltaFCs_ax5,'Fontsize',14)
%ylim([0 45])
title([SetLabels{2*(j-1)+1}, ' without pv membrane'], 'FontSize', 14)
outpath = [plotdir_nopv, filesep, 'StdErrDeltaFCs.png'];
saveas(RawDeltaFCFig5,outpath);

%%
close all
RawDeltaFCFig5 = figure;
set(RawDeltaFCFig5,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax5 = axes(RawDeltaFCFig5);

   scatter(WildTypeProfiles{1,j}.Times, WildTypeProfiles{1,j}.SigmaT, 20,'filled',...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k');%,'LineStyle', 'none')
  






  grid on
  hold off
  
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j)) 'Average Profile'];

  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);



xlabel('Time into cycle 14 (minutes)')
ylabel('$\sigma_{t}\,\, (m)$','interpreter','latex')
xlim([0, 80])
set(RawDeltaFCs_ax5,'Fontsize',14)
%ylim([0 45])
title([SetLabels{2*(j-1)+1}, ' without pv membrane'], 'FontSize', 14)
outpath = [plotdir_nopv, filesep, 'SigmaTDeltaFCs.png'];
saveas(RawDeltaFCFig5,outpath);


%% Start with data that includes perivitelline membrane
saveVars = {};


%%
if ~isdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements')
mkdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements');
end
plotdir = ['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/Figures',filesep, DirLabels{(j)*2}];
if ~isdir(plotdir)
mkdir(plotdir);
end
savedir = ['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements',filesep, DirLabels{(j)*2}];
if ~isdir(savedir)
mkdir(savedir);
end
%%

close all
RawDeltaFCFig = figure(1);
set(RawDeltaFCFig,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax = axes(RawDeltaFCFig);
RawDeltaFCs_ax.FontSize = 14;
for i = 1:NSets
    plot(NC14FrameTimes{2,j}{i}, DeltaFCs{2, j}{i}, '-o', 'Color', colors(i,:), 'MarkerSize', 5,...
         'MarkerFaceColor', colors(i,:));
hold on 
end
grid on 
hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j))];
hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)', 'FontSize', 14)
ylabel('Membrane Furrow Depth (microns)', 'FontSize', 14)
%xlim([0, 65])
ylim([0 50])
title([SetLabels{2*(j)}], 'FontSize', 14)
set(RawDeltaFCs_ax,'Fontsize',14)
outpath = [plotdir, filesep, 'RawDeltaFCs_noErrorBar.png'];
saveas(RawDeltaFCFig,outpath);

%%
close all
RawDeltaFCFig2 = figure(2);
set(RawDeltaFCFig2,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax2 = axes(RawDeltaFCFig2);
RawDeltaFCs_ax2.FontSize = 14;
for i = 1:NSets
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
  errorbar(NC14FrameTimes{2,j}{i}, DeltaFCs{2,j}{i},SEDeltaFCs{2,j}{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
 
 
 hold on 
end
grid on 
hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j))];
hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)', 'FontSize', 16)
ylabel('Membrane Furrow Depth (microns)', 'FontSize', 16)
%xlim([0, 65])
ylim([0 50])
title([SetLabels{2*(j)}], 'FontSize', 14)
set(RawDeltaFCs_ax2,'Fontsize',14)
outpath = [plotdir, filesep, 'RawDeltaFCs.png'];
saveas(RawDeltaFCFig2,outpath);
%%
MaxDeltaShift =1;
MaxTShift = 10*4;
LengthNC14s = zeros(1, NSets);
LengthNC14FrameTimes = zeros(1,NSets);
MaxDeltaMeasureFrames = zeros(1,NSets);
MaxDeltaMeasureTimes = zeros(1,NSets);
for i = 1:NSets
    GoodMeasurements = ~isnan(DeltaFCs{2,j}{i});
    DeltaFCs{2,j}{i} = DeltaFCs{2,j}{i}(GoodMeasurements);
    NC14FrameTimes{2,j}{i} = NC14FrameTimes{2,j}{i}(GoodMeasurements);
    LengthNC14s(i) = length(DeltaFCs{2,j}{i});
    LengthNC14FrameTimes(i) = max(NC14FrameTimes{2,j}{i});

   
end
   
l = max(LengthNC14s);



MaxT_wt = floor(max(LengthNC14FrameTimes));
HalfT = MaxT_wt/2;
InterpolatedTimes{2,j} = 0:0.25:MaxT_wt;


InterpolatedDeltaFCs{2,j} = cell(1,NSets);
InterpolatedSEs{2,j} = cell(1,NSets);
Interpolated_dDelta_dts{2,j} =cell(1,NSets);
Interpolated_sigma_ts{2,j} = cell(1,NSets);

for i = 1:NSets
    InterpolatedDeltaFCs{2,j}{i} = interp1(NC14FrameTimes{2,j}{i},DeltaFCs{2,j}{i}, InterpolatedTimes{2,j});
    InterpolatedSEs{2,j}{i} = NaN(1,length(InterpolatedDeltaFCs{2,j}{i}));
    for k = 1:length(InterpolatedTimes{2,j})
        t = InterpolatedTimes{2,j}(k);
        Tindex = find(round(NC14FrameTimes{2,j}{i},5) == t);
        if isempty(Tindex)
            Tlow = find(NC14FrameTimes{2,j}{i} < t, 1, 'last');
            Thigh = find(NC14FrameTimes{2,j}{i} > t, 1);
            if ~isempty(Tlow) & ~isempty(Thigh)
                t1 = NC14FrameTimes{2,j}{i}(Tlow);
                s1  = SEDeltaFCs{2,j}{i}(Tlow);
                t2 = NC14FrameTimes{2,j}{i}(Thigh);
                s2  = SEDeltaFCs{2,j}{i}(Thigh);
                dfdt1 = (t2-t)/(t2-t1);
                dfdt2 = (t-t1)/(t2-t1);
                InterpolatedSEs{2,j}{i}(k) = sqrt(dfdt1^2*s1^2+dfdt2^2*s2^2);
            end
        else 
           InterpolatedSEs{2,j}{i}(k) =  SEDeltaFCs{2,j}{i}(Tindex);
        end
    end
    
    deriv = diff(InterpolatedDeltaFCs{2,j}{i})./diff(InterpolatedTimes{2,j});
    Interpolated_dDelta_dts{2,j}{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    Interpolated_sigma_ts{2,j}{i} = InterpolatedSEs{2,j}{i}./abs(Interpolated_dDelta_dts{2,j}{i});
    Interpolated_sigma_ts{2,j}{i}(Interpolated_sigma_ts{2,j}{i} > 10) = NaN;
    
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
  errorbar(InterpolatedTimes{2,j}, InterpolatedDeltaFCs{2,j}{i},InterpolatedSEs{2,j}{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
 
 
 hold on 
end
grid on 
hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j))];
hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
x_limit = 5*ceil(max(InterpolatedTimes{2,j})/5);
xlim([0, x_limit])
ylim([0 50])
set(RawDeltaFCs_ax3,'Fontsize',14)
ylim([0 50])
title([SetLabels{2*(j)}], 'FontSize', 14)
outpath = [plotdir, filesep, 'InterpolatedDeltaFCs.png'];
saveas(RawDeltaFCFig3,outpath);

%%
lInt = length(InterpolatedTimes{2,j});
MinTimeToUse = find(InterpolatedTimes{2,j} >= HalfT/2, 1);
MinTimeToUse = 1;
y_temp = NaN(NSets, lInt);
for i = 1:NSets
    y_temp(i,MinTimeToUse:length(InterpolatedDeltaFCs{2,j}{i})) = InterpolatedDeltaFCs{2,j}{i}(MinTimeToUse:end);
end
%%
MatSize = [NSets,(NSets-1)*lInt-2];
LengthY = lInt;
fun = @(x)DeltaFC_Chi2(x,y_temp, MatSize, LengthY);
x0 = zeros(1, 2*(NSets-1));
lb = [-40*ones(1, NSets-1) -2*ones(1, NSets-1)];
ub = -1*lb;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1:(NSets-1);
x_wt= ga(fun,2*(NSets-1), A,b,Aeq,beq,lb,ub, nonlcon, intcon);
x_wts{2,j} = x_wt;

saveVars = [saveVars, 'x_wt'];
WildTypeProfiles{2,j} = {};
MinTime = 0;
%%
wt_Tshifts = zeros(1, NSets);
wt_ProfShifts = zeros(1, NSets);
for i = 1:NSets-1
    wt_Tshifts(i) = x_wt(i);
    wt_ProfShifts(i) = x_wt((NSets-1)+i);
end
AverageDeltaTs(2,j) = mean(wt_Tshifts);
DiffDeltaT = (AverageDeltaTs(2,j))*0.25;
if j == 2
    DiffDeltaT = 0;
end
MinTime = (min(wt_Tshifts))*0.25 ;
MaxTime =  MaxT_wt+(max(wt_Tshifts))*0.25; 
RefTimes = MinTime:0.25:MaxTime;

AverageDeltaProfs(2,j) = mean(wt_ProfShifts);
ShiftedProfiles = NaN(NSets,length(RefTimes));
ShiftedProfileSEs = NaN(NSets,length(RefTimes));

zeroTimeLocation = find(round(RefTimes, 6) == 0);

ShiftedProfiles(1,zeroTimeLocation:zeroTimeLocation + length( InterpolatedDeltaFCs{2,j}{1})-1) =...
    InterpolatedDeltaFCs{2,j}{1}-AverageDeltaProfs(2,j);
ShiftedProfileSEs(1,zeroTimeLocation:zeroTimeLocation + length( InterpolatedDeltaFCs{2,j}{1})-1) =...
    InterpolatedSEs{2,j}{1};
for i = 2:NSets
    firstSetTimeLocation = find(round(RefTimes, 6) == wt_Tshifts(i-1)*0.25);
    ShiftedProfiles(i,firstSetTimeLocation:firstSetTimeLocation + length( InterpolatedDeltaFCs{2,j}{i})-1) =...
        InterpolatedDeltaFCs{2,j}{i}+x_wt(i-1+NSets-1)-AverageDeltaProfs(2,j);
ShiftedProfileSEs(i,firstSetTimeLocation:firstSetTimeLocation + length( InterpolatedDeltaFCs{2,j}{i})-1) =...
   InterpolatedSEs{2,j}{i};

end


WildTypeProfiles{2,j}.Times =RefTimes-DiffDeltaT;

WTProfile = mean(ShiftedProfiles,1,'omitnan');
WildTypeProfiles{2,j}.Profile = WTProfile;
WTProfileSE = std(ShiftedProfiles,1,'omitnan')./sqrt(sum(~isnan(ShiftedProfiles),1));
WildTypeProfiles{2,j}.SE =WTProfileSE;
WTProfileSD = std(ShiftedProfiles,1,'omitnan');
WildTypeProfiles{2,j}.SD =WTProfileSD;
WildTypeProfiles{2,j}.Counts = sum(~isnan(ShiftedProfiles),1);
%WildTypeProfile.Times 

deriv = diff(WTProfile)./diff(WildTypeProfiles{2,j}.Times);
    WTDiff  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    WTProfileSigmaT = WTProfileSE./abs(WTDiff);
    WTProfileSigmaT(WTProfileSigmaT > 10) = NaN;
   WildTypeProfiles{2,j}.SigmaT = WTProfileSigmaT; 
   
   
   BadTimePoints = WildTypeProfiles{2,j}.Counts  < 2 | WildTypeProfiles{2,j}.Times < 0 | isnan(WildTypeProfiles{2,j}.Profile < 0) ;
   WildTypeProfiles{2,j}.Times = WildTypeProfiles{2,j}.Times(~BadTimePoints);
   WildTypeProfiles{2,j}.Profile = WildTypeProfiles{2,j}.Profile(~BadTimePoints);
   WildTypeProfiles{2,j}.SE = WildTypeProfiles{2,j}.SE(~BadTimePoints);
   WildTypeProfiles{2,j}.SD = WildTypeProfiles{2,j}.SD(~BadTimePoints);
   WildTypeProfiles{2,j}.Counts = WildTypeProfiles{2,j}.Counts(~BadTimePoints);
   WildTypeProfiles{2,j}.SigmaT = WildTypeProfiles{2,j}.SigmaT(~BadTimePoints);
   
   WildTypeProfile = WildTypeProfiles{2,j};
   AverageDeltaT = AverageDeltaTs(2,j)*0.25;
   AverageProfShift = AverageDeltaProfs(2,j);

   saveVars = [saveVars, 'WildTypeProfile','AverageDeltaT', 'AverageProfShift'];
   save([savedir, filesep, 'MembraneFurrowProfile.mat'], saveVars{:});
%%
%close all
RawDeltaFCFig4 = figure;
set(RawDeltaFCFig4,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax4 = axes(RawDeltaFCFig4);
for i = 1:NSets
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
    if i == 1
  errorbar(WildTypeProfiles{2,j}.Times, ShiftedProfiles(i,~BadTimePoints),ShiftedProfileSEs(i,~BadTimePoints),'.-',...
      'MarkerSize', 10,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    else
        errorbar(WildTypeProfiles{2,j}.Times, ShiftedProfiles(i,~BadTimePoints),...
            ShiftedProfileSEs(i,~BadTimePoints),'.-',...
      'MarkerSize', 10,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    end
 
 
 
end





   errorbar(WildTypeProfiles{2,j}.Times, WildTypeProfiles{2,j}.Profile, WildTypeProfiles{2,j}.SE,'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  
  grid on
  hold off
  
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j)) 'Average Profile'];
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 80])
set(RawDeltaFCs_ax4,'Fontsize',14)
ylim([0 50])
title([SetLabels{2*(j)}], 'FontSize', 14)
outpath = [plotdir, filesep, 'ShiftedDeltaFCs.png'];
saveas(RawDeltaFCFig4,outpath);

%%
close all
RawDeltaFCFig5 = figure;
set(RawDeltaFCFig5,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax5 = axes(RawDeltaFCFig5);

   scatter(WildTypeProfiles{2,j}.Times, WildTypeProfiles{2,j}.SE, 20,'filled',...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k');%,'LineStyle', 'none')
  






  grid on
  hold off
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j)) 'Average Profile'];


  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('$\sigma_{\delta_{FC}}\,\, (\mu m)$','interpreter','latex')
xlim([0, 80])
set(RawDeltaFCs_ax5,'Fontsize',14)
title([SetLabels{2*(j)}], 'FontSize', 14)
outpath = [plotdir, filesep, 'StdErrDeltaFCs.png'];
saveas(RawDeltaFCFig5,outpath);

%%
close all
RawDeltaFCFig5 = figure;
set(RawDeltaFCFig5,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax5 = axes(RawDeltaFCFig5);

   scatter(WildTypeProfiles{2,j}.Times, WildTypeProfiles{2,j}.SigmaT, 20,'filled',...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k');%,'LineStyle', 'none')
  






  grid on
  hold off
  
  Plot_Labels = [Prefixes(GroupStartIndices(j):GroupEndIndices(j)) 'Average Profile'];

  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);



xlabel('Time into cycle 14 (minutes)')
ylabel('$\sigma_{t}\,\, (m)$','interpreter','latex')
xlim([0, 80])
set(RawDeltaFCs_ax5,'Fontsize',14)
title([SetLabels{2*(j)}], 'FontSize', 14)
outpath = [plotdir, filesep, 'SigmaTDeltaFCs.png'];
saveas(RawDeltaFCFig5,outpath);


%%
save([savedir, filesep, 'MembraneFurrowProfile.mat'], saveVars{:});
end


%%
if ~isdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements')
mkdir('S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements');
end
Labels = cell(2, 4);
Labels{1,1} = SetLabels{1};
Labels{2,1} = SetLabels{2};
Labels{1,2} = SetLabels{3};
Labels{2,2} = SetLabels{4};
Labels{1,3} = SetLabels{5};
Labels{2,3} = SetLabels{6};
Labels{1,4} = SetLabels{7};
Labels{2,4} = SetLabels{8};
saveVars = {'Labels', 'WildTypeProfiles', 'DubuisTimes', 'DubuisMeanProfile','DubuisFixedMeanProfile', 'DubuisFurrowData', 'x_wts', 'AverageDeltaTs', 'AverageProfShift'};
save(['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements', filesep, 'All25CMembraneFurrowProfile.mat'], saveVars{:});
compplotdir = ['S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/Figures',filesep, 'T25C_CompPlots'];
if ~isdir(compplotdir)
mkdir(compplotdir);
end
colors = brewermap(min(NSets, 4),'Spectral');
close all
RawDeltaFCFig4 = figure;
set(RawDeltaFCFig4,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax4 = axes(RawDeltaFCFig4);
for i = 1:4
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
    if i == 1
  errorbar( WildTypeProfiles{1,i}.Times, WildTypeProfiles{1,i}.Profile,...
      WildTypeProfiles{1,i}.SE,'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    else
        errorbar(WildTypeProfiles{1,i}.Times,...
            WildTypeProfiles{1,i}.Profile,...
            WildTypeProfiles{1,i}.SE,'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    end
 
 
 hold on 
end





   errorbar(DubuisTimes, DubuisMeanProfile,zeros(1, length(DubuisMeanProfile)),'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  
  grid on
  hold off
  
  Plot_Labels = {'Squished yw 25ºC','All HisRFP 25ºC', 'Squished HisRFP 25ºC', 'Unsquished HisRFP 25ºC', 'Dubuis Profile'};
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 80])
set(RawDeltaFCs_ax4,'Fontsize',14)
ylim([0 50])
title([SetLabels{2*(j)}], 'FontSize', 14)
outpath = [compplotdir, filesep, 'CompT25CShiftedDeltaFCs_NoPV.png'];
saveas(RawDeltaFCFig4,outpath);

%%
close all
RawDeltaFCFig4 = figure;
set(RawDeltaFCFig4,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax4 = axes(RawDeltaFCFig4);
for i = 1:4
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
    if i == 1
  errorbar( WildTypeProfiles{2,i}.Times, WildTypeProfiles{2,i}.Profile,...
      WildTypeProfiles{2,i}.SE,'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    else
        errorbar(WildTypeProfiles{2,i}.Times,...
            WildTypeProfiles{2,i}.Profile,...
            WildTypeProfiles{2,i}.SE,'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    end
 
 
 hold on 
end





   errorbar(DubuisTimes, DubuisMeanProfile,zeros(1, length(DubuisMeanProfile)),'.-','MarkerSize', 20,...
      'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
      'Color', 'k');%,'LineStyle', 'none')
  
  grid on
  hold off
  
  Plot_Labels = {'Squished yw 25ºC','All HisRFP 25ºC', 'Squished HisRFP 25ºC', 'Unsquished HisRFP 25ºC', 'Dubuis Profile'};
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 80])
set(RawDeltaFCs_ax4,'Fontsize',14)
ylim([0 50])
title([SetLabels{2*(j)}], 'FontSize', 14)
outpath = [compplotdir, filesep, 'CompT25CShiftedDeltaFCs.png'];
saveas(RawDeltaFCFig4,outpath);

%%

AllProfiles = WildTypeProfiles;
MaxDeltaShift =1;
MaxTShift = 40;
LengthNC14s = zeros(1, 5);
LengthNC14FrameTimes = zeros(1,5);
MaxDeltaMeasureFrames = zeros(1,5);
MaxDeltaMeasureTimes = zeros(1,5);
for i = 1:5
    if i < 5
    GoodMeasurements = ~isnan(AllProfiles{1,i}.Profile) & (AllProfiles{1,i}.Times >= 0);
    AllProfiles{1,i}.Profile = AllProfiles{1,i}.Profile(GoodMeasurements);
    AllProfiles{1,i}.Times = AllProfiles{1,i}.Times(GoodMeasurements);
    AllProfiles{1,i}.SE = AllProfiles{1,i}.SE(GoodMeasurements);
    AllProfiles{1,i}.SD = AllProfiles{1,i}.SD(GoodMeasurements);
    AllProfiles{1,i}.SigmaT = AllProfiles{1,i}.SigmaT(GoodMeasurements);
    LengthNC14s(i) = length(AllProfiles{1,i}.Profile);
    LengthNC14FrameTimes(i) = max(AllProfiles{1,i}.Times);
    else
        GoodMeasurements = ~isnan(DubuisMeanProfile) & (DubuisTimes > 0) ;
    AllProfiles{1,i}.Profile= DubuisMeanProfile(GoodMeasurements);
    AllProfiles{1,i}.Times = DubuisTimes(GoodMeasurements);
    AllProfiles{1,i}.SE = NaN(1, length(AllProfiles{1,i}.Times));
    AllProfiles{1,i}.SD = NaN(1, length(AllProfiles{1,i}.Times));
    AllProfiles{1,i}.SigmaT = NaN(1, length(AllProfiles{1,i}.Times));
    LengthNC14s(i) = length(DubuisMeanProfile);
    LengthNC14FrameTimes(i) = max(DubuisTimes);
    end
   
end


   
l = max(LengthNC14s);


MaxT_wt = floor(max(LengthNC14FrameTimes));
HalfT = MaxT_wt/2;
AllTimes = 0:1:MaxT_wt;


AllInterpolatedDeltaFCs = cell(1,4);
AllInterpolatedSEs = cell(1,4);
AllInterpolated_dDelta_dts =cell(1,4);
AllInterpolated_sigma_ts = cell(1,4);


for i = 1:5
    if i < 5
    AllInterpolatedDeltaFCs{i} = interp1(AllProfiles{1,i}.Times,AllProfiles{1,i}.Profile, AllTimes);
    AllInterpolatedSEs{i} = NaN(1,length(AllInterpolatedDeltaFCs{i}));
    for k = 1:length(AllTimes)
        t = AllTimes(k);
        Tindex = find(round(AllProfiles{1,i}.Times,5) == t);
        if isempty(Tindex)
            Tlow = find(AllProfiles{1,i}.Times < t, 1, 'last');
            Thigh = find(AllProfiles{1,i}.Times > t, 1);
            if ~isempty(Tlow) & ~isempty(Thigh)
                t1 = AllProfiles{1,i}.Times(Tlow);
                s1  = AllProfiles{1,i}.SE(Tlow);
                t2 = AllProfiles{1,i}.Times(Thigh);
                s2  = AllProfiles{1,i}.SE(Thigh);
                dfdt1 = (t2-t)/(t2-t1);
                dfdt2 = (t-t1)/(t2-t1);
                AllInterpolatedSEs{i}(k) = sqrt(dfdt1^2*s1^2+dfdt2^2*s2^2);
            end
        else 
           AllInterpolatedSEs{i}(k) =  AllProfiles{1,i}.SE(Tindex);
        end
    end
    
    deriv = diff(AllInterpolatedDeltaFCs{i})./diff(AllTimes);
    AllInterpolated_dDelta_dts{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    AllInterpolated_sigma_ts{i} = AllInterpolatedSEs{i}./abs(AllInterpolated_dDelta_dts{i});
    AllInterpolated_sigma_ts{i}(AllInterpolated_sigma_ts{i} > 10) = NaN;
    else
        AllInterpolatedDeltaFCs{i} = interp1(AllProfiles{1,i}.Times,AllProfiles{1,i}.Profile, AllTimes);
    AllInterpolatedSEs{i} = NaN(1,length(AllInterpolatedDeltaFCs{i}));
    for k = 1:length(AllTimes)
        t = AllTimes(k);
        Tindex = find(round(AllProfiles{1,i}.Times,5) == t);
        if isempty(Tindex)
            Tlow = find(AllProfiles{1,i}.Times < t, 1, 'last');
            Thigh = find(AllProfiles{1,i}.Times > t, 1);
            if ~isempty(Tlow) & ~isempty(Thigh)
                t1 = AllProfiles{1,i}.Times(Tlow);
                s1  = AllProfiles{1,i}.SE(Tlow);
                t2 = AllProfiles{1,i}.Times(Thigh);
                s2  = AllProfiles{1,i}.SE(Thigh);
                dfdt1 = (t2-t)/(t2-t1);
                dfdt2 = (t-t1)/(t2-t1);
                AllInterpolatedSEs{i}(k) = sqrt(dfdt1^2*s1^2+dfdt2^2*s2^2);
            end
        else 
           AllInterpolatedSEs{i}(k) =  AllProfiles{1,i}.SE(Tindex);
        end
    end
    
    deriv = diff(AllInterpolatedDeltaFCs{i})./diff(AllTimes);
    AllInterpolated_dDelta_dts{i}  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    AllInterpolated_sigma_ts{i} = AllInterpolatedSEs{i}./abs(AllInterpolated_dDelta_dts{i});
    AllInterpolated_sigma_ts{i}(AllInterpolated_sigma_ts{i} > 10) = NaN;
    end
    
end
%%
close all
RawDeltaFCFig3 = figure(3);
set(RawDeltaFCFig3,'units', 'normalized', 'position',[0.05, 0.05, 0.8, 0.5]);
    set(gcf,'color','w');
RawDeltaFCs_ax3 = axes(RawDeltaFCFig3);
for i = 1:5
%  errorbar(NC14FrameTimes{3}, DeltaFCs{3},SEDeltaFCs{3},SEDeltaFCs{3},...
%      sigma_ts{3}, sigma_ts{3},'MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
if i < 5
  errorbar(AllTimes, AllInterpolatedDeltaFCs{i},AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
else
    errorbar(AllTimes, AllInterpolatedDeltaFCs{i},AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor','k',...
      'Color', 'k', 'LineWidth', 1.5);%,'LineStyle', 'none')
end
 
 
 hold on 
end
grid on 
hold off
  Plot_Labels = {'Squished yw 25ºC','All HisRFP 25ºC', 'Squished HisRFP 25ºC', 'Unsquished HisRFP 25ºC', 'Dubuis Profile'};
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);


xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
x_limit = 5*ceil(max(InterpolatedTimes{2,j})/5);
xlim([0, 80])
ylim([0 50])
set(RawDeltaFCs_ax3,'Fontsize',14)
ylim([0 50])
title('All InterpolatedProfiles', 'FontSize', 14)
outpath = [compplotdir, filesep, 'CompPlot_InterpolatedDeltaFCs.png'];
saveas(RawDeltaFCFig3,outpath);

%%
lInt = length(AllTimes);
MinTimeToUse = find(AllTimes >= HalfT, 1);
MinTimeToUse = 1;
y_temp = NaN(4, lInt);
for i = 1:5
    y_temp(i,MinTimeToUse:length(AllInterpolatedDeltaFCs{i})) = AllInterpolatedDeltaFCs{i}(MinTimeToUse:end);
end
%%
NSets = 5;
MatSize = [NSets,(NSets-1)*lInt-2];
LengthY = lInt;
fun = @(x)DeltaFC_Chi2(x,y_temp, MatSize, LengthY);
x0 = zeros(1, 2*(NSets-1));
lb = [-15*ones(1, NSets-1) -2*ones(1, NSets-1)];
ub = -1*lb;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1:(NSets-1);
x_wt= ga(fun,2*(NSets-1), A,b,Aeq,beq,lb,ub, nonlcon, intcon);

MinTime = 0;
%%
CompProfiles = {};

wt_Tshifts = zeros(1, NSets);
wt_ProfShifts = zeros(1, NSets);
for i = 1:NSets-1
    wt_Tshifts(i) = x_wt(i);
    wt_ProfShifts(i) = x_wt((NSets-1)+i);
end
AverageDeltaT = round(mean(wt_Tshifts));
StartMinDeltaT = min([wt_Tshifts 0]);
DiffDeltaT = AverageDeltaT - StartMinDeltaT;
MaxTime =  MaxT_wt+(max(wt_Tshifts)-min(wt_Tshifts)); 
CompProfiles.Times = 0:MaxTime;
RefTimes = AllTimes-min([wt_Tshifts 0]);

AverageDeltaProf = mean(wt_ProfShifts);
ShiftedProfiles = NaN(NSets,length(CompProfiles.Times));
ShiftedProfileSEs = NaN(NSets,length(CompProfiles.Times));
ShiftedProfiles(1,RefTimes+1) = AllInterpolatedDeltaFCs{1}-AverageDeltaProf;
ShiftedProfileSEs(1,RefTimes+1) = AllInterpolatedSEs{1};
for i = 2:NSets
    ShiftedProfiles(i,RefTimes+1+x_wt(i-1)) = AllInterpolatedDeltaFCs{i}+x_wt(i-1+NSets-1)-AverageDeltaProf;
    ShiftedProfileSEs(i,RefTimes+1+x_wt(i-1)) = AllInterpolatedSEs{i};
end


CompProfiles.Times = CompProfiles.Times-DiffDeltaT;

WTProfile = mean(ShiftedProfiles,1,'omitnan');
CompProfiles.Profile = WTProfile;
WTProfileSE = std(ShiftedProfiles,1,'omitnan')./sqrt(sum(~isnan(ShiftedProfiles),1));
CompProfiles.SE =WTProfileSE;
WTProfileSD = std(ShiftedProfiles,1,'omitnan');
CompProfiles.SD =WTProfileSD;
%WildTypeProfile.Times 

deriv = diff(WTProfile)./diff(CompProfiles.Times);
    WTDiff  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    WTProfileSigmaT = WTProfileSE./abs(WTDiff);
    WTProfileSigmaT(WTProfileSigmaT > 10) = NaN;
   CompProfiles.SigmaT = WTProfileSigmaT; 
   WildTypeProfile =CompProfiles;
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
  errorbar(AllTimes-min(wt_Tshifts)-DiffDeltaT, AllInterpolatedDeltaFCs{i}-AverageDeltaProf,...
      AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    elseif i < 5
        errorbar(AllTimes-min(wt_Tshifts)-DiffDeltaT+x_wt(i-1),...
            AllInterpolatedDeltaFCs{i}+x_wt(i+NSets-2)-AverageDeltaProf,...
            AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    else
        errorbar(AllTimes-min(wt_Tshifts)-DiffDeltaT+x_wt(i-1),...
            AllInterpolatedDeltaFCs{i}+x_wt(i+NSets-2)-AverageDeltaProf,...
            AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor','k',...
      'Color', 'k', 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    end
 
 
 hold on 
end



% 
% 
%    errorbar(WildTypeProfiles{2,j}.Times, WTProfile,WTProfileSE,'.-','MarkerSize', 20,...
%       'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
%       'Color', 'k');%,'LineStyle', 'none')
  
  grid on
  hold off
  
  Plot_Labels = {'Squished yw 25ºC', 'All HisRFP 25ºC','Squished HisRFP 25ºC', 'Unsquished HisRFP 25ºC', 'Dubuis Profile'};
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 80])
set(RawDeltaFCs_ax4,'Fontsize',14)
ylim([0 50])
title('Shifted Mean Profiles', 'FontSize', 14)
outpath = [compplotdir, filesep, 'ShiftedMeanDeltaFCs.png'];
saveas(RawDeltaFCFig4,outpath);

%%

close all
EmbryoSizeFig = figure;
set(EmbryoSizeFig,'units', 'normalized', 'position',[0.05, 0.05, 0.6, 0.5]);
set(gcf,'color','w');
EmbryoSizeAx = axes(EmbryoSizeFig);
colors2 = brewermap(11,'Spectral');
colors2 = colors2([1 3 8 9 11], :);
for i = 1:3

     scatter(APLengths{i}, DVLengths{i}, 100, 'o', 'MarkerFaceColor', colors2(i,:), 'MarkerEdgeColor', 'k');
     hold on 
end
grid on 
hold off

  
  Plot_Labels = {'Squished yw 25ºC', 'Squished HisRFP 25ºC', 'Unsquished HisRFP 25ºC'};
  hlegend = legend(Plot_Labels, 'Location', 'northeast',...
            'FontSize', 14);

xlabel('AP Length (\mum)')
ylabel('DV Length (\mum)')
xlim([500, 600])
set(EmbryoSizeAx,'Fontsize',14)
ylim([150 250])
title('Embryo Dimensions for 25ºC Membrane Movies', 'FontSize', 16)
outpath = [compplotdir, filesep, 'EmbryoDimensions.png'];
saveas(EmbryoSizeFig,outpath);

%%

lInt = length(AllTimes);
MinTimeToUse = find(AllTimes >= HalfT, 1);
MinTimeToUse = 1;
y_temp = NaN(5, lInt);
for i = 1:5
    y_temp(i,MinTimeToUse:length(AllInterpolatedDeltaFCs{i})) = AllInterpolatedDeltaFCs{i}(MinTimeToUse:end);
end
%%
NSets = 5;
MatSize = [NSets,(NSets-1)*lInt-2];
LengthY = lInt;
fun = @(x)DeltaFC_Chi2(x,y_temp, MatSize, LengthY);
x0 = zeros(1, 2*(NSets-1));
lb = [-20*ones(1, NSets-1) -2*ones(1, NSets-1)];
ub = -1*lb;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1:(NSets-1);
x_wt= ga(fun,2*(NSets-1), A,b,Aeq,beq,lb,ub, nonlcon, intcon);

MinTime = 0;
%%
CompProfiles = {};

wt_Tshifts = zeros(1, NSets);
wt_ProfShifts = zeros(1, NSets);
for i = 1:NSets-1
    wt_Tshifts(i) = x_wt(i);
    wt_ProfShifts(i) = x_wt((NSets-1)+i);
end
AverageDeltaT = round(mean(wt_Tshifts));
StartMinDeltaT = min([wt_Tshifts 0]);
DiffDeltaT = AverageDeltaT - StartMinDeltaT;
MaxTime =  MaxT_wt+(max(wt_Tshifts)-min(wt_Tshifts)); 
CompProfiles.Times = 0:MaxTime;
RefTimes = AllTimes-min([wt_Tshifts 0]);

AverageDeltaProf = mean(wt_ProfShifts);
ShiftedProfiles = NaN(NSets,length(CompProfiles.Times));
ShiftedProfileSEs = NaN(NSets,length(CompProfiles.Times));
ShiftedProfiles(1,RefTimes+1) = AllInterpolatedDeltaFCs{1}-AverageDeltaProf;
ShiftedProfileSEs(1,RefTimes+1) = AllInterpolatedSEs{1};
for i = 2:NSets
    ShiftedProfiles(i,RefTimes+1+x_wt(i-1)) = AllInterpolatedDeltaFCs{i}+x_wt(i-1+NSets-1)-AverageDeltaProf;
    ShiftedProfileSEs(i,RefTimes+1+x_wt(i-1)) = AllInterpolatedSEs{i};
end


CompProfiles.Times = CompProfiles.Times-DiffDeltaT;

WTProfile = mean(ShiftedProfiles,1,'omitnan');
CompProfiles.Profile = WTProfile;
WTProfileSE = std(ShiftedProfiles,1,'omitnan')./sqrt(sum(~isnan(ShiftedProfiles),1));
CompProfiles.SE =WTProfileSE;
WTProfileSD = std(ShiftedProfiles,1,'omitnan');
CompProfiles.SD =WTProfileSD;
%WildTypeProfile.Times 

deriv = diff(WTProfile)./diff(CompProfiles.Times);
    WTDiff  = [deriv(1), (deriv(1:end-1)+deriv(2:end))/2, deriv(end)];
    WTProfileSigmaT = WTProfileSE./abs(WTDiff);
    WTProfileSigmaT(WTProfileSigmaT > 10) = NaN;
   CompProfiles.SigmaT = WTProfileSigmaT; 
   WildTypeProfile =CompProfiles;
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
  errorbar(AllTimes-min(wt_Tshifts)-DiffDeltaT, AllInterpolatedDeltaFCs{i}-AverageDeltaProf,...
      AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    elseif i < 5
        errorbar(AllTimes-min(wt_Tshifts)-DiffDeltaT+x_wt(i-1),...
            AllInterpolatedDeltaFCs{i}+x_wt(i+NSets-2)-AverageDeltaProf,...
            AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor',colors(i,:),...
      'Color', colors(i,:), 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    else
        errorbar(AllTimes-min(wt_Tshifts)-DiffDeltaT+x_wt(i-1),...
            AllInterpolatedDeltaFCs{i}+x_wt(i+NSets-2)-AverageDeltaProf,...
            AllInterpolatedSEs{i},'.-',...
      'MarkerSize', 20,...%'MarkerEdgeColor', 'k',...
      'MarkerFaceColor','k',...
      'Color', 'k', 'LineWidth', 1.5);%,'LineStyle', 'none')
  hold on 
    end
 
 
 hold on 
end



% 
% 
%    errorbar(WildTypeProfiles{2,j}.Times, WTProfile,WTProfileSE,'.-','MarkerSize', 20,...
%       'MarkerEdgeColor', 'k','MarkerFaceColor','k',...
%       'Color', 'k');%,'LineStyle', 'none')
  
  grid on
  hold off
  
  Plot_Labels = {'Squished yw 25ºC', 'All HisRFP 25ºC','Squished HisRFP 25ºC', 'Unsquished HisRFP 25ºC', 'Dubuis Profile'};
  hlegend = legend(Plot_Labels, 'Location', 'eastoutside',...
            'FontSize', 12);

xlabel('Time into cycle 14 (minutes)')
ylabel('Membrane Furrow Depth (microns)')
xlim([0, 70])
set(RawDeltaFCs_ax4,'Fontsize',14)
ylim([0 50])
title('Shifted Mean Profiles', 'FontSize', 14)
outpath = [compplotdir, filesep, 'ShiftedMeanDeltaFCs.png'];
saveas(RawDeltaFCFig4,outpath);

