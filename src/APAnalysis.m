function APAnalysis(dataset, varargin)
%APAnalysis(dataset, options)
%
%DESCRIPTION
%Performs analyses on AP profiles generated from time traces of enhancer
%datasets. Should find general use
%for any AP-dependent enhancer data and can be modified with little hassle
%for DV enhancers. This uses the negative control data from
%controlset to perform corrections on the data given in dataset.
%
%PARAMETERS
%dataset: This is a string that is identical to the name of the tab in
%dataStatus.xlsx that you wish to analyze.
%
%
%OPTIONS
% 'nc', nc : Specify the nuclear cycle of interest (12,13,or 14).
% 'justMeans' :
% 'noLoading' :
% 'savePath', path :
%
%OUTPUT
%Does not return anything. Does generate graphs.
%
%Author (contact): Armando Reimer (areimer@berkeley.edu), Yang Joon Kim(yjkim90@berkeley.edu)
%Created: 6/3/2017
%Last Updated: 2/4/19
%
%To do:
%        3) Separate out graphs into functions
%        5) Make zeros in cumulative graph actually all zeros. PRIORITY
%        6) Make sure integration periods are consistent with APDiv times
%        7) Make duration graphs subfunction
%        8) Add ability to plot multiple data sets on same graphs. loop?
%%
nc = 2; % the default nuclear cycle shown is nc13
justMeans = 0;
savePath = '';
noLoading = 0;

for i=1:length(varargin)
    
    if strcmpi(varargin{i},'nc')
        nc = varargin{i+1} - 11; %In CompiledParticles, nc12 is indexed as 1, nc13 as 2, and nc14 as 3.
    elseif strcmpi(varargin{i}, 'justMeans')
        justMeans = 1;
    elseif strcmpi(varargin{i}, 'noLoading')
        noLoading = 1;
    elseif strcmpi(varargin{i}, 'savePath')
        savePath = varargin{i+1};
    end
end

[rawDataPath, ProcPath, DropboxFolder, MS2CodePath,...
    PreProcPath, configValues, movieDatabasePath] = DetermineLocalFolders;

warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
close all;

%Set up treatment data for plotting and analysis
if ischar(dataset)
    data = LoadMS2Sets(dataset);
else
    data = dataset; %this allows you to pass the data directly to the script,
    %bypassing the need to call loadms2sets
end

nSets = length(data);
ap = data(1).APbinID;
numAPBins = length(ap);
Prefix = cell(1, nSets);
channel = 1; %no support for 2 channel data at the moment.

accumIntensity(data, nc, justMeans);

%%
%Experiment fraction on
figure('Name', 'Individual Movie Fraction On')
totalEllipses =  zeros(numAPBins, 1);
ellipsesOn = zeros(numAPBins, 1);
rateSum = zeros(1,numAPBins);
%     tSum = zeros(numAPBins, 1);

g = zeros(numAPBins, 1);
g2 = zeros(numAPBins, 1);
n = zeros(numAPBins, 1);
fSet = zeros(numAPBins, nSets);
for dataSet = 1:nSets
    try
        Prefix{dataSet} = data(dataSet).SetName;
    catch
        Prefix{dataSet} = ['data set: ',int2str(dataSet)];
    end
    d = data(dataSet);
    
    if iscell(data(dataSet).EllipsesOnAP)
        data(dataSet).EllipsesOnAP = data(dataSet).EllipsesOnAP{channel};
    end
    if iscell(data(dataSet).TotalEllipsesAP)
        data(dataSet).TotalEllipsesAP = data(dataSet).TotalEllipsesAP{channel};
    end
    
    f = data(dataSet).EllipsesOnAP(:,nc)./data(dataSet).TotalEllipsesAP(:,nc);
    
    totalEllipses = totalEllipses + data(dataSet).TotalEllipsesAP(:,nc);
    ellipsesOn = ellipsesOn + data(dataSet).EllipsesOnAP(:,nc);
    
    rateAr = data(dataSet).rateOnAP{channel}(:,nc);
    b = (rateAr.*data(dataSet).EllipsesOnAP(:,nc))';
    c = vertcat(rateSum,b);
    rateSum = nansum(c);
    %         tSum = tSum + data(dataSet).timeOnOnAP{channel}.*data(dataSet).EllipsesOnAP(:,nc);
    
    %         rateStd = data(dataSet).rateOnAP{channel}.*data(dataSet).EllipsesOnAP
    
    nonanf = f;
    nonanf(isnan(nonanf))=0;
    g = g + nonanf;
    g2 = g2 + nonanf.^2;
    n = ~isnan(f) + n;
    if ~justMeans
        plot(ap,f,'-o','DisplayName',Prefix{dataSet});
        hold on
    end
    fSet(:, dataSet) = f;
end
rateSum = rateSum';
rateMean = rateSum ./ ellipsesOn;
%     timeMean = tSum ./ ellipsesOn;

rateChiSq = zeros(1,numAPBins);

for dataSet = 1:nSets
    %note to AR tomorrow. This step probably is not right syntax. want
    %to do elementwise subtraction
    rateCellAr = (data(dataSet).rateOnAPCell{channel}(:,nc))';
    
    for k = 1:length(rateCellAr)
        if ~isempty(rateCellAr{k})
            deviation = rateCellAr{k} - rateMean(k);
            rateChiSq(k) = rateChiSq(k) + sum(deviation.^2);
        end
    end
end

rateChiSq = rateChiSq';
rateStD = sqrt(rateChiSq ./ ellipsesOn);
rateStE = rateStD ./ sqrt(ellipsesOn);

n(~n) = 1;
fmean = g./n;
fstde = sqrt(g2./n - fmean.^2) ./ sqrt(n);
if nSets > 1
    %         e = errorbar(ap, fmean, fstde,'DisplayName', 'mean $\pm$ std. error');
    %         fMeanNaN = nanmean(fSet, 2);
    %         e = plot(ap, fMeanNaN,'DisplayName', 'mean');
end
legend('show')
hold off

figure('Name','Combined Fraction of On')
%     ellipsesOn(totalEllipses < 3) = NaN;
%     totalEllipses(totalEllipses < 3) = NaN;
%     fstde(totalEllipses < 3) = NaN;
idx = ~any(isnan(totalEllipses),2);
errorbar(ap(idx),ellipsesOn(idx)./totalEllipses(idx), fstde(idx));
xlim([0.275 0.7])
ylim([0, .4])
lgd2 = legend('mean $\pm$ std. error');
set(lgd2, 'Interpreter', 'latex');
title(['fraction of actively transcribing nuclei, nuclear cycle ',num2str(nc+11)]);
xlabel('fraction embryo length');
ylabel('fraction on');
standardizeFigure(gca, legend('show'),'red');
e.Color = [213,108,85]/256;

figure('Name','rate vs ap')
%     ellipsesOn(totalEllipses < 3) = NaN;
%     totalEllipses(totalEllipses < 3) = NaN;
%     fstde(totalEllipses < 3) = NaN;
idx = ~any(isnan(rateMean),2);
apidx = ap(idx);rateMeanidx = rateMean(idx);rateStEidx=rateStE(idx);
e = errorbar(apidx(apidx<=.55),rateMeanidx(apidx<=.55), rateStEidx(apidx<=.55));
xlim([.275, .7])
%     ylim([0, 1000])
lgd2 = legend('mean $\pm$ std. error');
set(lgd2, 'Interpreter', 'latex');
title(['single trace loading rates of ellipses flagged as on, nuclear cycle ',num2str(nc+11)]);
xlabel('fraction embryo length');
ylabel('pol II loading rate (a.u./min)');
standardizeFigure(gca, legend('show'),'red');
e.Color = [213,108,85]/256;
%     e.MarkerSize= 40;
%     e.LineWidth= 4;


%% Combined Ellipses Count
% This is for a quick visual check of the fraction on plots
figure('Name','Combined Ellipses Count')
plot(ap,totalEllipses,'DisplayName','Total')
hold on
plot(ap,ellipsesOn,'DisplayName','On')
title(['on ellipses and total ellipses, nuclear cycle ',num2str(nc+11)]);
xlabel('fraction embryo length');
ylabel('ellipses')
legend('show')
standardizeFigure(gca, legend('show'),'red');
hold off
%%

% %rate plots
% g = zeros(numAPBins, 1);
%     g2 = zeros(numAPBins, 1);
%     n = zeros(numAPBins, 1);
%     fSet = zeros(numAPBins, nSets);
%
% %      f =
%         nonanf = f;
%         nonanf(isnan(nonanf))=0;
%         g = g + nonanf;
%         g2 = g2 + nonanf.^2;
%         n = ~isnan(f) + n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%Experiment number on
figure('Name', 'number_spots')
g = zeros(numAPBins, 1);
g2 = zeros(numAPBins, 1);
n = zeros(numAPBins, 1);
for dataSet = 1:nSets
    f = data(dataSet).EllipsesOnAP(:,nc);
    nonanf = f;
    nonanf(isnan(nonanf))=0;
    g = g + nonanf;
    g2 = g2 + nonanf.^2;
    n = ~isnan(f) + n;
    plot(ap, f, '-o','DisplayName', Prefix{dataSet});
    hold on
end
n(~n) = 1;
fmean = g./n;
fstde = sqrt(g2./n - fmean.^2) ./ sqrt(n);
if nSets > 1
    errorbar(ap, fmean, fstde,'DisplayName', 'mean $\pm$ std. error');
end

hold off
lgd2 = legend('show');
set(lgd2, 'Interpreter', 'Latex');
xlim([.1, .8])
if max(f) ~= 0
    ylim([0, max(f)*1.1]);
else
end
title(['number of actively transcring nuclei, nuclear cycle ',num2str(nc+11)]);
xlabel('fraction embryo length');
ylabel('number on');
standardizeFigure(gca, legend('show'), 'red');

%%
if ~noLoading
    try
        %             pol2LoadingAnalysis(dataset);
    catch
    end
end

analyzeContiguity(d);
plotWindowTimings(d);

%saving every figure
if ~isempty(savePath)
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = [date, '_', dataset, '_', get(FigHandle, 'Name')];
        savefig(FigHandle, [savePath,filesep, APAnalysis, filesep, FigName, '.fig']);
        saveas(FigHandle, [savePath,filesep,APAnalysis, FigName, '.png']);
    end
end

end