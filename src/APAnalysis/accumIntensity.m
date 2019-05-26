function accumIntensity(data, nc, justMeans)
    %Do some analysis and plotting of treatment data
    
     nSets = length(data);
    ap = data(1).APbinID;
    numAPBins = length(ap);
    Prefix = cell(1, nSets);
    channel = 1; %no support for 2 channel data at the moment. 
    
    cum = zeros(0,numAPBins);
    rate = zeros(0, numAPBins);
    integrationFrames = [1, 0];
    timeIntegrationFig = figure('Name', 'time integrated intensity');
    timeIntAx = axes(timeIntegrationFig);
    rateFig = figure('Name', 'rate');
    rateAx = axes(rateFig);
    
    
    dataSets = 1:nSets;
    for dataSet = dataSets
        try
            Prefix{dataSet} = data(dataSet).SetName;
        catch 
            Prefix{dataSet} = ['data set: ',int2str(dataSet)];
        end
        d = data(dataSet);
        if nc==1
            integrationFrames = d.nc12:d.nc13-1;
        elseif nc==2
            integrationFrames = d.nc13:d.nc14-1;
        elseif nc==3
            integrationFrames = d.nc14:length(d.ElapsedTime);
        end
        
        if integrationFrames(1) == 0
            integrationFrames(1) = 1;
        end

        
        fluo = [];
        if ~isempty(d.MeanVector3DAP)
            if iscell(d.MeanVector3DAP)
                d.MeanVector3DAP = d.MeanVector3DAP{channel};
            end
            fluo = d.MeanVector3DAP(integrationFrames,:);
            fluo(isnan(fluo)) = 0;
        end
        
        for APBin = 1:numAPBins
            if ~isempty(fluo)
                cum(dataSet,APBin) = trapz(d.ElapsedTime(integrationFrames),fluo(:,APBin));
                rate(dataSet, APBin) = max(fluo(:,APBin));
            else
                cum(dataSet,APBin) = 0;
            end
        end
    end
    
    cummean = zeros(1,numAPBins);
    cumstd = zeros(1,numAPBins);
    ratemean = zeros(1, numAPBins);
    
    for APBin = 1:numAPBins
        cummean(1, APBin) = nanmean(cum(:,APBin));
        cumstd(1, APBin) = nanstd(cum(:,APBin));
        ratemean(1, APBin) = nanmean(rate(:,APBin));
        
        if cummean(APBin) == 0
            cummean(APBin) = NaN;
        end
         if cumstd(APBin) == 0
            cumstd(APBin) = NaN;
         end
         if ratemean(APBin) == 0
            ratemean(APBin) = NaN;
         end
         cumstde(APBin) = cumstd(APBin) /  sqrt(sum(cum(:,APBin) ~= 0));
    end
    if ~justMeans
        for dataSet = 1:nSets
            plot(timeIntAx, ap,cum(dataSet,:),'-o','DisplayName',Prefix{dataSet});
            hold(timeIntAx, 'on');
            plot(rateAx, ap,rate(dataSet,:),'-o','DisplayName',Prefix{dataSet});
            hold(rateAx, 'on');
        end
    end
    if nSets > 1
        timeIntegrationErrorPlot = errorbar(timeIntAx,ap, cummean, cumstde, 'DisplayName', 'mean $\pm$ std. error');
        ratePlot = plot(rateAx,ap, ratemean, 'DisplayName', 'mean $\pm$ std. error');

    end
        
    hold(timeIntAx, 'off');
    lgd1 = legend(timeIntAx,'show');
    set(lgd1, 'Interpreter', 'Latex');
    xlim(timeIntAx,[.275, .7])
%     ylim(timeIntAx, [0, 100]);
try
    ylim(timeIntAx, [0, max(cummean)]);
end
%     try
%         ylim(timeIntAx,[0, max([cummean+abs(cumstde), cum(:)']).*1.1 ])
%     catch
%         %
%     end
    title(timeIntAx,{'total average nuclear intensity across';['anterior-posterior axis, nuclear cycle ',num2str(nc+11)]});
    xlabel(timeIntAx,'fraction embryo length');
    ylabel(timeIntAx,'intensity (a.u.)');
    standardizeFigure(timeIntAx, legend('show'), 'red');
    timeIntegrationErrorPlot.Color = [213,108,85]/256;
    
    standardizeFigure(rateAx, legend('show'), 'red');
    ratePlot.Color = [213,108,85]/256;