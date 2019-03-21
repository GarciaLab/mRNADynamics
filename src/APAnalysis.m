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
%controlset: This is the same string as the name of the tab
%in datastatus.xlsx that contains your negative control data. Note that
%the datastatus.xlsx doesn't actually have to be the same file as the one 
%used for dataset since LoadMS2Sets searches all similar files for the tab 
%with the right name.
%
%OPTIONS
% 'control', control : 
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
%        4) Put control stuff in another script or subfunctions
%        5) Make zeros in cumulative graph actually all zeros. PRIORITY
%        6) Make sure integration periods are consistent with APDiv times
%        7) Make duration graphs subfunction
%        8) Add ability to plot multiple data sets on same graphs. loop?
%% 
    control = '';
    nc = 2; % the default nuclear cycle shown is nc13  
    justMeans = 0;
    savePath = '';
    noLoading = 0;

    for i=1:length(varargin)
        if strcmpi(varargin{i},'control')
            control = varargin{i+1};
        elseif strcmpi(varargin{i},'nc')
            nc = varargin{i+1} - 11; %In CompiledParticles, nc12 is indexed as 1, nc13 as 2, and nc14 as 3.
        elseif strcmpi(varargin{i}, 'justMeans')
            justMeans = 1;
        elseif strcmpi(varargin{i}, 'noLoading')
            noLoading = 1;
        elseif strcmpi(varargin{i}, 'savePath')
            savePath = varargin{i+1};
        end
    end
        
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
    to = 0;
    tf = 0;
    channel = 1; %no support for 2 channel data at the moment. 
    
    %% 
    
    %Do some analysis and plotting of treatment data
    cum = zeros(0,numAPBins);
    for dataSet = 1:nSets
        try
            Prefix{dataSet} = data(dataSet).SetName;
        catch 
            Prefix{dataSet} = ['data set: ',int2str(dataSet)];
        end
        d = data(dataSet);
        lastFrame = length(d.ElapsedTime);
        if nc==1
            to = d.nc12;
            tf = d.nc13-1;
        elseif nc==2
            to = d.nc13;
            tf = d.nc14-1;
        elseif nc==3
            to = d.nc14;
            tf = lastFrame;
        end
        if to == 0
            to = 1; %to fix an indexing issue 
        end
        if iscell(d.MeanVectorAP)
            d.MeanVectorAP = d.MeanVectorAP{channel};
        end
        fluo = d.MeanVectorAP(to:tf ,:);
        fluo2 = fluo;
        fluo(isnan(fluo)) = 0;
        for APBin = 1:numAPBins
            cum(dataSet,APBin) = trapz(d.ElapsedTime(to:tf),fluo(:,APBin));
        end
    end
    cummean = zeros(1,numAPBins);
    cumstd = zeros(1,numAPBins);
    for APBin = 1:numAPBins
        cummean(1, APBin) = nanmean(cum(:,APBin));
        cumstd(1, APBin) = nanstd(cum(:,APBin));
        if ~cummean(APBin)
            cummean(APBin) = NaN;
        end
         if ~cumstd(APBin)
            cumstd(APBin) = NaN;
         end
         cumstde(APBin) = cumstd(APBin) /  sqrt(sum(cum(:,APBin) ~= 0));
    end
    figure('Name', 'intensity')
    if ~justMeans
        for dataSet = 1:nSets
            plot(ap,cum(dataSet,:),'-o','DisplayName',Prefix{dataSet});
            hold on
        end
    end
    if nSets > 1
        e = errorbar(ap, cummean, cumstde, 'DisplayName', 'mean $\pm$ std. error');
    end
        
    hold off
    lgd1 = legend('show');
    set(lgd1, 'Interpreter', 'Latex');
    xlim([0, 1])
    ylim([0, max([cummean+abs(cumstde), cum(:)']).*1.1 ])
    title({'total average nuclear intensity across';['anterior-posterior axis, nuclear cycle ',num2str(nc+11)]});
    xlabel('fraction embryo length');
    ylabel('intensity (a.u.)');
    standardizeFigure(gca, legend('show'), 'red');
    e.Color = [213,108,85]/255;
    
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
    xlim([0.15 0.6])%xlim([.275, .65])
    ylim([0, 1.1])
    lgd2 = legend('mean $\pm$ std. error');
    set(lgd2, 'Interpreter', 'latex');
    title(['fraction of actively transcribing nuclei, nuclear cycle ',num2str(nc+11)]);
    xlabel('fraction embryo length');
    ylabel('fraction on');
    standardizeFigure(gca, legend('show'),'red');
    e.Color = [213,108,85]/255;
    
 figure('Name','rate vs ap')
%     ellipsesOn(totalEllipses < 3) = NaN;
%     totalEllipses(totalEllipses < 3) = NaN;
%     fstde(totalEllipses < 3) = NaN;
    idx = ~any(isnan(rateMean),2);
    apidx = ap(idx);rateMeanidx = rateMean(idx);rateStEidx=rateStE(idx);
    e = errorbar(apidx(apidx<=.55),rateMeanidx(apidx<=.55), rateStEidx(apidx<=.55));
    xlim([.25, .7])
    ylim([0, 1000])
    lgd2 = legend('mean $\pm$ std. error');
    set(lgd2, 'Interpreter', 'latex');
    title(['single trace loading rates of ellipses flagged as on, nuclear cycle ',num2str(nc+11)]);
    xlabel('fraction embryo length');
    ylabel('pol II loading rate (a.u./min)');
    standardizeFigure(gca, legend('show'),'red');
    e.Color = [213,108,85]/255;
    e.MarkerSize= 40;
    e.LineWidth= 4;

        
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
            pol2LoadingAnalysis(dataset);
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
              savefig(FigHandle, [savePath,filesep, FigName, '.fig']);
              saveas(FigHandle, [savePath,filesep, FigName, '.png']);
        end
    end

end
%%
function plotWindowTimings(movie)
 
    channel = 1; %no support for 2 channel
    %movie is the data set we want to look into. 
    compiledParticles = movie.CompiledParticles;
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{channel};
    end
    nParticles = length(compiledParticles);
    onTimes = [];
    offTimes = [];
    durations = [];
    
    for i = 1:nParticles
        frames = compiledParticles(i).Frame; 
        onTimes(i) = frames(1);
        offTimes(i) = frames(end);
        duration(i) = frames(end) - frames(1);
    end
    
    figure('Name', 'timings')
    subplot(1, 3, 1)
    h = histogram(onTimes);
    title('on times')
    xlabel('on time (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(1, 3, 2)
    h = histogram(offTimes);
    title('off times')
    xlabel('off time (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(1, 3, 3)
    h = histogram(duration);
    title('duration')
    xlabel('duration of transcription (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')


end
    %%
function analyzeContiguity(movie)

    %movie is the data set we want to look into. 
    
    compiledParticles = movie.CompiledParticles;
    nParticles = length(compiledParticles);
    contiguityLong = [];
    contiguity2Long = [];
    contiguity = zeros(1, nParticles);
    contiguity2 = zeros(1, nParticles);
    channel = 1; %no support for 2 channels at the moment
    
    for i = 1:nParticles
        if iscell(compiledParticles)
            compiledParticles = compiledParticles{channel};
        end
        frames = compiledParticles(i).Frame; len = length(frames); contiguity(i) = len; contiguity2(i) = len; missed = 0; missed2 = 0;
       
        if len > 1
           for j = 2:len
               frameInterval = frames(j) - frames(j-1);
               missed = missed + frameInterval - (len + 1);
               if frameInterval > 1 
                    missed2 = missed2 + 1; %counting gaps. 
               end
           end
       end
       
       contiguity(i) = (contiguity(i) - missed)/len; %every gap weighted by its duration
       contiguity2(i) = (contiguity2(i) - missed2)/len; %this counts the gaps and normalizes by full trace length

       if len > 1
           contiguityLong = [contiguityLong,contiguity(i)]; %every gap weighted by its duration but traces one frame long excluded
           contiguity2Long = [contiguity2Long,contiguity2(i)]; %analagous to above
       end
    end

    figure('Name', 'contiguity')
    subplot(2, 2, 1)
    h = histogram(contiguity);
    title({'contiguity of traces relative to';' trace length weighted by'; 'length of gaps'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(2, 2, 2)
    h = histogram(contiguityLong);
    title({'contiguity of traces > 1 frame';'relative to trace length';'weighted by lengths of gaps'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(2, 2, 3)
    h = histogram(contiguity2);
    title({'contiguity of traces';'relative to trace length'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(2, 2, 4)
    h = histogram(contiguity2Long);
    title({'contiguity of traces > 1 frame';' relative to trace length'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    
end

%     %% 
%     %Set up control data 
%     if ~isempty(control)
%         data0 = LoadMS2Sets(controlset);
%         if numAPBins ~= length(data0(1).APbinID)
%             error('Experiment and negative control have different number of AP bins. Can''t make graphs.')
%         end   
%         nsets0 = length(data0);
%         Prefix0 = cell(1, nsets0);
%         for dataSet = 1:nsets0
%             Prefix0{dataSet} = data0(dataSet).SetName;
%         end
%     end
%     
%     %% 
%     if ~isempty(control)
%         %Control fraction on analysis
%         figure()
%         clf('reset')
%         g = zeros(numAPBins, 1);
%         g2 = zeros(numAPBins, 1);
%         n = zeros(numAPBins, 1);
%         for i = 1:nsets0
%             f = data0(i).EllipsesOnAP(:,nc)./data0(i).TotalEllipsesAP(:,nc); %nc13
%             nonanf = f;
%             nonanf(isnan(nonanf))=0;
%             g = g + nonanf;
%             g2 = g2 + nonanf.^2; 
%             n = ~isnan(f) + n;
%             scatter(ap, f,'DisplayName', Prefix0{i})
%             hold on
%         end
%         n(~n) = 1;
%         fmean = g./n;
%         fstde = sqrt(g2./n - fmean.^2) ./ sqrt(n);
%         errorbar(ap, fmean, fstde,'DisplayName', 'Mean $\pm$ std. error')
%         hold off
%         lgd2 = legend('show');
%         set(lgd2, 'Interpreter', 'Latex')
%         xlim([.1, .8])
%         ylim([0, 1.3])
%         title('Fraction On Control')
%         xlabel('fraction embryo length')
%         ylabel('fraction on')
%         standardizeFigure(gca, legend('show'))
% 
%         %This figure will correct fraction on from the experiment using the negative control's  fraction on. 
%         figure()
%         clf('reset')
%         g = zeros(numAPBins, 1);
%         g2 = zeros(numAPBins, 1);
%         n = zeros(numAPBins, 1);
%         for i = 1:nsets0
%             f = data0(i).EllipsesOnAP(:,nc)./data0(i).TotalEllipsesAP(:,nc);
%             nonanf = f;
%             nonanf(isnan(nonanf))=0;
%             g = g + nonanf;
%             g2 = g2 + nonanf.^2; 
%             n = ~isnan(f) + n;
%     %         scatter(ap, f, 'LineWidth', 1, 'DisplayName', Prefix0{i}) %these
%     %         really need to be individually corrected by the control mean that's
%     %         calculated below. i just didn't want to have to make a separate l
%     %         oop. 
%             hold on
%         end
%         n(~n) = 1;
%         fmean0 = g./n;    
%         fstde0 = sqrt(g2./n - fmean0.^2) ./ sqrt(n);
% 
%         fcorr = fmean - fmean0; %first order approximation in the limit of no enhancer coupling
%     %%    
%         %Here's the full correction 
% 
%         % KD(f - fo) + f*alpha - r*alpha = 0
% 
%         correctionEqn = @(params) params(1)*(fmean - fmean0) + params(2)*fmean -...
%             params(2)*params(3); %params(1) = KD. params(2) = alpha. params(3) = r
%         initParams = [1, 1, 1]; %will figure out units later   
%         lsqOptions=optimset('Display','none','maxfunevals',10000,'maxiter',10000);
%         [fit, res1, residual, exitflag, output, ~, jacobian] = lsqnonlin(correctionEqn, ...
%             initParams,zeros(1,3),inf(1,3), lsqOptions);
%         KD = fit(1);
%         r = fit(3);
%         confidence_intervals = nlparci(fit,residual,'jacobian',jacobian);
%         errors = zeros(1, length(fit));
%         for i = 1:length(confidence_intervals)
%             errors(i) = abs((abs(confidence_intervals(i, 1)) - abs(confidence_intervals(i, 2)))/2);
%         end
%         relative_errors = abs(errors./fit);  
%         %%
% 
%         totalerr = sqrt(fstde0.^2 + fstde.^2);
% 
%         errorbar(ap, fcorr, totalerr, 'LineWidth', 5, 'DisplayName', 'Mean $\pm$ std. error')
%         hold off
%         lgd3 = legend('show');
%         set(lgd3, 'Interpreter', 'Latex')
%         xlim([.1, .8])
%         ylim([-1.5, 1.5])
%         title('Fraction On Corrected')
%         xlabel('Fraction EL')
%         ylabel('Fraction on corrected')
%         standardizeFigure(gca, lgd3)
% 
%         %%
%         %Enhancer-specific modelling
%         BcdMax = 120; %nM
%         lambda = .2;
%         Bcd = BcdMax*exp(-ap/lambda); 
%         f1A3 = @(constants) constants(1)*Bcd./(Bcd+constants(2)); %Unitless fraction active nuclei
%         rateCoeff = 1; %mRNA / unit time
%         rate1A3 =  @(constants) rateCoeff*constants(1)*Bcd./(Bcd+constants(2)); %dmRNA/dt as a function of activator concentration
% 
%         if contains(dataset, 'v7')
%             %linear fit to Langmuir's isotherm to get thermodynamic parameters
%             %directly. 1/f = 1/r + (KD/r)Bcd 
%             f1 = (1./fmean(8:end))'; %8th bin is roughly 20% AP where terminal system drops off
%             k = ~isinf(f1);
%             Bcd1 = 1./Bcd(8:end);
%             [p, S] = polyfit(Bcd1(k),f1(k),1);        
%             r = 1/p(1);
%             KD = p(2)*r;
%             fFitLinear = r*(Bcd./(Bcd+KD));
%             figure()
%             scatter(ap,fFitLinear)
%             title('v7 Fit')
%             xlabel('Fraction EL')    
%             hold on
%             scatter(ap, fmean)
%             legend('fit', 'data')
% 
%             %not sure why the linear fit results were bad. gonna try nonlinear
%             fv7 = @(constants) constants(1)*Bcd(8:end)./(Bcd(8:end)+constants(2)) - fmean(8:end); %Unitless fraction active nuclei
%             initParams = [.5, .3];  
%             lsqOptions=optimset('Display','none','maxfunevals',10000,'maxiter',10000);
%             [fit, res1, residual, exitflag, output, ~, jacobian] = lsqnonlin(fv7, ...
%                 initParams,[0,.9],[1,120]);
%             r_nonlin = fit(1);
%             disp(r_nonlin);
%             KD_nonlin = fit(2);
%             disp(KD_nonlin);
%             fFitNonLinearBicoid = r_nonlin*(Bcd./(Bcd+KD_nonlin));    
%             fFitNonLinearAP = r_nonlin*(BcdMax*exp(-ap/lambda)./(BcdMax*exp(-ap/lambda)+KD_nonlin));     
%             figure()
%             scatter(Bcd,fFitNonLinearBicoid)
%             title('v7 nonlinear fit')
%             xlabel('Bcd (nM)')  
%             ylabel('Fraction Active')
%             figure()
%             scatter(ap,fFitNonLinearAP, 'LineWidth', 5)
%             title('v7 nonlinear fit')
%             xlabel('Fraction EL')  
%             ylabel('Fraction Active')
%             hold on
%             scatter(ap, fmean, 'LineWidth', 5)
%             legend('fit', 'data')
%             standardizeFigure(gca, legend('show'))
% 
%         else
%             figure()
%             scatter(ap,f1A3([KD, r]))
%             title('Fitted 1A3 Fraction On')
%             xlabel('Fraction EL')
%             figure()
%             scatter(ap, rate1A3([KD, r]))
%             title('Fitted 1A3 dmRNA/dt')
%             xlabel('Fraction EL')
%         end
%     end
%         %%
% %     %plotting the time projection(s) 
% %     for i = 1:nsets
% %         figure()
% %         stringTemp = Prefix{i};
% %         currentPrefix = stringTemp(11:length(stringTemp)-1);
% %         [maxTProj, medianTProj] = timeProjection(currentPrefix,'Justnc13');
% %         imshow(maxTProj,[0 80],'Border','Tight')
% %         title(Prefix{i})
% %     end
% end
%% Averaging multiple datasets (This is written by Yang Joon Kim)
function AverageDatasets(DataType,varargin)

% DESCRIPTION
% This function has input of datatype in DataStatus.xls, grabs all datasets
% in that tab,and calculates 
% 1) Weighted Sum of averaged MS2 spot fluorescence,Standard
% Deviation, and the total number of MS2 spots from multiple embryos in
% nc13 and nc14.
% 2) Accumulated mRNA (Accumulated fluorescence over nc13 and nc14)
% This is edited from Meghan's CombineMultipleEmbryos.m script

% Note. This is assuming that you're interested in nc13 and nc14, thus you
% can just edit this to include nc12, but you should make it work for the
% datasets that doesn't have the whole nc12, probably making it as an
% option.

% To do : if the start of nc12 was caught, it can combine nc12 data as well.
%
% OPTIONS
% No Options (This can be edited later, like for nc12 or sth else)

% PARAMETERS
% DataType: This is a string that is identical to the name of the tab in
% dataStatus.xlsx that you wish to analyze.
%
% OUTPUT
% Variables for plotting, or more detailed analysis with the Averaged spot
% fluorescence over time. Save as 'Name of the DataType'.mat file
% (nc12, nc13, nc14, NParticlesAP,MeanVectorsAP, SDVectorAP, ElapsedTime) 
% corresponding to the embryos combined

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Data = LoadMS2Sets(DataType);

% Save path option
savePath = '';

    for i=1:length(varargin)
        if strcmpi(varargin{i}, 'savePath')
            savePath = varargin{i+1};
        end
    end

numEmbryos=length(Data);

%Find the total number of frames for each embryo
numFrames = zeros(1, numEmbryos);
maxAPIndex = zeros(1, numEmbryos);
maxTime = zeros(1, numEmbryos);
for i = 1:numEmbryos
    numFrames(i) = size(Data(i).ElapsedTime, 2);
    nc13(i) = Data(i).nc13;
    lengthNC1314(i) = numFrames(i) - nc13(i)+1; % length of frames from nc13 to the last frame
    maxAPIndex(i) = Data(i).MaxAPIndex;
    maxTime(i) = Data(i).ElapsedTime(numFrames(i));
end

%Store the number of AP bins (this should always be 41).
numAPBins = maxAPIndex(1);

%Store all variables to be combined in a single structure. 
combinedData = struct('ElapsedTime',{},'NParticlesAP',{},...
                        'MeanVectorAP',{},'SDVectorAP',{});

%% Synchronize the vectors as the beginning of the nc 13
% This should be edited to include the nc12 or even earlier in the future.
% For now, nc13 and nc14 might be good enough.

% Define the new ElapsedTime vector for the combined embryo. 
% The new ElapsedTime should start with the beginning of nc13, and also has
% the length of the frames of nc13 + nc14 (of the longest dataset), all the
% empty values can be plugged with Nans.

% This ElapsedTime variable has evenly spaced time points estimated from diff(ElapsedTime)
% (This is assumption that we took the data with negligible time between serieses, which is pretty fair)
 
% Calculate the length of the new ElapsedTime, and also from which dataset
[NewFrameLength, Index] = max(lengthNC1314);

% Define an empty matrices (filled with Nans)
MeanVectorAP = NaN(NewFrameLength,numAPBins,numEmbryos);
SDVectorAP = NaN(NewFrameLength,numAPBins,numEmbryos);
NParticlesAP = NaN(NewFrameLength,numAPBins,numEmbryos);

% Synchornize all fields as all of them starts from nc 13
for i=1:numEmbryos
    if nc13(i)==0
        error('Check the Movie if it really does not start from nc13, then you should edit this code or make that dataset as an exception')
    else
        MeanVectorAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).MeanVectorAP(nc13(i):numFrames(i),:);
        SDVectorAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).SDVectorAP(nc13(i):numFrames(i),:);
        NParticlesAP(1:numFrames(i)-nc13(i)+1,:,i) = Data(i).NParticlesAP(nc13(i):numFrames(i),:);
    end
end
        
% % Take the most frequent value of dT from the ElapsedTime. It's because the
% % dT can be different in case we stopped and restarted the movie.
deltaT = mode(diff(Data(1).ElapsedTime)); 
ElapsedTime = deltaT*(0:NewFrameLength-1);


%% Average all fields at each time point

% Make Nans as zeros
MeanVectorAP(isnan(MeanVectorAP)) = 0;
SDVectorAP(isnan(SDVectorAP)) = 0;
NParticlesAP(isnan(NParticlesAP)) = 0;

sumMean = zeros(NewFrameLength,numAPBins);
sumSD = zeros(NewFrameLength,numAPBins);
sumNParticles = zeros(NewFrameLength,numAPBins);

for i=1:numEmbryos
    sumMean = sumMean + squeeze(MeanVectorAP(:,:,i).*NParticlesAP(:,:,i));
    sumSD = sumSD + squeeze(SDVectorAP(:,:,i).^2.*NParticlesAP(:,:,i));
    sumNParticles = sumNParticles + squeeze(NParticlesAP(:,:,i));
end
    
MeanVectorAPTemp = sumMean./sumNParticles;
SDVectorAPTemp = sqrt(sumSD./sumNParticles);
NParticlesAPTemp = sumNParticles;

MeanVectorAP = MeanVectorAPTemp;
SDVectorAP = SDVectorAPTemp;
SEVectorAP = SDVectorAP/sqrt(numEmbryos); % Standard error of mean (SD / sqrt(number of observation)
NParticlesAP = NParticlesAPTemp;
ElapsedTime = ElapsedTime;

%% Accumulate mRNA over time (This can be made as an optional)
% I will calculate the Integrated mRNA from the MeanVectorAP
NFrames = length(ElapsedTime);
nAPbins = max(maxAPIndex);

AccumulatedmRNA = zeros(NFrames,nAPbins);
AccumulatedmRNA_SD =  zeros(NFrames,nAPbins);
MeanVectorAP(isnan(MeanVectorAP))=0;
SDVectorAP(isnan(SDVectorAP))=0;

for i=1:maxAPIndex
    for j=2:length(ElapsedTime)
        AccumulatedmRNA(j,i) = trapz(ElapsedTime(1:j),MeanVectorAP(1:j,i));
        AccumulatedmRNA_SD(j,i) = sqrt(trapz(ElapsedTime(1:j),SDVectorAP(1:j,i).^2));
    end
end

%% Save the fields in .mat file
    if ~isempty(savePath)
        save([savePath,filesep,DataType,'.mat'],...
            'MeanVectorAP','SDVectorAP','SEVectorAP','NParticlesAP','ElapsedTime',...
            'AccumulatedmRNA','AccumulatedmRNA_SD')
    else
        warning('Define the File Path in the argument above')
    end
end



