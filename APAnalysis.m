function APAnalysis(dataset, varargin)
%APAnalysis(dataset, controlset)
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
%No options
%
%OUTPUT
%Does not return anything. Does generate graphs. 
%
%Author (contact): Armando Reimer (areimer@berkeley.edu)
%Created: 6/3/2017
%Last Updated: 7/5/17
%
%To do: 
%        1) Save graphs somewhere automatically 
%        2) Fix stde error bars
%        3) Separate out graphs into functions
%        4) Put control stuff in another script or subfunctions
%        5) Make zeros in cumulative graph actually all zeros
%        6) Make sure integration periods are consistent with APDiv times
%         7) Make duration graphs subfunction
%% 
    control = '';
    nc = 2;

    for i=1:length(varargin)
        if strcmpi(varargin{i},'control')
            control = varargin{i+1};
        elseif strcmpi(varargin{i},'nc')
            nc = varargin{i+1} - 11; %Because in the CompiledParticles, nc12 is indexed as 1, nc13 as 2, etc.
        end
    end
        
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
    close all force;

    %Set up treatment data for plotting and analysis
    data = LoadMS2Sets(dataset);
    nSets = length(data);
    ap = data(1).APbinID;
    numAPBins = length(ap);
    Prefix = cell(1, nSets);
    to = 0;
    tf = 0;
    
    %% 
    
    %Do some analysis and plotting of treatment data
    cum = zeros(0,numAPBins);
    for dataSet = 1:nSets
        Prefix{dataSet} = data(dataSet).SetName;
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
    figure(1)
    clf('reset')
    for dataSet = 1:nSets
        s = scatter(ap, cum(dataSet, :), 100, 'k', 'filled', 'DisplayName', Prefix{dataSet});
        hold on
    end
    if nSets > 1
        e = errorbar(ap, cummean, cumstde, 'DisplayName', 'Mean $\pm$ std. error');
    end
        
    hold off
    lgd1 = legend('show');
    set(lgd1, 'Interpreter', 'Latex');
%     xlim([0, 1])
%     ylim([0, max([cummean+abs(cumstde), cum(:)']).*1.1 ])
    title({'total average nuclear intensity across';['anterior-posterior axis, nuclear cycle',num2str(nc+11)]});
    xlabel('fraction embryo length');
    ylabel('intensity (a.u.))');
    standardizeFigure(gca, legend('show'), 'scatter', s, 'red');
    
    %% 
    %Experiment fraction on
    figure(2)
    clf('reset')
    g = zeros(numAPBins, 1);
    g2 = zeros(numAPBins, 1);
    n = zeros(numAPBins, 1);
    for dataSet = 1:nSets
        f = data(dataSet).EllipsesOnAP(:,nc)./data(dataSet).TotalEllipsesAP(:,nc);
        nonanf = f;
        nonanf(isnan(nonanf))=0;
        g = g + nonanf;
        g2 = g2 + nonanf.^2; 
        n = ~isnan(f) + n;
        s = scatter(ap, f, 100, 'k', 'filled', 'DisplayName', Prefix{dataSet});
        hold on
    end
    n(~n) = 1;
    fmean = g./n;
    fstde = sqrt(g2./n - fmean.^2) ./ sqrt(n);
    if nSets > 1
        e = errorbar(ap, fmean, fstde,'DisplayName', 'Mean $\pm$ std. error');
    end
        
    hold off
    lgd2 = legend('show');
    set(lgd2, 'Interpreter', 'Latex');
%     xlim([.1, .8])
%     ylim([0, 1.3])
    title(['fraction of actively transcring nuclei, nuclear cycle',num2str(nc+11)]);
    xlabel('fraction embryo length');
    ylabel('fraction on');
    standardizeFigure(gca, legend('show'), 'scatter', s, 'red');
    
    
    %% 
 
    %Experiment number on
    figure(3)
    clf('reset')
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
        s = scatter(ap, f,100, 'k', 'filled', 'DisplayName', Prefix{dataSet});
        hold on
    end
    n(~n) = 1;
    fmean = g./n;
    fstde = sqrt(g2./n - fmean.^2) ./ sqrt(n);
    if nSets > 1
        e = errorbar(ap, fmean, fstde,'DisplayName', 'Mean $\pm$ std. error');
    end
        
    hold off
    lgd2 = legend('show');
    set(lgd2, 'Interpreter', 'Latex');
%     xlim([.1, .8])
%     ylim([0, 1.3])
    title(['numer of actively transcring nuclei, nuclear cycle',num2str(nc+11)]);
    xlabel('fraction embryo length');
    ylabel('number on');
    standardizeFigure(gca, legend('show'), 'scatter', s, 'red');    
  
%%
    analyzeContiguity(d);
    plotWindowTimings(d);

end
%%
function plotWindowTimings(movie)
 
    %movie is the data set we want to look into. 
    
    compiledParticles = movie.CompiledParticles;
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
    
    figure()
    subplot(1, 3, 1)
    h = histogram(onTimes);
    title('on times')
    xlabel('on time (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'bar', h, 'red')
    subplot(1, 3, 2)
    h = histogram(offTimes);
    title('off times')
    xlabel('off time (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'bar', h, 'red')
    subplot(1, 3, 3)
    h = histogram(duration);
    title('duration')
    xlabel('duration of transcription (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'bar', h, 'red')


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
    
    for i = 1:nParticles
       
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

    figure()
    subplot(2, 2, 1)
    h = histogram(contiguity);
    title({'contiguity of traces relative to';' trace length weighted by'; 'length of gaps'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'bar', h, 'red')
    subplot(2, 2, 2)
    h = histogram(contiguityLong);
    title({'contiguity of traces > 1 frame';'relative to trace length';'weighted by lengths of gaps'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'bar', h, 'red')
    subplot(2, 2, 3)
    h = histogram(contiguity2);
    title({'contiguity of traces';'relative to trace length'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'bar', h, 'red')
    subplot(2, 2, 4)
    h = histogram(contiguity2Long);
    title({'contiguity of traces > 1 frame';' relative to trace length'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'bar', h, 'red')
    
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
