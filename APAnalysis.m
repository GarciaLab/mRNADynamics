function APAnalysis(dataset, controlset)
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
%To do:  1) Allow control set to be optional. 
%        2) Save graphs somewhere automatically 
%        3) Fix stde error bars
 
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
    close all force
    
    data0 = LoadMS2Sets(controlset);
    data = LoadMS2Sets(dataset);
    ap = data(1).APbinID;
    len = length(ap);
    if len ~= length(data0(1).APbinID)
        error('Experiment and negative control have different number of AP bins. Can''t make graphs.')
    end
    nsets = length(data);
    nsets0 = length(data0);
    Prefix = cell(1, nsets);
    Prefix0 = cell(1, nsets0);
    
    cum = zeros(0,len);
    for j = 1:nsets
        Prefix{j} = data(j).SetName;
        fluo = data(j).MeanVectorAP;
        fluo2 = fluo;
        fluo(isnan(fluo)) = 0;
        for i = 1:size(fluo,2)
            cum(j,i) = trapz(data(j).ElapsedTime,fluo(:,i));
        end
    end
    
    for j = 1:nsets0
        Prefix0{j} = data0(j).SetName;
    end

    cummean = zeros(1,len);
    cumstd = zeros(1,len);
    for i = 1:length(cum)
        cummean(1, i) = nanmean(cum(:,i));
        cumstd(1, i) = nanstd(cum(:,i));
        if ~cummean(i)
            cummean(i) = NaN;
        end
         if ~cumstd(i)
            cumstd(i) = NaN;
         end
         cumstde(i) = cumstd(i) /  sqrt(sum(cum(:,i) ~= 0));
    end

    figure()
    clf('reset')
    for i = 1:nsets
        plot(ap, cum(i, :), 'LineWidth', 1, 'DisplayName', Prefix{i})
        hold on
    end
    errorbar(ap, cummean, cumstde, 'LineWidth', 5, 'DisplayName', 'Mean $\pm$ std. error')

    hold off
    lgd1 = legend('show');
    set(lgd1, 'Interpreter', 'Latex')
    xlim([.1, .8])
    ylim([0, max([cummean+abs(cumstde), cum(:)']).*1.1 ])
    title('Total average nuclear intensity across AP')
    xlabel('Fraction EL')
    ylabel('Intensity (A.U.)')
    standardizeFigure(gca, lgd1)
    
    
    %Control fraction on
    figure()
    clf('reset')
    g = zeros(len, 1);
    g2 = zeros(len, 1);
    n = zeros(len, 1);
    for i = 1:nsets0
        f = data0(i).EllipsesOnAP(:,2)./data0(i).TotalEllipsesAP(:,2);
        nonanf = f;
        nonanf(isnan(nonanf))=0;
        g = g + nonanf;
        g2 = g2 + nonanf.^2; 
        n = ~isnan(f) + n;
        plot(ap, f, 'LineWidth', 1, 'DisplayName', Prefix0{i})
        hold on
    end
    n(~n) = 1;
    fmean = g./n;
    fstde = sqrt(g2./n - fmean.^2) ./ sqrt(n);
    errorbar(ap, fmean, fstde, 'LineWidth', 5, 'DisplayName', 'Mean $\pm$ std. error')
    hold off
    lgd2 = legend('show');
    set(lgd2, 'Interpreter', 'Latex')
    xlim([.1, .8])
    ylim([0, 1.3])
    title('Fraction On Control')
    xlabel('Fraction EL')
    ylabel('Fraction on')
    standardizeFigure(gca, lgd2)
    
    %Experiment fraction on
    figure()
    clf('reset')
    g = zeros(len, 1);
    g2 = zeros(len, 1);
    n = zeros(len, 1);
    for i = 1:nsets
        f = data(i).EllipsesOnAP(:,2)./data(i).TotalEllipsesAP(:,2);
        nonanf = f;
        nonanf(isnan(nonanf))=0;
        g = g + nonanf;
        g2 = g2 + nonanf.^2; 
        n = ~isnan(f) + n;
        plot(ap, f, 'LineWidth', 1, 'DisplayName', Prefix{i})
        hold on
    end
    n(~n) = 1;
    fmean = g./n;
    fstde = sqrt(g2./n - fmean.^2) ./ sqrt(n);
    errorbar(ap, fmean, fstde, 'LineWidth', 5, 'DisplayName', 'Mean $\pm$ std. error')
    hold off
    lgd2 = legend('show');
    set(lgd2, 'Interpreter', 'Latex')
    xlim([.1, .8])
    ylim([0, 1.3])
    title('Fraction On')
    xlabel('Fraction EL')
    ylabel('Fraction on')
    standardizeFigure(gca, lgd2)
    
    
    %This figure will correct fraction on from the experiment using the negative control's  fraction on. 
    figure()
    clf('reset')
    g = zeros(len, 1);
    g2 = zeros(len, 1);
    n = zeros(len, 1);
    for i = 1:nsets0
        f = data0(i).EllipsesOnAP(:,2)./data0(i).TotalEllipsesAP(:,2);
        nonanf = f;
        nonanf(isnan(nonanf))=0;
        g = g + nonanf;
        g2 = g2 + nonanf.^2; 
        n = ~isnan(f) + n;
%         plot(ap, f, 'LineWidth', 1, 'DisplayName', Prefix0{i}) %these
%         really need to be individually corrected by the control mean that's
%         calculated below. i just didn't want to have to make a separate l
%         oop. 
        hold on
    end
    n(~n) = 1;
    fmean0 = g./n;    
    fstde0 = sqrt(g2./n - fmean0.^2) ./ sqrt(n);
    
    fcorr = fmean - fmean0; %first order approximation in the limit of no enhancer coupling
%%    
    %Here's the full correction 
    
    % KD(f - fo) + f*alpha - r*alpha = 0
    
    correctionEqn = @(params) params(1)*(fmean - fmean0) + params(2)*fmean -...
        params(2)*params(3); %params(1) = KD. params(2) = alpha. params(3) = r
    initParams = [1, 1, 1]; %will figure out units later   
    lsqOptions=optimset('Display','none','maxfunevals',10000,'maxiter',10000);
    [fit, res1, residual, exitflag, output, ~, jacobian] = lsqnonlin(correctionEqn, ...
        initParams,zeros(1,3),inf(1,3), lsqOptions);
    KD = fit(1);
    r = fit(3);
    confidence_intervals = nlparci(fit,residual,'jacobian',jacobian);
    errors = zeros(1, length(fit));
    for i = 1:length(confidence_intervals)
        errors(i) = abs((abs(confidence_intervals(i, 1)) - abs(confidence_intervals(i, 2)))/2);
    end
    relative_errors = abs(errors./fit);  
    %%
    
    totalerr = sqrt(fstde0.^2 + fstde.^2);
    
    errorbar(ap, fcorr, totalerr, 'LineWidth', 5, 'DisplayName', 'Mean $\pm$ std. error')
    hold off
    lgd3 = legend('show');
    set(lgd3, 'Interpreter', 'Latex')
    xlim([.1, .8])
    ylim([-1.5, 1.5])
    title('Fraction On Corrected')
    xlabel('Fraction EL')
    ylabel('Fraction on corrected')
    standardizeFigure(gca, lgd3)
    
    %%
    %Enhancer-specific modelling
    BcdMax = 120; %nM
    lambda = .2;
    Bcd = BcdMax*exp(-ap/lambda); 
    f1A3 = @(constants) constants(1)*Bcd./(Bcd+constants(2)); %Unitless fraction active nuclei
    rateCoeff = 1; %mRNA / unit time
    rate1A3 =  @(constants) rateCoeff*constants(1)*Bcd./(Bcd+constants(2)); %dmRNA/dt as a function of activator concentration
    
    if contains(dataset, 'v7')
        %linear fit to Langmuir's isotherm to get thermodynamic parameters
        %directly. 1/f = 1/r + (KD/r)Bcd 
        f1 = (1./fmean(8:end))'; %8th bin is roughly 20% AP where terminal system drops off
        k = ~isinf(f1);
        Bcd1 = 1./Bcd(8:end);
        [p, S] = polyfit(Bcd1(k),f1(k),1);        
        r = 1/p(1);
        KD = p(2)*r;
        fFitLinear = r*(Bcd./(Bcd+KD));
        figure()
        plot(ap,fFitLinear)
        title('v7 Fit')
        xlabel('Fraction EL')    
        hold on
        plot(ap, fmean)
        legend('fit', 'data')

        %not sure why the linear fit results were bad. gonna try nonlinear
        fv7 = @(constants) constants(1)*Bcd(8:end)./(Bcd(8:end)+constants(2)) - fmean(8:end); %Unitless fraction active nuclei
        initParams = [.5, .3];  
        lsqOptions=optimset('Display','none','maxfunevals',10000,'maxiter',10000);
        [fit, res1, residual, exitflag, output, ~, jacobian] = lsqnonlin(fv7, ...
            initParams,[0,.9],[1,120]);
        r_nonlin = fit(1)
        KD_nonlin = fit(2)
        fFitNonLinearBicoid = r_nonlin*(Bcd./(Bcd+KD_nonlin));    
        fFitNonLinearAP = r_nonlin*(BcdMax*exp(-ap/lambda)./(BcdMax*exp(-ap/lambda)+KD_nonlin));     
        figure()
        plot(Bcd,fFitNonLinearBicoid)
        title('v7 nonlinear fit')
        xlabel('Bcd (nM)')  
        ylabel('Fraction Active')
        figure()
        plot(ap,fFitNonLinearAP, 'LineWidth', 5)
        title('v7 nonlinear fit')
        xlabel('Fraction EL')  
        ylabel('Fraction Active')
        hold on
        plot(ap, fmean, 'LineWidth', 5)
        legend('fit', 'data')
        standardizeFigure(gca, legend('show'))

    else
        figure()
        plot(ap,f1A3([KD, r]))
        title('Fitted 1A3 Fraction On')
        xlabel('Fraction EL')
        figure()
        plot(ap, rate1A3([KD, r]))
        title('Fitted 1A3 dmRNA/dt')
        xlabel('Fraction EL')
    end
    %%
   %Plotting the time projection(s) 
%    for i = 1:nsets
%        figure()
%        stringTemp = Prefix{i};
%        currentPrefix = stringTemp(11:length(stringTemp)-1);
%        [maxTProj, medianTProj] = timeProjection(currentPrefix,'Justnc13');
%        imshow(maxTProj,[0 80],'Border','Tight')
%        title(Prefix{i})
%    end
   
   
   
end
