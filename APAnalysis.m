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
%Last Updated: 6/3/17
 
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
    
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
    cumsum = zeros(1,len);
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

    figure(1)
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
    

    figure(2)
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
    
    %This figure will correct figure (2) using the negative control's  fraction on. 
    figure(3)
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
    
   %Plotting the time projection(s) 
   for i = 1:nsets
       figure()
       stringTemp = Prefix{i};
       currentPrefix = stringTemp(11:length(stringTemp)-1);
       [maxTProj, medianTProj] = timeProjection(currentPrefix,'Justnc13');
       imshow(maxTProj,[0 80],'Border','Tight')
       title(Prefix{i})
   end
   
end
