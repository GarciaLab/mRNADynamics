function simulateSingleTraces(params,varargin)
%
% simulateSingleTraces(varargin)
%
% DESCRIPTION
% Generates noisy single traces of a hunchback-like promoter and averages them. 
%
%
% ARGUMENTS
%
% 'displaySingle': display the single traces
% 'displayAverage': display the average traces
%
% OUTPUT
% none
%
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 7/21/2018
% Last Updated: 7/21/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)


displaySingle = 0;
displayAverage = 0;

    for i=1:length(varargin)
        if strcmpi(varargin{i},'displaySingle')
            displaySingle = 1;
        elseif strcmpi(varargin{i},'displayAverage')
            displayAverage = 1;
        end
    end
       
   
    %perform the simulation
    
    %iterate over loading rates
%     paramsIter = params;
% %     rate = .1:.05:20;
%     ap = 0:.01:10;
%     rate = 10*sigmf(ap, [3 6]);
% %     plot(ap*1000, rate)
%     for i = 1:length(rate)
%         paramsIter.MEAN_POL_LOADING_RATE = rate(i);
%         paramsIter.MEAN_POL_OFFLOADING_RATE = rate(i);
%         [fluo(i), trunc(i)] = averageSimulatedTraces(paramsIter, 1);
%     end
%     scatter(ap, fluo)
%     hold on
%     scatter(ap, trunc)
% %     set(gca, 'YScale', 'log')
%     legend('detection-unlimited', 'detection-limited')
%     ylabel('cumulative fluorescence (au)')
%     xlabel('loading rate (pols/min)')
%     standardizeFigure(gca,[])
    
      [fluo, trunc] = averageSimulatedTraces(params, displaySingle, displayAverage);
end

function [timeElapsed, fluos, singleFluosTrunc] = simulateTrace(params, displaySingle)
    
          
    %calibrations and conversions
    timePerFrame = params.TIME_PER_FRAME/60; %min
    t_elon = params.GENE_LENGTH/params.ELON_RATE; %min
    numFrames = params.INTERPHASE / timePerFrame;    
    f_init = params.T_INIT/timePerFrame; f_elon = t_elon/timePerFrame; f_fin = params.T_FIN/timePerFrame;
    frames = 1:numFrames;
    num_steps = 1000;
    time_steps = 1:num_steps;
    timeElapsed = (frames-1)*timePerFrame;
   
    %initialization
    fluos = [];
    fluo = 0;
    pols = [];
    pol = 0;
%     occ = {};
%     occupation_state = zeros(length(frames), 24, params.N_POLS);
%     for f = time_steps
    for f = frames
            
        loadingNoise = params.NOISE_SCALE_BIO*params.MEAN_POL_LOADING_RATE*randn;
        offLoadingNoise = params.NOISE_SCALE_BIO*params.MEAN_POL_OFFLOADING_RATE*randn;
        loadingRate = params.MEAN_POL_LOADING_RATE + loadingNoise;
        offLoadingRate = params.MEAN_POL_OFFLOADING_RATE + offLoadingNoise;          
        scopeNoise = fluo*params.NOISE_SCALE_EXP*randn();
        backGround = params.MEAN_BACKGROUND*randn();
        
        if f < f_init
            %do nothing
        elseif f_init <= f && f <= f_init+f_elon %period of constant rise
            pol = pol + loadingRate;
              %load a new pol; 
              
%               for i = 1:length(pols)
%                   %update each pol
%                   if pols{i} < 24
%                     pols{i} = [pols{i}, 1];
%                   end
%                   for k = 1:length(pols{i})
%                       %update each loop 
%                       if ~pols{i}(k)
%                           mcp_binding_event = rand() > .5; %change this to exponential pdf
%                           pols{i}(k) = pols{i}(k) + mcp_binding_event;
%                       end
%                   end
%               end
%               loading_success = rand() > .5; %change this to exponential pdf
%               if loading_success 
%                   pols{end+1} = [];                  
%               end
%               
              
              
        elseif f_init + f_elon < f && f < f_fin
            pol = pol + loadingRate - offLoadingRate;
            %do nothing. this is the trapezoid plateau
        elseif f >= f_fin
                pol = pol - offLoadingRate; %period of fluorescence decrease
        end
        
        if pol < 0
            pol = 0;
        end
        
        %NOTE FOR AR 8/30/2018- make this calibration more complicated by
        %adding a poissonian random variable that decides if an mcp lands
        %on a loop and converts pol to fluo. 
        
        fluo = pol*params.FLUO_CALIBRATION + scopeNoise + backGround - params.MEAN_BACKGROUND;
        
        if fluo < 0 
            fluo = 0;
        end
        
        pols(f) = pol;
        fluos(f) = fluo;
       
    end
    
%     fluosMask = fluos==0; %keep track of where the zeroes were
    singleFluosTrunc = fluos.*(fluos>params.DETECTION_LIMIT); %set values below limit to zero
    singleFluosTrunc(singleFluosTrunc==0) = NaN; %make those values above NaNs to avoid averaging them
%     singleFluosTrunc(fluosMask == 1) = 0; %get back the original zeros for proper averaging
%     if sum(isnan(singleFluosTrunc)) == length(singleFluosTrunc)
%         disp 'lost trace'
%     end

    if displaySingle
        figure()
        plot(timeElapsed, fluos)
        xlim([0, timeElapsed(end)])
        ylim([0, max(fluos)*1.20])
        ylabel('fluorescence (a.u.)')
        xlabel('time since division (min)')
        title('single trace- detection unlimited')
        standardizeFigure(gca, [])

        figure()
        plot(timeElapsed, singleFluosTrunc)
        xlim([0, timeElapsed(end)])
        ylim([0, max(fluos)*1.20])
        ylabel('fluorescence (a.u.)')
        xlabel('time since division (min)')
        title('single trace- detection limited')
        standardizeFigure(gca, [])
        pause(3) %leave 6 seconds to view the figure
    end
    
    
end

function [intFluo, intMeanOfTrunc] =averageSimulatedTraces(params, displaySingle, displayAverage)
   
   nTraces = params.nTraces;
   fractionActive = params.FRACTION_ACTIVE;
   detectionLimit = params.DETECTION_LIMIT;
   
   nActive = round(fractionActive*nTraces);
   reportedActive = nActive;
   
   times = [];
   fluos = [];
   
   for i = 1:nActive
       [times(i,:), fluos(i,:), singleFluosTrunc(i,:)] = simulateTrace(params, displaySingle);
       if sum(isnan(singleFluosTrunc(i,:))) == length(singleFluosTrunc(i,:))
           reportedActive = reportedActive - 1;
       end
   end
   
   reportedFractionActive = reportedActive/nTraces;
 
   timeElapsed = times(1,:); %the time elapsed is the same for every row so just grab the first
      
   %get the detection-unlimited information
   meanFluo = nanmean(fluos,1);
   intFluo= nansum(fluos(:));
   
   %get the detection limited traces
   meanOfTrunc = nanmean(singleFluosTrunc, 1);
   intMeanOfTrunc = nansum(meanOfTrunc(:));
   
   %truncate the mean of the detection unlimited traces
   meanFluoMask = meanFluo==0;
   truncOfMean = meanFluo.*(meanFluo>detectionLimit);
   truncOfMean(truncOfMean==0) = NaN;
%    truncOfMean(meanFluoMask == 1) = 0;
   intTruncOfMean =  nansum(truncOfMean(:));
   
   
   if displayAverage
       figure()
       plot(timeElapsed, meanFluo)
       xlim([0, timeElapsed(end)])
       ylim([0, max(meanFluo)*1.20])
       ylabel('mean fluorescence (a.u.)')
       xlabel('time since division (min)')
       title('mean trace- detection unlimited')
       legend(['cumulative fluorescence: ', num2str(intFluo)])
       standardizeFigure(gca, [])

       figure()
       plot(timeElapsed, meanOfTrunc)
       xlim([0, timeElapsed(end)])
       try
        ylim([0, max(meanOfTrunc)*1.20])
       catch
       end
       ylabel('mean fluorescence (a.u.)')
       xlabel('time since division (min)')
       title('mean of detection-limited traces')
       legend(['cumulative fluorescence: ', num2str(intMeanOfTrunc)])
       standardizeFigure(gca, [])

       figure()
       plot(timeElapsed, truncOfMean)
       xlim([0, timeElapsed(end)])
       try
        ylim([0, max(truncOfMean)*1.20])
       catch
           %
       end
       ylabel('mean fluorescence (a.u.)')
       xlabel('time since division (min)')
       title('detection-limited mean of traces')
       legend(['cumulative fluorescence: ', num2str(intTruncOfMean)])
       standardizeFigure(gca, [])

       figure()
       plot(timeElapsed, meanFluo)
       hold on
       plot(timeElapsed, meanOfTrunc)
       hold on
       plot(timeElapsed, truncOfMean)
       hold on
       plot(timeElapsed, detectionLimit*ones(length(timeElapsed)));
       hold off
       xlim([0, timeElapsed(end)])
       try
         ylim([0, max(meanOfTrunc)*1.20])
       catch
       end
       ylabel('mean fluorescence (a.u.)')
       xlabel('time since division (min)')
       title('effect of detection limit on trace averaging')
       legend(['mean of detection-unlimited traces \newline cumulative fluorescence: ', num2str(intFluo)],...
           ['mean of detection-limited traces \newline cumulative fluorescence: ', num2str(intMeanOfTrunc),...
           '\newline fraction active measured/true fraction active: ', num2str(reportedFractionActive),'/',num2str(fractionActive)],...
           ['detection-limited mean of traces \newline cumulative fluorescence: ', num2str(intTruncOfMean)],...
           'detection limit')
       standardizeFigure(gca, [])
       leg = findobj(gcf, 'Type', 'Legend');
       leg.FontSize = 11;
       lines = findobj(gcf, 'Type', 'Line');
       for i = 1:length(lines)
           lines(i).LineWidth = 3;
           lines(i).Marker = '.';
           lines(i).MarkerSize = 12;
       end
   
    end
   
end

