

set(0, 'DefaultAxesFontSize', 20)

% Define the dynamic range of the offseted DV axis
MaxFakeDV = ceil(max([schnitzcells.DVpos]));
MinFakeDV = floor(min([schnitzcells.DVpos]));
TotalFrames = length(FrameInfo);
PixelSize = FrameInfo(1).PixelSize;

% Make an empty 2 x N x T empty array to fill with each frame's data.
% N = maximum number of nuclei per frame. Take from CompiledNuclei ->
% NParticlesAll
% T = time in frames.
% the first dimension has a size of two because it's a collection of at
% most N pairs of values: DVpos and Fluo
DlFluoGaussians = nan(2,nanmax(NParticlesAll),TotalFrames);
Start = 1;
End = TotalFrames;

% extract nuclear fluorescence
for f = Start:End
    f/length([Start:End])
    %waitbar(f/length([Start:End]),['frame #' f]);
    count = 0;
    for s = 1:length(schnitzcells)        
        schnitzFrames = schnitzcells(s).frames;
        schnitzFluo = ExtractDlFluo(schnitzcells(s).Fluo,0.5);
        CurrentFrame = find(schnitzFrames==f); %find the position in the schnitzFrames vector of the current frame
        if CurrentFrame
            count = count+1;
            CurrentDV = round(schnitzcells(s).DVpos,0);
            CurrentDV = CurrentDV(CurrentFrame);
            CurrentFluo = schnitzFluo(CurrentFrame);
            FluoGaussians(1,count,f) = CurrentDV; 
            FluoGaussians(2,count,f) = CurrentFluo;
        end
    end
    'fin'
end
%% Visualizing
for f = 1:End
    plot(FluoGaussians(1,:,f),FluoGaussians(2,:,f),'o','MarkerFaceColor',[.8 .7 1]);
    ylim([300,3500])
    xlim([MinFakeDV,MaxFakeDV])
    title(num2str(f))
    waitforbuttonpress
end

%% Fitting
for frame = 1:TotalFrames%[119:137]
    DataX = FluoGaussians(1,:,frame);
    %DataX = DataX(DataX ~= 0);
    DataY = FluoGaussians(2,:,frame);
    %DataY = DataY(DataX ~= 0);
    LastPoint = find(DataY,1,'last'); %this is necessary to get rid of zeros
    DataX = DataX(1:LastPoint)
    DataY = DataY(1:LastPoint)
    
    if length(DataY) > 10
          
        % Sort based on DV position
        DataXY(1,:) = DataX;
        DataXY(2,:) = DataY;
        SortedDataXY = sortrows(DataXY', 1);
        DataX = SortedDataXY(:,1);
        DataY = SortedDataXY(:,2);
        SmoothDataY = smooth(DataY,5);
        clear DataXY
        
        % Interpolate over NaNs 
        DataY(isnan(DataY)) = interp1(find(~isnan(DataY)), DataY(~isnan(DataY)), find(isnan(DataY)), 'cubic');         
        
        % now do fitting stuff
        Gaussfun = @(z)z(1) * exp(-(DataX-z(2)).^2./(2*(z(3)^2))) + z(4) - SmoothDataY;

        % make guesses about z(1), z(2), z(3) and z(4) (amplitude, center, spread and baseline)
        z0 = [2000,-200,200,500];
        % define upper and lower bounds for each
        lb = [200,-550,50,100];
        ub = [5000,550,500,700];

        Guess = lsqnonlin(Gaussfun,z0,lb,ub);

        FitGauss = Guess(1)*exp(-(DataX-Guess(2)).^2./(2*(Guess(3)^2))) + Guess(4) ;

        Centers(frame) = Guess(2) * PixelSize;

        %Error(rep) = round((b-Guess(2))/c,3); %metric for how far the guessed center is in terms of one standard deviation
        
        figure(1)
        plot(DataX * PixelSize,DataY,'o','MarkerFaceColor',[1 .8 .8],'MarkerEdgeColor',[1 .8 .8])
        hold on
        plot(DataX * PixelSize,SmoothDataY,'o','MarkerFaceColor',[1 .3 .3])%,'LineWidth',1.8)
        plot(DataX * PixelSize,FitGauss,'-','Color',[.6 .6 1],'LineWidth',3)
        plot([Guess(2)* PixelSize Guess(2)* PixelSize],[300,3500],'k')
        hold off
        ylim([300,5500])
        xlim([MinFakeDV* PixelSize,MaxFakeDV* PixelSize])
        title(['frame' num2str(frame)])
        legend('Dl-Venus nuclear fluo (raw)','Dl-Venus nuclear fluo (smoothened)','Gaussian fit','Fitted center')

        %waitforbuttonpress
    end
end

%
%Centers = Centers * PixelSize; %convert from pixels to microns
figure(2)
subplot(1,2,1)
hist(Centers,TotalFrames/3)
title('Gaussian centers')
xlabel('DV distance (\mum)')
ylabel('counts')
subplot(1,2,2)
plot(Centers,'o')
hold on
plot([0,length(Centers)],[mean(Centers(Centers~=0)),mean(Centers(Centers~=0))],'k-')
plot([0,length(Centers)],[median(Centers(Centers~=0)),median(Centers(Centers~=0))],'r-')
ylabel('DV distance (\mum)')
xlabel('frame')
hold off
legend('fitted centers','mean','median')

    