function HistD = AutoDoGThreshold(Prefix, InitialThreshold, SkipASLFlag, SkipFinalAnalysisFlag, NC_R, LogXBinC, XTickR)

% The idea is to automatically determine the optimum threshold to used
% based on the separation between legitimate transcription spots and noise.


%% Information about about folders and parameters
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);


%% More parameters
if ~exist('InitialThreshold', 'var') || isempty(InitialThreshold), InitialThreshold = 30; end
% [[[Optimize InitialThreshold once you've run the script on a couple of
% sets]]]
if ~exist('SkipASLFlag', 'var') || isempty(SkipASLFlag), SkipASLFlag = 1; end
if ~exist('SkipFinalAnalysisFlag', 'var') || isempty(SkipFinalAnalysisFlag), SkipFinalAnalysisFlag = 0; end
if ~exist('LogXBinC', 'var') || isempty(LogXBinC)
    LogXBinC = (linspace(log(InitialThreshold), log(10^2 * InitialThreshold), 32))';
end
XBinC = exp(LogXBinC);
if ~exist('NC_R', 'var') || isempty(NC_R), NC_R = [12 13 14]; end
% This is the list of which nuclear cycles to investigate
if ~exist('XTickR', 'var') || isempty(XTickR)
    XTickR = [10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 80, 90,...
        100, 120, 140, 160, 180, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000];
end
% These are the spot brightnesses that the code will bin over.
MaxDoGColorR = [0.5 0.5 1];
MaxDoGFitColorR = [0 0 0.75];
MaxThirdColorR = [1 0.5 0.5];
MaxThirdFitColorR = [0.75 0 0];


%% Initial analysis run

% We are using the FindAPAxis and nuclear tracking information from the
% original data set.

if ~SkipASLFlag
    cd([FISHPath, filesep, '\Analysis']);
    analyzeShawnLibrary('fad',@(x)tagged(x,'id',[Prefix, '_']),'params_mRNADynamics', InitialThreshold)
end

if ~SkipFinalAnalysisFlag
    
    % Delete existing files
    if exist([DropboxFolder,filesep,Prefix, filesep, 'Particles.mat'], 'file')
        delete([DropboxFolder,filesep,Prefix, filesep, 'Particles.mat']);
    end
    if exist([DropboxFolder,filesep,Prefix, filesep, 'CompiledParticles.mat'], 'file')
        delete([DropboxFolder,filesep,Prefix, filesep, 'CompiledParticles.mat']);
    end
    
    % Track the particles
    fprintf('\nTracking %s at a threshold of %d\n', Prefix, InitialThreshold);
    TrackmRNADynamics(Prefix, InitialThreshold, InitialThreshold);
    
    %Do the check. This is necesary for CompileParticles to run.
    CheckParticleTracking(Prefix, 'ForCompileAll');
    
    %Compile the particles. In this case we are approving all particles.
    CompileParticles(Prefix, 'ApproveAll', 'ForceAP', 'SkipTraces', 'SkipFluctuations', 'SkipFits');
    
    %Rename the Particles.mat file
    movefile([DropboxFolder,filesep,Prefix, filesep, 'Particles.mat'],...
        [DropboxFolder,filesep,Prefix, filesep, 'Particles-',num2str(InitialThreshold),'.mat']);
    
    %Rename the CompiledParticles.mat file
    movefile([DropboxFolder,filesep,Prefix, filesep, 'CompiledParticles.mat'],...
        [DropboxFolder,filesep,Prefix, filesep, 'CompiledParticles-',num2str(InitialThreshold),'.mat']);
end


%% Plotting and optimum threshold determination

% Loading particle data
DataDoG = load([DropboxFolder,filesep,Prefix, filesep, 'CompiledParticles-',num2str(InitialThreshold),'.mat']);
DataDoGOriginal = load([DropboxFolder,filesep,Prefix, filesep, 'Particles-',num2str(InitialThreshold),'.mat']);
MaxDoG_RCA{14} = [];
MaxThirdRCA{14} = [];
MeanAPPosRCA{14} = [];
% The entries of MaxDoG_RCA and MaxThirdRCA are nuclear cycles.
for lNC = NC_R
    [MaxDoG_RCA{lNC}, MaxThirdRCA{lNC}, MeanAPPosRCA{lNC}] =...
        BrightnessAndThreshold_FindMaxDoGAndThirds(DataDoG, DataDoGOriginal, lNC);
end

FigHandle = figure('WindowStyle', 'docked');
if find(NC_R == 12)
    Handle = subplot(2, 3, 1);
    PlotNCAllBins(Handle, 'nc12', MaxDoG_RCA{12}, 'Brightest slice', MaxDoGColorR, MaxDoGFitColorR);
    Handle = subplot(2, 3, 4);
    ThirdNC12_C =...
        PlotNCAllBins(Handle, 'nc12', MaxThirdRCA{12}, 'Dimmest adjacent slice', MaxThirdColorR, MaxThirdFitColorR);
else
    ThirdNC12_C = [];
end
if find(NC_R == 13)
    Handle = subplot(2, 3, 2);
    PlotNCAllBins(Handle, 'nc13', MaxDoG_RCA{13}, 'Brightest slice', MaxDoGColorR, MaxDoGFitColorR);
    Handle = subplot(2, 3, 5);
    ThirdNC13_C =...
        PlotNCAllBins(Handle, 'nc13', MaxThirdRCA{13}, 'Dimmest adjacent slice', MaxThirdColorR, MaxThirdFitColorR);
else
    ThirdNC13_C = [];
end
if find(NC_R == 14)
    Handle = subplot(2, 3, 3);
    PlotNCAllBins(Handle, 'nc14', MaxDoG_RCA{14}, 'Brightest slice', MaxDoGColorR, MaxDoGFitColorR);
    Handle = subplot(2, 3, 6);
    ThirdNC14_C =...
        PlotNCAllBins(Handle, 'nc14', MaxThirdRCA{14}, 'Dimmest adjacent slice', MaxThirdColorR, MaxThirdFitColorR);
else
    ThirdNC14_C = [];
end
if ~exist([DropboxFolder, filesep, Prefix, filesep, 'Threshold'], 'dir'),
    mkdir([DropboxFolder, filesep, Prefix, filesep, 'Threshold']);
end
saveas(FigHandle, [DropboxFolder, filesep, Prefix, filesep, 'Threshold', filesep, 'AutoDoGThreshold.fig']);
saveas(FigHandle, [DropboxFolder, filesep, Prefix, filesep, 'Threshold', filesep, 'AutoDoGThreshold.png']);
    
HistD = struct('XBinC', XBinC, 'ThirdNC12_C', ThirdNC12_C, 'ThirdNC13_C', ThirdNC13_C, 'ThirdNC14_C', ThirdNC14_C);

%% Plotting histograms as a function of AP bin
if find(NC_R == 12)
    PlotNCOneBin('nc12', MeanAPPosRCA{12}, MaxThirdRCA{12}, 'Dimmest adjacent slice', MaxThirdColorR);
end
if find(NC_R == 13)
    PlotNCOneBin('nc13', MeanAPPosRCA{13}, MaxThirdRCA{13}, 'Dimmest adjacent slice', MaxThirdColorR);
end
if find(NC_R == 14)
    PlotNCOneBin('nc14', MeanAPPosRCA{14}, MaxThirdRCA{14}, 'Dimmest adjacent slice', MaxThirdColorR);
end


%% Final analysis run

% {{{}}}


%% FitGauss2
    function [GaussLowC, GaussHighC, Crossing] = FitGauss2(HistC)
        
        FitTwoGauss = fit(LogXBinC, HistC, 'gauss2');
        
        %Find which one is to the right and left
        if FitTwoGauss.b1>FitTwoGauss.b2
            aRight=FitTwoGauss.a1;
            bRight=FitTwoGauss.b1;
            cRight=FitTwoGauss.c1;
            
            aLeft=FitTwoGauss.a2;
            bLeft=FitTwoGauss.b2;
            cLeft=FitTwoGauss.c2;
        else
            aRight=FitTwoGauss.a2;
            bRight=FitTwoGauss.b2;
            cRight=FitTwoGauss.c2;
            
            aLeft=FitTwoGauss.a1;
            bLeft=FitTwoGauss.b1;
            cLeft=FitTwoGauss.c1;
        end
        
        GaussLowC = aLeft*exp(-((LogXBinC-bLeft)/cLeft).^2);
        GaussHighC = aRight*exp(-((LogXBinC-bRight)/cRight).^2);
        DifferenceC = abs(GaussHighC-GaussLowC);
        
        [~, iGaussLowMax] = max(GaussLowC);
        [~, iGaussHighMax] = max(GaussHighC);
        DifferenceC(1:iGaussLowMax-1) = NaN;
        DifferenceC(iGaussHighMax+1:end) = NaN;
        
        [~, iCrossing] = nanmin(DifferenceC);
        LogCrossing = LogXBinC(iCrossing);
        Crossing = exp(LogCrossing);
        
    end


%% PlotNCAllBins
    function HistC = PlotNCAllBins(PlotHandle, TitleS, BrightnessesR, SetS, BarColorR, LineColorR)
        
        hold(PlotHandle, 'all');
        
        HistC = (hist(log(BrightnessesR), LogXBinC))';
        BarHandle = bar(XBinC, HistC, 'histc');
        set(BarHandle, 'EdgeColor', 'none', 'FaceColor', BarColorR);
        h = findobj(gca,'Type','line');
        set(h,'Marker','none'); 
        
        [LowC, HighC, Crossing] = FitGauss2(HistC);
        plot(XBinC, [LowC, HighC], 'Color', LineColorR, 'LineWidth', 2);        
        
        title({[SetS, ', ', TitleS], sprintf('Threshold at crossing: %02.1f', Crossing)});
        %legend({'Brightest slide', 'Third-brightest slice'}, 'Location',
        %'SouthOutside');
        xlim([XBinC(1), XBinC(end)]);
        set(gca, 'FontSize', 8);
        set(gca, 'XScale', 'log');
        set(gca, 'XTick', XTickR);
        xlabel({'Threshold (counts)'});
        ylabel('Number of nuclear spots');
        YLimR = get(gca, 'YLim');
        set(gca, 'YLim', [-YLimR(2)/6, YLimR(2)]);
        %legend(SetS);
        
    end


%% PlotNCOneBin
    function PlotNCOneBin(TitleS, ParticleAPBinR, BrightnessesR, SetS, BarColorR)
        APBinR = 0:0.1:1;
        
        [APBinCountsR, APBinIndexR] = histc(ParticleAPBinR, APBinR);
        APBinHasNucleiLR = (APBinCountsR > 0);
        APBinHasNucleiR = find(APBinHasNucleiLR);
        NumPopulatedBins = sum(APBinHasNucleiLR);
        NumRows = ceil(sqrt(NumPopulatedBins));
        
        figure('WindowStyle', 'docked');
        for lPopulatedBin = 1:NumPopulatedBins
            subplot(NumRows, NumRows, lPopulatedBin);
            hold on;
            ParticlesInBinLR = (APBinIndexR == APBinHasNucleiR(lPopulatedBin));
            HistC = (hist(log(BrightnessesR(ParticlesInBinLR)), LogXBinC))';
            BarHandle = bar(XBinC, HistC, 'histc');
            set(BarHandle, 'EdgeColor', 'none', 'FaceColor', BarColorR);
            h = findobj(gca,'Type','line');
            set(h,'Marker','none'); 
            
            BinMin = APBinR(APBinHasNucleiR(lPopulatedBin));
            BinMax = APBinR(APBinHasNucleiR(lPopulatedBin) + 1);
            title({[SetS, ', ', TitleS, ', '], sprintf('%.1f to %.1f%%', 100*BinMin, 100*BinMax)});
            xlim([XBinC(1), XBinC(end)]);
            set(gca, 'FontSize', 8);
            set(gca, 'XScale', 'log');
            %set(gca, 'XTick', XTickR);
            xlabel({'Threshold (counts)'});
            ylabel('Number of nuclear spots');
            YLimR = get(gca, 'YLim');
            set(gca, 'YLim', [-YLimR(2)/6, YLimR(2)]);
        end
        
        saveas(gcf, [DropboxFolder, filesep, Prefix, filesep, 'Threshold', filesep,...
            'AutoDoGThreshold_', TitleS, '.fig']);
        saveas(gcf, [DropboxFolder, filesep, Prefix, filesep, 'Threshold', filesep,...
            'AutoDoGThreshold_', TitleS, '.png']);
        
    end

end