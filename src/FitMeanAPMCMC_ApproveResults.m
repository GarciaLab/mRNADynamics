function FitMeanAPMCMC_ApproveResults(varargin)
%Last updated: 7/11/19 by Jonathan Liu

%Analyzes saved MCMC results from FitMeanAPMCMC. The user has the option of
%approving or rejecting the results of each single nucleus fit. The
%approved/rejected results are saved in the same .mat file as the MCMC results.

%The MCMC results should be saved as a .mat file with 2 stuctures, MCMCplot
%and MCMCresults. Each should be a cell array of size 2, containing a
%structure of size 1xN, where N is the number of AP bins in the dataset.
%Each cell corresponds to a different nuclear cycle. Additionally, 
%each AP bin possesses a structure field called ApprovedFits. The values of ApprovedFits correspond
%to as follows:
%   1: approved
%   0: uncurated (can be approved if ApproveAll is used in post-analysis)
%   -1: rejected

%If desired, this script can also load the raw MCMC chains for detailed
%visualization. The raw chains are saved in a separate .mat file containing
%the structure MCMCchain.

%Variable input arguments:
%   Prefix: Prefix string. If none chosen, user has the option to select
%           using a dialog menu.
%   'RawChains': visualize raw chains from MCMC inference.

%% Input arguments
LoadPrefix = true; %By default, user selects which Prefix to load.
RawChains = false; %By default, don't look at raw chains.
for i=1:length(varargin)
    if strcmpi(varargin{i},'Prefix')
        Prefix = varargin{i+1};
        LoadPrefix = false;
    end
    if strcmpi(varargin{i},'RawChains')
        RawChains = true;
    end
end

%% Load results
%Get the default folders for Prefix loading
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders;

if LoadPrefix
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end

%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);
        
%Load MCMC results                                
m = load([DropboxFolder,filesep,Prefix,'\MeanFitsMCMC.mat']);
N_bins = length(m.MCMCplot{1}); %Number of AP bins in this dataset
N_cycles = length(m.MCMCplot); %Number of nuclear cycles in this dataset

%Extract the data into temporary variables.
MCMCplot = m.MCMCplot;
MCMCresults = m.MCMCresults;
if RawChains
    chains = load([DropboxFolder,filesep,Prefix,'\MeanFitsMCMC_RawChains.mat']);
    MCMCchain = chains.MCMCchain;
end

%% Approve/reject fits

%User keypress options
InputKey = {'a';'r';',';'.';'j';'n';'m';'x'};
Function = {'Approve this fit';'Reject this fit';'Previous AP bin';...
    'Next AP bin';'Jump to a specific AP bin';'Previous nuclear cycle';...
    'Next nuclear cycle';'Exit and save'};
t = table(InputKey,Function);

%Figure color map
colormap = {'red',[0.94 0.94 0.94],'green'}; %Approved/Uncurated/Rejected colors

running = true; %Keep running program until stopped

%Find first nc and AP bin that has data
i = 1; %Start with first indexed AP bin that has data
nc = 1; %Start with the first indexed nuclear cycle for now

firstflag = false;
while ~firstflag
    if ~isempty(MCMCplot{nc}(i).t_plot)
        firstflag = true;
    else
        if i < length(MCMCplot{nc})
            i = i+1;
        else
            if nc < length(MCMCplot)
                i = 1;
                nc = nc+ 1;
            else
                warning('These MCMC results appear to be empty.');
            end
        end
    end
end


close all
if RawChains
    f = figure('Name','Inference Results','Position',[50 100 500 600]);
    c = figure('Name','Parameter distributions','Position',[600 100 600 600]);
else
    f = figure('Name','Inference Results','Position',[200 100 800 600]);
    
end

while running
    
clf(f); %Clear figures
clf(c);

%Extract plotting variables from results
t_plot = MCMCplot{nc}(i).t_plot;
MS2_plot = MCMCplot{nc}(i).MS2_plot;
t_interp = MCMCplot{nc}(i).t_interp;
MS2_interp = MCMCplot{nc}(i).MS2_interp;
simMS2 = MCMCplot{nc}(i).simMS2;
nuclearcycle = MCMCplot{nc}(i).nc;

%Check to see if this AP bin/nuclear cycle has any results
emptyflag = false;
if isempty(t_plot)
    emptyflag = true;
end

%Extract inference results
APPos = MCMCresults{nc}(i).APPosition;
mean_R0 = MCMCresults{nc}(i).mean_R0;
sigma_R0 = MCMCresults{nc}(i).sigma_R0;
mean_dR = MCMCresults{nc}(i).mean_dR;
sigma_dR = MCMCresults{nc}(i).sigma_dR;
mean_ton = MCMCresults{nc}(i).mean_ton;
sigma_ton = MCMCresults{nc}(i).sigma_ton;
mean_dwelltime = MCMCresults{nc}(i).mean_dwelltime;
sigma_dwelltime = MCMCresults{nc}(i).sigma_dwelltime;
mean_basalfluor = MCMCresults{nc}(i).mean_MS2_basal;
sigma_basalfluor = MCMCresults{nc}(i).sigma_MS2_basal;

%Smoothed loading rate
rate_plot = mean_R0 + mean_dR;
rateerror_plot = sigma_R0 + sigma_dR;
span = 0.1;
ratesmooth_plot = smooth(rate_plot,span);

%Remove loading rates before inferred initiation time
remove_times = find(t_interp(1:end-1) < mean_ton);
rate_plot(remove_times) = nan;
rateerror_plot(remove_times) = nan;
ratesmooth_plot(remove_times) = nan;

%Raw chain results
if RawChains
    dwelltime_chain = MCMCchain{nc}(i).dwelltime_chain; %termination dwell time
    R0_chain = MCMCchain{nc}(i).R0_chain; %mean loading rate
    dR_chain = MCMCchain{nc}(i).dR_chain; %last loading rate fluctuation 
end

%% Plot fit results
figure(f);

subplot(2,1,1)
if ~emptyflag
    hold on
    plot(t_plot,MS2_plot,'go','DisplayName','Fluorescence data');
    plot(t_interp,MS2_interp,'g-','DisplayName','Fluorescence interpolation');
    plot(t_interp,simMS2,'k--','DisplayName','MCMC fit');
    hold off

    xlim([t_interp(1), t_interp(end)*1.3]);
    %ylim([0,  max(MS2_plot) * 1.2]);
    line(xlim,mean_basalfluor*[1,1],'Color','magenta','LineStyle','--',...
        'DisplayName','Inferred basal fluorescence');
    line(mean_ton*[1,1],ylim,'Color','blue','LineStyle','--',...
        'DisplayName','Inferred time on');
end
xlabel('Time since nuclear cycle start (min)');
ylabel('Fluorescence (AU)');
legend('Location','southeast');
title({['MCMC fit: AP Position ',num2str(APPos), ', Nuclear cycle ',num2str(nuclearcycle)]...
    ['Termination dwell time = ',num2str(mean_dwelltime), ' +/- ',num2str(sigma_dwelltime), 'min']});
f.Color = colormap{MCMCresults{nc}(i).ApprovedFits+2}; %Set color depending on approval status

%Plot inferred rate
if ~emptyflag
    subplot(2,1,2);
    hold on
    errorbar(t_interp(1:end-1),rate_plot,rateerror_plot,'r.-','CapSize',0,...
        'DisplayName','Inferred loading rate');
    plot(t_interp(1:end-1),ratesmooth_plot,'k--','DisplayName','Smoothed loading rate');
    line(xlim,mean_R0(1)*[1,1],'Color','black','LineStyle','-','DisplayName',...
        'Inferred mean loading rate');
    hold off

    xlim([t_interp(1), t_interp(end)*1.3]);
    %ylim([0,  max(MS2_plot) * 1.2]);
    line(mean_ton*[1,1],ylim,'Color','blue','LineStyle','--',...
        'DisplayName','Inferred time on');
end
    xlabel('Time since nuclear cycle start (min)');
    ylabel('Loading rate (AU/min)');
    legend('Location','southeast');

%Plot raw chains if desired
if RawChains
    if ~emptyflag
        figure(c);
        hold on;

        subplot(3,2,1)
        histogram(dwelltime_chain);
        xlabel('Dwell time (min)');

        subplot(3,2,2)
        plot(dwelltime_chain,'b.');
        ylabel('Dwell time(min)');

        subplot(3,2,3)
        histogram(R0_chain);
        xlabel('Mean loading rate (AU/min)');

        subplot(3,2,4)
        plot(R0_chain,'b.');
        ylabel('Mean loading rate (AU/min)');

        subplot(3,2,5)
        histogram(dR_chain(:,end));
        xlabel('Last loading rate fluctuation (AU/min)');

        subplot(3,2,6)
        plot(dR_chain(:,end),'b.');
        ylabel('Last loading rate fluctuation (AU/min)');
    end
end
%% User options (approve/reject, change nucleus)
disp(t); %Display options
exitflag = false; %Loop keypress query until valid exit keypress

while ~exitflag
figure(f);
waitforbuttonpress; %User input to press a key
key = f.CurrentCharacter; %Last pressed key

if strcmp(key,'a')
    MCMCresults{nc}(i).ApprovedFits = 1;
    disp('Approved');
elseif strcmp(key,'r')
    MCMCresults{nc}(i).ApprovedFits = -1;
    disp('Rejected');
elseif strcmp(key,',')
    if i > 1
        i = i - 1;
        disp('Switching to previous AP position');
        exitflag = true;
    elseif i == 1
        disp('Already at first inferred AP position!');
    end
elseif strcmp(key,'.')
    if i < N_bins
        i = i + 1;
        disp('Switching to next AP position');
        exitflag = true;
    elseif i == N_bins
        disp('Already at last inferred AP position!');
    end
elseif strcmp(key,'n')
    if nc > 1
        nc = nc - 1;
        disp('Switching to previous nuclear cycle');
        exitflag = true;
    elseif nc == 1
        disp('Already at first inferred nuclear cycle!');
    end
elseif strcmp(key,'m')
    if nc < N_cycles
        nc = nc + 1;
        disp('Switching to next nuclear cycle');
        exitflag = true;
    elseif nc == N_cycles
        disp('Already at last inferred nuclear cycle!');
    end
elseif strcmp(key,'j')
    j = input('Enter in AP bin index to jump to:');
    if j >= 1 && j <= N_bins
        i = j;
        disp(['Switching to AP bin ',num2str(i)]);
        exitflag = true;
    else
        disp(['Error: please enter an integer between 1 and ',num2str(N_bins)]);
    end
elseif strcmp(key,'x')
    disp('Exiting and saving results')
    exitflag = true;
    running = 0;
end
f.Color = colormap{MCMCresults{nc}(i).ApprovedFits+2}; %Set color depending on approval status
end

end

%Save approved fits results
m.MCMCresults = MCMCresults;
save([DropboxFolder,filesep,Prefix,'\MeanFitsMCMC.mat'],'MCMCresults','MCMCplot','Prefix');

disp('Results saved.');

close all;

end