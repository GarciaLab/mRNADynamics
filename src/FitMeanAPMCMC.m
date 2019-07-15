function FitMeanAPMCMC(varargin)

%Last updated on: 7/13/19 by Jonathan Liu

%
%
%Fits constant elongation model with dwell time to 1-color data using MCMCstat package. The script
%loads a single analyzed dataset and finds a probability distribution for the parameters
%using Markov Chain Monte Carlo, for each MeanVectorAP over AP positions
%over nuclear cycles 13 and/or 14.

%The main relevant parameters are the time-dependent RNAP loading rate and
%the termination dwell time of each polymerase (i.e. the time the polymerase
%spends on after elongating the gene). The elongation rate is fixed, and
%can be specified by the user.

%The initial parameter values are randomly chosen
%within a reasonable range. This script assumes a fixed elongation rate and
%and basal activity levels (i.e. noise floor) for each signal. The inferred
%parameters, along with relevant dataset descriptions, are saved in a .mat
%file for each dataset, called MeanFitsMCMC.mat, to be used for later
%analysis.

%For inferring the loading rate as a function of time, the script assumes a
%constant mean loading rate and allows for time-dependent fluctuations, but
%punishes the fluctuations if they deviate too far from zero. The scale of
%deviation (i.e. Bayesian prior) is hardcoded, but can also be set by the user.

%The user has the option of specifying certain properties, or can use the
%default settings. These properties are saved in the parameter log
%structure param_log, inside the MeanFitsMCMC.mat file.

%Variable inputs (Name/Value pairs):
%   'Prefix', Prefix: Input a Prefix string. If not chosen, user
%                           selects Prefix in a dialog box.
%   'MCMCsteps', [n_steps, n_burn]: MCMC step/burn-in counts. Default is
%                           n_steps = 30000, n_burn = 15000.
%   'NCWindow', [nc13start, nc13end; nc14start, nc14end]: Start and end
%                           times of fitting for cycles 13 and 14.
%                           Pass this in as a 2x2 array of doubles.
%                           Read as [time after nc13 start, time after nc14 start, time
%                           after nc14 start, time after nc14 start].
%                           Default values are [1.5, -5; 1.5, 18].
%   'Construct', construct: Name of construct used to get total length of
%                           gene. Definitions of constructs must first be
%                           given in FitMeanAPMCMC_GetFluorFromPolPos
%                           function. Default construct is 'P2P-MS2v5-LacZ'.
%   'ElongRate', v_elon: Elongation rate in kb/min. Default value is
%                           1.8kb/min.
%   'MeanRateDeviation', prior_sigma: set scale of loading rate fluctuations, where
%                           the fluctuations have a Gaussian prior centered
%                           around zero with variance s. Default value is
%                           10.
%   'ChooseNC', {ncstrings}: select which nuclear cycles to analyze by
%                            passing in a cell array of strings. Currently
%                            only 'nc13' and 'nc14' are supported.
%   'KeepPool', []: Keeps the parallel pool after running (no 2nd argument
%                   required)

%% Variable inputs

%Default settings
prior_sigma = 10; %Width of Gaussian prior in loading rate fluctuations
ncwindow = [1.5, -5; 1.5, 18]; %Time of MCMC fitting for nc13 and nc14
MCMCsteps = [30000 15000]; %Number of MCMC total steps and burn-in steps.
LoadPrefix = true; %User selection of Prefix by default.
construct = 'P2P-MS2v5-LacZ'; %Default construct.
v_elong = 1.8; %Default elongation rate.
ChooseNC = false; %By default, analyze both nc13 and nc14
KeepPool = false; %By default, shut down parallel pool after running code

%User specified settings
for i=1:length(varargin)
    if strcmpi(varargin{i},'Prefix')
        Prefix = varargin{i+1};
        LoadPrefix = false;
    end
    if strcmpi(varargin{i},'MCMCsteps')
        MCMCsteps = varargin{i+1};
    end
    if strcmpi(varargin{i},'NCWindow')
        ncwindow = varargin{i+1};
    end
    if strcmpi(varargin{i},'Construct')
        construct = varargin{i+1};
    end
    if strcmpi(varargin{i},'ElongRate')
        v_elong = varargin{i+1};
    end
    if strcmpi(varargin{i},'MeanRateDeviation')
        prior_sigma = varargin{i+1};
    end
    if strcmpi(varargin{i},'ChooseNC')
        ncstrings = varargin{i+1};
        ChooseNC = true;
    end
    if strcmpi(varargin{i},'KeepPool')
        KeepPool = true;
    end
end

%Extract MCMC settings and nuclear cycle fit windows.
n_steps = MCMCsteps(1);
n_burn = MCMCsteps(2);
firstnc13time = ncwindow(1,1);
lastnc13time = ncwindow(1,2);
firstnc14time = ncwindow(2,1);
lastnc14time = ncwindow(2,2);

%Save MCMC parameters into cell array
param_log.prior_sigma = prior_sigma;
param_log.ncwindow = ncwindow;
param_log.MCMCsteps = MCMCsteps;
param_log.construct = construct;
param_log.v_elong = v_elong;

%Extract which nuclear cycle(s) to analyze
if ChooseNC
    ncToAnalyze = [];
    for i = 1:length(ncstrings)
        if contains(ncstrings{i},'nc13')
            ncToAnalyze(end+1) = 1;
        elseif contains(ncstrings{i},'nc14')
            ncToAnalyze(end+1) = 2;
        end
    end
    ncToAnalyze = sort(ncToAnalyze); %Sort nc indices in ascending order
else
    ncToAnalyze = [1 2];
end

%Check to make sure n_burn and n_steps are more than zero
if n_steps < 1
    error('n_steps must be a positive integer.');
end
if n_burn < 1
    error('n_burn must be a positive integer.');
end

%% Load data
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
        
%Load the compiled particles and the division information                                    
data = load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat']);

% Extract the fields from the cell structure (This is for fields like MeanVectorAP
% that are saved inside {}.
channel = 1; %We're working with one-color data for now.

if iscell(data.MeanVectorAP)
    MeanVectorAP = data.MeanVectorAP{channel};
    SDVectorAP = data.SDVectorAP{channel};
else
    MeanVectorAP = data.MeanVectorAP;
    SDVectorAP = data.SDVectorAP;
end

%Figure out what type of experiment we have. Note: we name the var "DateFromDateColumn" to avoid shadowing previously defined "Date" var.
[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DropboxFolder);

%% Set some initial definitions

%Definitions
nc_string = {'nc13','nc14'}; %Cell array of nc strings
nc13 = data.nc13; %nc13 frame index
nc14 = data.nc14; %nc14 frame index
t = data.ElapsedTime; %time variable
dt = mean(t((nc14+1):end)-t(nc14:(end-1))); %Set a simulation timestep given by mean time resolution
scalefac = 10; %Arbitrary scale factor to reduce simulation size
bins = 0:APResolution:1; %AP resolution
N_bins = length(bins); %Total number of bins
N_cycles = 2; %Number of nuclear cycles

%% MCMC analysis loop through nuclear cycles and AP bins
nc_array = [nc13, nc14, length(t)]; %Array of relevant nuclear cycle info

w = waitbar(0,'Processing nuclear cycles');

for nc = ncToAnalyze
waitbar(nc/N_cycles,w,['Processing nuclear cycle ',num2str(nc),' of ',num2str(N_cycles)]);

%Setup data to be saved in structure. For parallel usage, save each nuclear
%cycle in a separate structure array, which we will combine into a master
%structure array later. Populate all the results with nan's for now, to
%preserve 1xN_bins structure.
clear MCMCchain_temp;
clear MCMCresults_temp;
clear MCMCplot_temp;

MCMCchain_temp(N_bins) = struct();
MCMCresults_temp(N_bins) = struct();
MCMCplot_temp(N_bins) = struct();

for AP = 1:N_bins
    MCMCchain_temp(AP).ton_chain = [];
    MCMCchain_temp(AP).MS2_basal_chain = [];
    MCMCchain_temp(AP).dwelltime_chain = [];
    MCMCchain_temp(AP).R0_chain = [];
    MCMCchain_temp(AP).dR_chain = [];
    MCMCchain_temp(AP).s2_chain = [];

    MCMCresults_temp(AP).mean_ton = nan;
    MCMCresults_temp(AP).sigma_ton = nan;
    MCMCresults_temp(AP).mean_MS2_basal = nan;
    MCMCresults_temp(AP).sigma_MS2_basal = nan;
    MCMCresults_temp(AP).mean_dwelltime = nan;
    MCMCresults_temp(AP).sigma_dwelltime = nan;
    MCMCresults_temp(AP).mean_R0 = nan;
    MCMCresults_temp(AP).sigma_R0 = nan;
    MCMCresults_temp(AP).mean_dR = nan;
    MCMCresults_temp(AP).sigma_dR = nan;
    MCMCresults_temp(AP).mean_sigma = nan;
    MCMCresults_temp(AP).sigma_sigma = nan;
    MCMCresults_temp(AP).APPosition = bins(AP);
    MCMCresults_temp(AP).ApprovedFits = 0; %Set approval status to uncurated (0).
    MCMCresults_temp(AP).nc = nc_array(nc);

    MCMCplot_temp(AP).t_plot = [];
    MCMCplot_temp(AP).MS2_plot = [];
    MCMCplot_temp(AP).t_interp = [];
    MCMCplot_temp(AP).MS2_interp = [];
    MCMCplot_temp(AP).simMS2 = [];
    MCMCplot_temp(AP).firstnc13time = [];
    MCMCplot_temp(AP).lastnc13time = [];
    MCMCplot_temp(AP).firstnc14time = [];
    MCMCplot_temp(AP).lastnc14time = [];
    MCMCplot_temp(AP).nc = nc_array(nc);
end

parfor AP = 1:N_bins
%Information about this AP bin
MCMCresults_temp(AP).APPosition = bins(AP); %AP position information
MCMCresults_temp(AP).ApprovedFits = 0; %Approval variable for later curation
MCMCresults_temp(AP).nc = nc_string{nc}; %Nuclear cycle info

%MS2 fluorescence data
MS2 = MeanVectorAP(:,AP)'./scalefac;

%Skip this AP bin if there are fewer than 10 datapoints in the nuclear
%cycle
if sum(~isnan(MS2(nc_array(nc):nc_array(nc+1)))) < 10
    disp(['skipping AP bin ',num2str(AP),' of ',num2str(N_bins), ' (no data)']);
    continue
end

%To regularize the fitting, we interpolate the traces to have constant time
%resolution, and to fill in nans.

%Find the frames for interpolation cutoff times for each nuclear cycle. For the
%first timepoint, find the time of the first non-nan datapoint. For the last
%timepoint, find the time of the last non-nan datapoint.

%First frame of interpolation (nc13)
firstnc13MS2 = find(and(t - t(nc13) > firstnc13time, ~isnan(MS2)),1);

%Last frame of interpolation (nc13)
lastnc13MS2 = find(and(t - t(nc14) > lastnc13time, ~isnan(MS2)),1);
if isempty(lastnc13MS2)
    lastnc13MS2 = find(~isnan(MS2(nc13:nc14)),1,'last'); %If no last nc13 time found, use last frame of MS2 data
end

%First frame of interpolation (nc14)
firstnc14MS2 = find(and(t - t(nc14) > firstnc14time, ~isnan(MS2)),1);

%Last frame of interpolation (nc14)
lastnc14MS2 = find(and(t - t(nc14) > lastnc14time, ~isnan(MS2)),1);
if isempty(lastnc14MS2)
    lastnc14MS2 = find(~isnan(MS2(nc14:end)),1,'last'); %If no last nc14 time found, use last frame of MS2 data
end

%Now, interpolate the data.
firstMS2 = [firstnc13MS2, firstnc14MS2];
lastMS2 = [lastnc13MS2, lastnc14MS2];

try
    t_interp = (t(nc_array(nc)):dt:t(lastMS2(nc))) - t(nc_array(nc));
catch
    warning(['Error in interpolating nc',num2str(nc+12),...
        ', AP position ', num2str(bins(AP)),'...skipping']);
    continue
end

%Fill in nans with pchip interpolation
MS2_pchip = fillmissing(MS2,'pchip');

try
    MS2_interp = interp1(t(nc_array(nc):lastMS2(nc))-t(nc_array(nc)),...
        MS2_pchip(nc_array(nc):lastMS2(nc)),t_interp);
catch
    warning(['Error in interpolating nc',num2str(nc+12),...
        ', AP position ', num2str(bins(AP)),'...skipping']);
    continue
end

%Define the subset of data to be analyzed. Analyze between the first
%and last datapoints in each nuclear cycle,
%but for the beginning cutoff leave the earlier times in, just with nan
%values. This is so we can still fit t_on, but we don't interpolate past
%the first datapoint.

if ~isempty(firstMS2(nc))
    MS2_interp(t_interp+t(nc_array(nc))<t(firstMS2(nc))) = nan;
end
if ~isempty(lastMS2(nc))
    MS2_interp(t_interp+t(nc_array(nc))>t(lastMS2(nc))) = nan;
end

%% Setup MCMC Fit

%Define the anonymous function for the sum of squares residual function.
ssfun = @(x,data) FitMeanAPMCMC_SumOfSquares(construct,v_elong,data,x);

% Initialize MCMC parameters with random values
ton0 = 4*rand; %Transcription turn on time (in min)
dwelltime0 = 4*rand+2; %Total dwell time (in min)
MS2_basal0 = 10; %Basal MS2 fluorescence (in AU)
R0 = 10+5*rand; %Mean loading rate (AU/min)
dR0 = normrnd(0,3,1,size(MS2_interp,2)-1); %Fluctuations in loading rate (AU/min)

x0 = [ton0,dwelltime0,MS2_basal0,R0,dR0]; %Initial parameter set

sigma2_0 = 1; %Initial guess for error variance

%This is the initial variance in the parameter proposal distribution. Change these
%to change the proposal acceptance rate, or if the convergence doesn't look
%good.
ton_step = dt;
dwelltime_step = 0.1;
MS2_basal_step = 1;
R_step = 0.5;
dR_step = 0.5*ones(size(dR0));

J0 = diag([ton_step,dwelltime_step,MS2_basal_step,R_step,dR_step]); %Initial covariance matrix

%% Setup MCMC parameters and options
%These parameters are expressed in the form:
%{parameter name, initial value, lower bound, upper boind}
params = {
    {'ton', x0(1), 0, 10}
    {'dwelltime', x0(2), 0, 20}
    {'MS2_basal', x0(3), 0, 50}
    {'R_0', x0(4), 0, 40}
    };

for i = 1:size(dR0,2)
    params{end+1} = {strcat('dR',num2str(i)), x0(4+i), -30, 30, 0, prior_sigma}; %Fluctuations with priors around zero
end

%MCMC data structure
data = [];
data.xdata = t_interp; %Time data
data.ydata = MS2_interp; %MS2 data

%Model definitions
model = [];
model.ssfun = ssfun;
model.sigma2 = sigma2_0;
model.N = length(data.ydata);

%MCMC options
options = [];
options.nsimu = n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
options.qcov = J0; %Initial covariance
options.burnintime = n_burn; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output

%Run the MCMC
[results,chain,s2chain] = mcmcrun(model,data,params,options);

%Extract chain results into individual parameters
ton_chain = chain(n_burn:end,1);
dwelltime_chain = chain(n_burn:end,2);
MS2_basal_chain = chain(n_burn:end,3);
R0_chain = chain(n_burn:end,4);
dR_chain = chain(n_burn:end,5:end);

%Change negative loading rates to zero
for i = 1:size(dR_chain,1)
    for j = 1:size(dR_chain,2)
        if dR_chain(i,j) + R0_chain(i) < 0
            dR_chain(i,j) = -R0_chain(i);
        end
    end
end

%Mean results
mean_ton = mean(ton_chain);
sigma_ton = std(ton_chain,1);
mean_dwelltime = mean(dwelltime_chain);
sigma_dwelltime = std(dwelltime_chain,1);
mean_MS2_basal = mean(MS2_basal_chain);
sigma_MS2_basal = std(MS2_basal_chain,1);
mean_R0 = mean(R0_chain,1);
sigma_R0 = std(R0_chain,1);
mean_dR = mean(dR_chain,1);
sigma_dR = std(dR_chain,1);
mean_sigma = sqrt(mean(s2chain));
sigma_sigma = std(sqrt(s2chain),1);

%Plotting variables
simMS2 = FitMeanAPMCMC_GetFluorFromPolPos(construct,FitMeanAPMCMC_ConstantElongationSim(v_elong,mean_ton,...
    mean_R0+mean_dR,data.xdata),v_elong,mean_MS2_basal,mean_dwelltime);
t_plot = t(nc_array(nc):nc_array(nc+1)) - t(nc_array(nc));
MS2_plot = MS2(nc_array(nc):nc_array(nc+1));

%% Save data
MCMCchain_temp(AP).ton_chain = ton_chain;
MCMCchain_temp(AP).MS2_basal_chain = MS2_basal_chain;
MCMCchain_temp(AP).dwelltime_chain = dwelltime_chain;
MCMCchain_temp(AP).R0_chain = R0_chain;
MCMCchain_temp(AP).dR_chain = dR_chain;
MCMCchain_temp(AP).s2chain = s2chain;

MCMCresults_temp(AP).mean_ton = mean_ton;
MCMCresults_temp(AP).sigma_ton = sigma_ton;
MCMCresults_temp(AP).mean_dwelltime = mean_dwelltime;
MCMCresults_temp(AP).sigma_dwelltime = sigma_dwelltime;
MCMCresults_temp(AP).mean_MS2_basal = mean_MS2_basal;
MCMCresults_temp(AP).sigma_MS2_basal = sigma_MS2_basal;
MCMCresults_temp(AP).mean_R0 = mean_R0;
MCMCresults_temp(AP).sigma_R0 = sigma_R0;
MCMCresults_temp(AP).mean_dR = mean_dR;
MCMCresults_temp(AP).sigma_dR = sigma_dR;
MCMCresults_temp(AP).mean_sigma = mean_sigma;
MCMCresults_temp(AP).sigma_sigma = sigma_sigma;

MCMCplot_temp(AP).t_plot = t_plot;
MCMCplot_temp(AP).MS2_plot = MS2_plot;
MCMCplot_temp(AP).t_interp = t_interp;
MCMCplot_temp(AP).MS2_interp = MS2_interp;
MCMCplot_temp(AP).simMS2 = simMS2;
MCMCplot_temp(AP).firstnc13time = firstnc13time;
MCMCplot_temp(AP).lastnc13time = lastnc13time;
MCMCplot_temp(AP).firstnc14time = firstnc14time;
MCMCplot_temp(AP).lastnc14time = lastnc14time;
MCMCplot_temp(AP).nc = nc_string{nc};
end

%Save nuclear cycle results into master structure array
MCMCchain{nc} = MCMCchain_temp;
MCMCresults{nc} = MCMCresults_temp;
MCMCplot{nc} = MCMCplot_temp;
end

%% Save data into .mat structure
loc = [DropboxFolder,filesep,Prefix,'\'];
filename = 'MeanFitsMCMC';
save([loc,filename,'.mat'],'MCMCresults','MCMCplot','Prefix','param_log');

filename = 'MeanFitsMCMC_RawChains';
save([loc,filename,'.mat'],'MCMCchain');

close(w);
disp(['MCMC analysis complete. Information stored in: ',loc]);

%Close the parallel pool if desired
if ~KeepPool
    poolobj = gcp('nocreate');
    delete(poolobj);
end
end