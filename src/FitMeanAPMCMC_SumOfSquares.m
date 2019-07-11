function SS = FitMeanAPMCMC_SumOfSquares(construct,v,data,x,varargin)
% Last updated: 7/5/19 by Jonathan Liu

% Log likelihood probability function for constant model elongation rate MCMC fitting.
% Currently the model infers values for the termination dwell time, the
% mean loading rate, time-dependent loading rate fluctiations, the
% transcription turn on time, the basal MS2 fluorescence level, and
% measurement noise.
% This is for use in the MCMCstat package and returns the sum-of-squares function.

% Parameters
% ----------
% construct: string.
%   String specifying the construct used.
% data: structure
%   Structure containing time data and fluorescence data.
% v: double
%   Value of elongation rate (kb/min).
% x: array
%   Array of parameters (t_on, dwelltime, fluor_basal, R0, dR)
% varargin: variable inputs
%   None so far

% Return
% ------
% SS: sum-of-squares residual

% Note that this assumes our sampling assumes Gaussian error with fixed
% variance/standard deviation for each measurement.

%% Extract fit parameters and data
ton = x(1);
dwelltime = x(2);
fluor_basal = x(3);
R0 = x(4);
dR = x(5:end);

R = R0 + dR; %Loading rate as a sum of mean rate and fluctuations

t = data.xdata;
fluorExp = data.ydata;

%% Set up Posterior function

%Calculate simulated fluorescence.
PolPos = FitMeanAPMCMC_ConstantElongationSim(v,ton,R,t);
fluorSim = FitMeanAPMCMC_GetFluorFromPolPos(construct,PolPos,v,fluor_basal,dwelltime);

%% Compute sum-of-squares
% Calculate the residuals of experimental and theoretical predictions
residuals = fluorExp - fluorSim;

% Compute the sum-of-squares function
SS = nansum(residuals.^2);
end