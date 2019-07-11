function MS2 = FitMeanAPMCMC_GetFluorFromPolPos(construct,PolPos,v_elong,fluor_basal,t_dwell)
%Last updated: 7/11/19 by Jonathan Liu

%Calculates an MS2 signal given a simulated position matrix of Pol
%II molecules over time on a gene. If a Pol II is past the end of a loop
%sequence, give it a unit fluorescence (or # count of Pol II). If it is
%partially done through a loop sequence, give it a fractional unit of
%fluorescence. Each Pol II pauses at the end of the gene for a fixed
%termination dwell time. Also, we assume a basal level of fluorescence.

%Inputs:
%   construct: string defining the construct used
%   PolPos: matrix giving positions of Pol II molecules over time (time x
%   space)
%   v_elong: elongation rate in kb/min
%   fluor_basal: basal level of fluorescence
%   t_dwell: termination pause/dwell time in min

% Define construct parameters. Incorporate the idea of a termination dwell
% time by extending the "effective" length of the construct by the dwell
% time multiplied by the elongation rate.
% Any additional constructs are to be defined here.
if strcmp(construct,'P2P-MS2v5-LacZ')
    L = 5.160+t_dwell*v_elong; %Length in kb of construct
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    
    MS2_loopn = 24; %Number of MS2 loops
    
    MS2 = 0;
end

%Calculate MS2 signal by making a fluorescence map
for i = 1:length(MS2_start)
    MS2_fluorval = MS2_loopn(i)/24;
    MS2map = zeros(size(PolPos));
    MS2map(and(PolPos > MS2_end(i), PolPos <L)) = MS2_fluorval; %Pol II's that are past the MS2 sequence.
    frac_ind = and(PolPos > MS2_start(i), PolPos < MS2_end(i)); %Indices for partially transcribed loops.
    MS2map(frac_ind) = (PolPos(frac_ind) - MS2_start(i)).*MS2_fluorval./(MS2_end(i) - MS2_start(i));

    MS2 = MS2 + sum(MS2map,2)';
    
    %Basal activity
    MS2(MS2<fluor_basal) = fluor_basal;
end
    
end