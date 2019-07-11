function x = FitMeanAPMCMC_ConstantElongationSim(v,ton,R,t)
%Simulation of Pol II elongation using constant elongation type of model. Each
%individual Pol II is simulated and its position is kept track of. The Pol
%II molecule has a drift term. At each timestep,
%R_i new Pol II molecules are initiated, where i is the timestep we are on.
%If R_i is non-integer, then it is rounded down to the nearest integer.
%All existing Pol II molecules are progressed forward, using the deterministic
%elongation rate v.

%More specifically, each timestep progresses a Pol II by the deterministic
%drift length (v * dt). The simulation does not
%include the finite gene length, and instead keeps track of the positions
%of all the Pol II molecules, to be used in any further analysis in a
%separate script. Also note that to be more easily compared with data, the
%simulation uses a time series t that may not have fixed timesteps in
%between.

%For the definitions below, n is the number of Pol II molecules and m is
%the number of timepoints.

%The inputs are:
%   v: elongation rate (the drift term) (scalar, units of bp/time)
%   ton: time on of initation (no loading allowed before this time)
%   R: loading rate of Pol II at each timepoint (1x(m-1) vector, where m is the number of timepoints)
%   t: timepoints of simulation (1xm vector, units of time)

%The output is:
%   x: position of each Pol II molecule as a function of time (mxn matrix)

%% Initialize variables

%If there are any negative rates, make them zero.
R(R<0) = 0;

%Number of timepoints and RNAP molecules
m = size(t,2);

%Timestep definition
dt = nan(1,m-1);
for i = 1:(m-1)
    dt(i) = t(i+1) - t(i);
end

n = floor(sum(R.*dt,2));

%Start all the Pol II molecules at the origin
x = zeros(m,n);

%Counter to keep track of how many Pol II molecules have been loaded
counter = 0;

%% Time evolution of RNAP molecules
for i = 1:(m-1)
    if t(i) < ton
        continue
    else
        counter = counter + R(i)*dt(i); %Update Pol II counter
        k = 1:floor(counter); %Indices of loaded Pol II's.

        %Only evolve Pol II molecules that have been loaded
        x(i+1,k) = x(i,k) + v*dt(i);
        x(x(i+1,k)<0)=0;
    end
end

end

