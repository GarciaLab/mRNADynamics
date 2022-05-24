function filteredParticles = filterParticlesForClusterAnalysis(...
                                          Particles, ncStart, minTraceLen)
%
% DESCRIPTION
% This function filters MS2 Particles structure based on user-input quality
% control metrics.
%
%
% ARGUMENTS
% Particles: Particles structure for the experiment you want to filter
%
% ncStart: frame corresponding to start time of nuclear cycle (nc) user
%          wants to analyze
%
% minTraceLength: minimum number of frames a particle trace must have
%
%
% OPTIONS
% NA
%
% OUTPUT
% filteredParticles: Particles structure with all particles that don't meet
%                    QC criteria removed
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2022/05/09
% Last Updated:


%% Enforce QC criteria
% Get trace length & trace start time info
traceLengths = arrayfun(@(x) size(Particles(x).Frame, 2), 1:numel(Particles));
traceStartFrames = arrayfun(@(x) Particles(x).Frame(1), 1:numel(Particles));

% Keep only traces that:
% (1) meet minimum MS2 particle trace length 
% AND
% (2) start during the specified nuclear cycle (nc)
qcFilter = (traceLengths >= minTraceLen) & (traceStartFrames >= ncStart);
keepParticles = Particles(qcFilter);

%% Check that we found at least one particle that meets QC criteria
filteredParticles = keepParticles;

if isempty(filteredParticles)
    error('No particles met your filtering criteria. Please adjust inputs and try again.')
end