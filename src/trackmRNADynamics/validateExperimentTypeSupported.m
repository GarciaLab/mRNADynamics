%subfunction for TrackmRNADynamics but could be repurposed too
function validateExperimentTypeSupported(ExperimentType)

  if ~(strcmpi(ExperimentType, '1spot') || strcmpi(ExperimentType, '2spot') || ...
      strcmpi(ExperimentType, 'inputoutput') || ...
      strcmpi(ExperimentType, '2spot2color') || ...
      strcmpi(ExperimentType, 'lattice'))

    error('Experiment type in MovieDatabase not recognized')

  end

end