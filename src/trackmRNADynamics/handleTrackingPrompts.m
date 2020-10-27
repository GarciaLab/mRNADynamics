function retrack =  handleTrackingPrompts(liveExperiment,Spots,retrack)

  if retrack && ~liveExperiment.hasParticlesFile
    error('No Particles structure found. Re-run without "retrack" option') 
    
  % prompt user if old Particles structure exists, but retracking not specified  
  elseif ~retrack && liveExperiment.hasParticlesFile
    retrack_str = input('Particles structure detected. Do you want to perform retracking? (y/n)','s');
    if strcmpi(retrack_str,'n') % check to see if spots were added
      disp('overwriting...')
      [~, SpotFilter] = getParticles(liveExperiment);
      if ~iscell(SpotFilter)
        SpotFilter = {SpotFilter};
      end
      newSpotsFlag = 0;
      for Ch = 1:length(SpotFilter)
        newSpotsFlag = newSpotsFlag + sum(1*SpotFilter{Ch}(:)==2);
      end
      if newSpotsFlag > 0
        reset_spots_str = input(['Particles will be overwritten. Do you wish to overwrite ' ...
              num2str(newSpotsFlag) ' manually added spot(s)? (y/n)'],'s');
        if strcmpi(reset_spots_str,'y')
          disp('removing user-added spots...')        
          for Ch = 1:length(SpotFilter)
            [delFrames, delIndices] = find(SpotFilter{Ch}==2);
            for i = 1:length(delFrames)
              Spots{Ch}(delFrames(i)).Fits = Spots{Ch}(delFrames(i)).Fits([1:delIndices(i)-1 delIndices(i)+1:end]);
            end
          end
          disp('saving Spots structure...')
          save([liveExperiment.resultsFolder 'Spots.mat'],'Spots')
        end
      end
      
    elseif strcmpi(retrack_str,'y')    
      retrack = true;
    end
  end

  if retrack
    disp('retracking...');
  end