% Order particles by the earliest frame they appear at. This makes the
% tracking a lot easier! Can also track by the number of spots in a trace
function [Particles] = sortParticles(Sort, sortByLength, NChannels, Particles)
  if Sort
    direction = 'ascend';
    
    for ChN = 1:NChannels
      nParticles = length(Particles{ChN});
      sortIndex0 = zeros(nParticles, 1);
      sortIndex1 = zeros(nParticles, 1);
      sortIndex2 = zeros(nParticles, 1);
      sortIndex3 = zeros(nParticles, 1); % let's use xPos as a sub index for stable sorting
      maxFlags = max([Particles{ChN}.flagsPerFrame]);
      maxUrgentFlags = max([Particles{ChN}.urgentFlagsPerFrame]);
      for i = 1:length(Particles{ChN})
        sortIndex0(i) = maxUrgentFlags-Particles{ChN}(i).urgentFlagsPerFrame;
        sortIndex1(i) = maxFlags-Particles{ChN}(i).flagsPerFrame;        
        if sortByLength %sort by most points in particle
          sortIndex0(i) = length(Particles{ChN}(i).Frame);
          direction = 'descend';
        else %Otherwise, sort by first frame as normal
          sortIndex2(i) = Particles{ChN}(i).Frame(1);
        end
        sortIndex3(i) = Particles{ChN}(i).xPos(1);
      end
      [~, Permutations] = sortrows([sortIndex0 sortIndex1 sortIndex2 sortIndex3], direction);
      Particles{ChN} = Particles{ChN}(Permutations);
    end

  end
end