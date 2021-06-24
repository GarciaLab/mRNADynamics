% Order particles by the earliest frame they appear at. This makes the
% tracking a lot easier! Can also track by the number of spots in a trace
function [Particles] = sortParticles(Sort, sortByLength, NChannels, Particles)
  if Sort
    direction = 'ascend';
    
    for ChN = 1:NChannels
      nParticles = length(Particles{ChN});
      sortIndex = zeros(1, nParticles);

      for i = 1:length(Particles{ChN})

        if sortByLength %sort by most points in particle
          sortIndex(i) = length(Particles{ChN}(i).Frame);
          direction = 'descend';
        else %Otherwise, sort by first frame as normal
          sortIndex(i) = Particles{ChN}(i).Frame(1);
        end

      end

      [~, Permutations] = sort(sortIndex, direction);
      Particles{ChN} = Particles{ChN}(Permutations);
    end

  end
end