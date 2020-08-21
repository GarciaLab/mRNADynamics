% Order particles by the earliest frame they appear at. This makes the
% tracking a lot easier! Can also track by the number of spots in a trace
function [Particles] = sortParticles(Sort, sortByLength, NChannels, Particles)
  if Sort
    direction = 'ascend';
    
    for ChN = 1:NChannels
      nParticles = length(Particles{ChN});
      sortIndex1 = zeros(nParticles, 1);
      sortIndex2 = zeros(nParticles, 1); % let's use xPos as a sub index for stable sorting

      for i = 1:length(Particles{ChN})

        if sortByLength %sort by most points in particle
          sortIndex1(i) = length(Particles{ChN}(i).Frame);
          direction = 'descend';
        else %Otherwise, sort by first frame as normal
          sortIndex1(i) = Particles{ChN}(i).Frame(1);
        end
        sortIndex2(i) = Particles{ChN}(i).xPos(1);
      end

      [~, Permutations] = sortrows([sortIndex1 sortIndex2], direction);
      Particles{ChN} = Particles{ChN}(Permutations);
    end

  end
end