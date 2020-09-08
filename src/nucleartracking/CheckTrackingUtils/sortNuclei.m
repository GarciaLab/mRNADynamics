% sortNuclei.m
% author: Gabriella Martini
% date created: 9/7/20
% date last modified: 9/7/20
% Does not support multiple protein channels
% Order particles by the earliest frame they appear at. This makes the
% tracking a lot easier! Can also track by the number of spots in a trace
function [schnitzcells] = sortNuclei(Sort, sortByLength, schnitzcells)
  if Sort
    direction = 'ascend';
    

  nNuclei = length(schnitzcells);
  sortIndex = zeros(1, nNuclei);

  for i = 1:nNuclei

    if sortByLength %sort by most points in particle
      sortIndex(i) = length(schnitzcells(i).frames);
      direction = 'descend';
    else %Otherwise, sort by first frame as normal
      sortIndex(i) = schnitzcells(i).frames(1);
    end

  end

  [~, Permutations] = sort(sortIndex, direction);
  schnitzcells = schnitzcells(Permutations);

  end
end