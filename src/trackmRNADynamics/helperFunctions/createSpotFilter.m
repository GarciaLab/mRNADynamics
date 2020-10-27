function SpotFilter = createSpotFilter(Spots)

  NCh = length(Spots);
  MaxSpots = cell(NCh);
  SpotFilter = cell(1,NCh);
  for Channel = 1:NCh
      %Determine the maximum number of spots in a given frame for the
      %whole movie
      MaxSpots{Channel} = 0;

      for i = 1:length(Spots{Channel})
          MaxSpots{Channel} = max([MaxSpots{Channel}, length(Spots{Channel}(i).Fits)]);
      end
   
      SpotFilter{Channel} = nan(length(Spots{Channel}), MaxSpots{Channel});
      % Populate the filter
      for i = 1:length(Spots{Channel})

          for j = 1:length(Spots{Channel}(i).Fits)

              % Initializes the filter as all 1s, since the Threshold has been removed
              SpotFilter{Channel}(i, j) = 1;

          end

      end

  end
    
end
