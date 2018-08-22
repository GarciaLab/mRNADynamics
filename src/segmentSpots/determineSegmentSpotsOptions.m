function [displayFigures, trackSpots, numFrames, numShadows, customFilter, highPrecision, filterType, intScale, nWorkers, keepPool, pool] = determineSegmentSpotsOptions(varargin)
  % Default options
  displayFigures = 0;
  trackSpots = 0;
  numFrames = 0;
  numShadows = 2;
  customFilter = 0;
  highPrecision = 0;
  filterType = 'Difference_of_Gaussian';
  intScale = 1;
  nWorkers = 8;
  keepPool = 0;
  pool = 1;
  varargin = varargin{1};
  for i = 1:length(varargin)

    if strcmpi(varargin{i}, 'displayFigures')
      displayFigures = 1;
    elseif strcmpi(varargin{i}, 'TrackSpots')
      trackSpots = 1;
    elseif strcmpi(varargin{i}, 'Shadows')

      if ~ isnumeric(varargin{i + 1}) || varargin{i + 1} > 2
        error('Wrong input parameters. After ''Shadows'' you should input number of shadows(0, 1 or 2)')
      else 
        numShadows = varargin{i + 1};
      end 

    elseif strcmp(varargin{i}, 'Frames')

      if ~ isnumeric(varargin{i + 1})
        error('Wrong input parameters. After ''Frames'' you should input the number of frames')
      else 
        numFrames = varargin{i + 1};
      end 

    elseif strcmpi(varargin{i}, 'keepPool')
      pool = 1;
    elseif strcmpi(varargin{i}, 'highPrecision')
      highPrecision = 1;
    elseif strcmpi(varargin{i}, 'intScale')
      intScale = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'nWorkers')
      nWorkers = varargin{i + 1};

      if nWorkers == 0
        pool = 0;
      end 

    elseif strcmpi(varargin{i}, 'customFilter')
      customFilter = 1;

      try 
        filterType = varargin{i + 1};
      catch 
        warning('Entered filter not recognized. Defaulting to DoG')
      end 

      if iscell(varargin{i + 2})
        sigmas = varargin{i + 2};

        if strcmp(filterType, 'Difference_of_Gaussian') || ...
          strcmp(filterType, 'Structure_largest') || ...
          strcmp(filterType, 'Structure_smallest')

          if length(sigmas) ~= 2
            error('DoG and Structure filters require two sigma values e.g.{lower_sigma,higher_sigma}')
          end 

        else 

          if length(sigmas) ~= 1
            error('All filters besides DoG and Structure require only 1 sigma value')
          end 

        end 

      else 
        error('Entered sigma(s) not recognized. Make sure the sigma(s) are entered as numbers in a cell {}')
      end 

    end 

  end 

end 
