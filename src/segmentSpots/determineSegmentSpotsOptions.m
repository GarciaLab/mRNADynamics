function [displayFigures, numFrames, numShadows, intScale, nWorkers, keepPool, pool, autoThresh] = determineSegmentSpotsOptions(varargin)

  varargin = varargin{1};
  
  % Default options
  displayFigures = 0;
  numFrames = 0;
  numShadows = 2;
  intScale = 1;
  nWorkers = 8;
  keepPool = 0;
  pool = 1;
  use_integral_center = 0;
  autoThresh = 0;

  for i = 1:length(varargin)

    if strcmpi(varargin{i}, 'displayFigures')
      displayFigures = 1;
    elseif strcmpi(varargin{i}, 'Shadows')

      
      if (i + 1) > length(varargin)|| ~ isnumeric(varargin{i + 1}) || varargin{i + 1} > 2
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
    elseif strcmpi(varargin{i}, 'IntegralZ')
      use_integral_center = 1;
    elseif strcmpi(varargin{i}, 'nWorkers')
      nWorkers = varargin{i + 1};

      if nWorkers == 0
        pool = 0;
      end 
    elseif strcmpi(varargin{i}, 'autoThresh')
        autoThresh = 1;

  end 

end 
