function [displayFigures, numFrames, numShadows, intScale, nWorkers, keepPool, pool, autoThresh, initialFrame, useIntegralCenter] = determineSegmentSpotsOptions(varargin)

  varargin = varargin{1};

  % Default options
  displayFigures = 0;
  numFrames = 0;
  numShadows = 2;
  intScale = 1;
  nWorkers = 8;
  keepPool = 0;
  pool = 1;
  autoThresh = 0;
  % Default is 1
  useIntegralCenter = 1;
  initialFrame = 1;

  for i = 1:length(varargin)

    if strcmpi(varargin{i}, 'displayFigures')
      displayFigures = 1;

    elseif strcmpi(varargin{i}, 'Shadows')

      if (i + 1) > length(varargin) || ~ isnumeric(varargin{i + 1}) || varargin{i + 1} > 2
        error('Wrong input parameters. After ''Shadows'' you should input number of shadows(0, 1 or 2)')
      else
        numShadows = varargin{i + 1};
      end

    elseif strcmp(varargin{i}, 'Frames') || strcmpi(varargin{i}, 'LastFrame')

      if ~ isnumeric(varargin{i + 1})
        error('Wrong input parameters. After ''Frames'' you should input the number of frames')
      else
        numFrames = varargin{i + 1};
      end

    elseif strcmpi(varargin{i}, 'InitialFrame')

      if ~ isnumeric(varargin{i + 1}) || varargin{i + 1} < 1
        error('Wrong input parameter for initial frame.')
      else
        initialFrame = varargin{i + 1};
      end

    elseif strcmpi(varargin{i}, 'keepPool')
      pool = 1;
      keepPool = 1;
    elseif strcmpi(varargin{i}, 'intScale')
      intScale = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'noIntegralZ')
      useIntegralCenter = 0;
    elseif strcmpi(varargin{i}, 'nWorkers')
      nWorkers = varargin{i + 1};

      if nWorkers == 0
        pool = 0;
      end

    elseif strcmpi(varargin{i}, 'autoThresh')
      autoThresh = 1;

    elseif strcmpi(varargin{i}, 'tifs')
      error('Tifs generation is no longer supported from segmentSpotsML, try filterMovie(Prefix, ''Tifs'') instead.');
    else

      if ~ isnumeric(varargin{i})
        error('Input parameters not recognized. Check spelling and case.')
      end

    end

  end
