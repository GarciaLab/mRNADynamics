function [displayFigures, numFrames, initialFrame, customFilter, highPrecision, filterType, keepPool, sigmas, nWorkers, app, kernelSize, weka, justTifs, ignoreMemoryCheck, classifierFolder, classifierPathCh1] = determineFilterMovieOptions(varargin)

  varargin = varargin{1};

  % Default options
  displayFigures = 0;
  numFrames = 0;
  customFilter = 0;
  highPrecision = 0;
  keepPool = 0;
  filterType = 'Difference_of_Gaussian';
  nWorkers = 8;
  sigmas = {};
  app = {};
  kernelSize = [];

  weka = false;
  justTifs = false;
  ignoreMemoryCheck = false;
  initialFrame = 1;
  %Added new argument to specify a preferred classifier name and enable automatic testing
  classifierPathCh1 = [];
  classifierFolder = [];

  for i = 1:length(varargin)

    if strcmpi(varargin{i}, 'displayFigures')
      displayFigures = 1;

    elseif strcmp(varargin{i}, 'Frames') || strcmpi(varargin{i}, 'LastFrame')

      if ~ isnumeric(varargin{i + 1})
        error('Wrong input parameters. After ''Frames'' you should input the number of frames')
      else
        numFrames = varargin{i + 1};
      end

    elseif strcmpi(varargin{i}, 'InitialFrame')

      if ~ isnumeric( varargin{i + 1} || (varargin{i + 1} < 1) )
        error('Wrong input parameter for initial frame.')
      else
        initialFrame = varargin{i + 1};
      end

    elseif strcmpi(varargin{i}, 'highPrecision')
      highPrecision = 1;
    elseif strcmpi(varargin{i}, 'keepPool')
      keepPool = 1;
    elseif strcmpi(varargin{i}, 'app')
      app{1} = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'nWorkers')

      nWorkers = varargin{i + 1};

    elseif strcmpi(varargin{i}, 'kernelSize')

      kernelSize = varargin{i + 1};

    elseif strcmpi(varargin{i}, 'ignoreMemoryCheck')
      ignoreMemoryCheck = true;

    elseif strcmpi(varargin{i}, 'tifs')
      justTifs = true;

    elseif strcmpi(varargin{i}, 'weka')
      weka = true;

    elseif isobject(varargin{i}) && isa(varargin{i}, 'ClassifierForTest')
      ClassifierForTest = varargin{i};
      classifierFolder = ClassifierForTest.classifierFolder;
      classifierPathCh1 = ClassifierForTest.classifierPathCh1;      

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
