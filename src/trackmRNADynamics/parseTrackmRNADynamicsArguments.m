function [Prefix, app, retrack, optionalResults, displayFigures] = parseTrackmRNADynamicsArguments(DefaultDropboxFolder, varargin)
  app = {};
  retrack = 1;
  optionalResults = '';
  displayFigures = 0;
  
  %Look at the input parameter and use defaults if missing
  if isempty(varargin)
    Prefix = promptUserForFolder(DefaultDropboxFolder);
  elseif ~ischar(varargin{1})
    Prefix = promptUserForFolder(DefaultDropboxFolder);
  else
    Prefix = varargin{1};

    for i = 2:length(varargin)

      if strcmpi(varargin{i}, 'app')
        app{1} = varargin{i + 1};
        app{2} = varargin{i + 2};
      elseif strcmpi(varargin{i}, 'noRetracking')
          retrack = 0;
      elseif strcmpi(varargin{i}, 'optionalResults')
          optionalResults = varargin{i+1};
      elseif strcmpi(varargin{i}, 'displayFigures')
          displayFigures = 1;
      end
      

    end

  end

end

function Prefix = promptUserForFolder(DefaultDropboxFolder)
  FolderTemp = uigetdir(DefaultDropboxFolder, 'Select folder with data to analyze');
  Dashes = strfind(FolderTemp, filesep);
  Prefix = FolderTemp((Dashes(end) + 1):end);
end