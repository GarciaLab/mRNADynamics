function [Prefix, app] = parseTrackmRNADynamicsArguments(DefaultDropboxFolder, varargin)
  app = {};
  
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
      end

    end

  end

end

function Prefix = promptUserForFolder(DefaultDropboxFolder)
  FolderTemp = uigetdir(DefaultDropboxFolder, 'Select folder with data to analyze');
  Dashes = strfind(FolderTemp, filesep);
  Prefix = FolderTemp((Dashes(end) + 1):end);
end