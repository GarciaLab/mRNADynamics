function [Prefix, Threshold1, Threshold2, app, bypassUserPrompt] = parseTrackmRNADynamicsArguments(DefaultDropboxFolder, varargin)
  app = {};
  bypassUserPrompt = false;
  
  %Look at the input parameter and use defaults if missing
  if isempty(varargin)
    Prefix = promptUserForFolder(DefaultDropboxFolder);

    % This is the first threshold we apply. We'll use it for
    % the initial search and then switch to Threshold2.
    Threshold1 = 100;
    Threshold2 = 30;
  elseif ~ischar(varargin{1})
    Prefix = promptUserForFolder(DefaultDropboxFolder);

    %Thresholds
    Threshold1 = varargin{1};
    Threshold2 = varargin{2};
  else
    Prefix = varargin{1};

    %Thresholds

    %2 spot channel input should be like TrackmRNADynamics(Prefix,
    %[thresh1, thresh2], [thresh3, thresh4])
    Threshold1 = varargin{2};
    Threshold2 = varargin{3};

    for i = 4:length(varargin)

      if strcmpi(varargin{i}, 'app')
        app{1} = varargin{i + 1};
        app{2} = varargin{i + 2};
      elseif strcmpi(varargin{i}, 'bypassUserPrompt')
        bypassUserPrompt = true;
      end

    end

  end

end

function Prefix = promptUserForFolder(DefaultDropboxFolder)
  FolderTemp = uigetdir(DefaultDropboxFolder, 'Select folder with data to analyze');
  Dashes = strfind(FolderTemp, filesep);
  Prefix = FolderTemp((Dashes(end) + 1):end);
end