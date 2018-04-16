function [Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
	Channel1, Channel2,Objective, Power,  DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3]...
	= getExperimentDataFromMovieDatabase(Prefix, movieDatabaseFolder)

  movieDatabasePath = [movieDatabaseFolder, '/MovieDatabase.csv'];

  movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
  movieDatabaseHeaderRow = movieDatabase(1, :);

  [DropboxFolder, PrefixRow] = getDropboxFolderFromMovieDatabase(movieDatabasePath, Prefix, '[\\\\/-]');

  Date = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Date');
  ExperimentType = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'ExperimentType');
  ExperimentAxis = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'ExperimentAxis');
  CoatProtein = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'CoatProtein');
  StemLoop = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'StemLoop');
  APResolution = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'APResolution'));
  Channel1 = { getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Channel1') };
  Channel2 = { getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Channel2') };
  Objective = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Objective');
  Power = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Power');
  DataFolder = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'DataFolder');
  DropboxFolderName = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'DropboxFolder');
  Comments = getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Comments');
  nc9 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'nc9'));
  nc10 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'nc10'));
  nc11 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'nc11'));
  nc12 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'nc12'));
  nc13 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'nc13'));
  nc14 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'nc14'));
  CF = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'CF'));
  % For Channel3, make this as an optional
  if ~isempty(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Channel3'))
      Channel3 = { getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Channel3') };
  else 
      Channel3 = 'DoesNotExist';
  end

end