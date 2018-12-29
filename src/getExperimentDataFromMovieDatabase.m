function [Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
	Channel1, Channel2,Objective, Power,  DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3,prophase,metaphase]...
	= getExperimentDataFromMovieDatabase(Prefix, movieDatabaseFolder)

  movieDatabasePath = [movieDatabaseFolder, '/MovieDatabase.csv'];

  movieDatabase = csv2cell(movieDatabasePath, 'fromfile');

  [~, PrefixRow] = getDropboxFolderFromMovieDatabase(movieDatabasePath, Prefix, '[\\\\/-]');

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
  anaphase = [nc9, nc10, nc11, nc12, nc14];
  CF = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'CF'));
  
  % For Channel3, make this optional
  try ~isempty(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Channel3'));
      Channel3 = { getValueFromMovieDatabase(movieDatabase, PrefixRow, 'Channel3') };
  catch 
      Channel3 = {'DoesNotExist'};
  end
  
  % Making prophase and metaphase time points optional
  % assuming that nuclear cycles included are 9-14
  % for prophase
  try ~isempty(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'p9'));
      p9 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'p9'));
      p10 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'p10'));
      p11 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'p11'));
      p12 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'p12'));
      p13 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'p13'));
      p14 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'p14'));
      prophase = [p9, p10, p11, p12, p13, p14];
  catch
      prophase = [];
  end
  
  % doing the same for metaphase

  try ~isempty(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'm9'));
      m9 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'm9'));
      m10 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'm10'));
      m11 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'm11'));
      m12 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'm12'));
      m13 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'm13'));
      m14 = str2num(getValueFromMovieDatabase(movieDatabase, PrefixRow, 'm14'));
      metaphase = [m9, m10, m11, m12,m13,m14];
  catch
      metaphase = [];
  end

end