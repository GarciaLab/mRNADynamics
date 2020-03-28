function [anaphaseFrames, anaphaseFile] = retrieveAnaphaseFrames(Prefix)

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

anaphaseFile = [DropboxFolder,filesep,Prefix,filesep, 'anaphaseFrames.mat'];

if exist(anaphaseFile, 'file')
    
	load(anaphaseFile, 'anaphaseFrames')

else
    
	try

	[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~,Channel1,Channel2,~,...
    Channel3, ~, movieDatabaseFolder, movieDatabase]...
    = readMovieDatabase(Prefix);
    
    [   ~, ~, ~, ~, ~, ~, ~,~,~, ~,  ~, ~, ~,...
            nc9, nc10, nc11, nc12, nc13, nc14,]...
            = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);
			
        anaphaseFrames = [nc9; nc10; nc11; nc12; nc13; nc14];
    
	catch anaphaseFrames = [0; 0; 0; 0; 0; 1]; end
    
end

if numel(anaphaseFrames) < 6
    anaphaseFrames = vertcat(anaphaseFrames, nan(6-numel(anaphaseFrames), 1));
end

end