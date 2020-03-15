function [anaphaseFrames, anaphaseFile] = retrieveAnaphaseFrames(Prefix)

   [~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);
    
    [~,~,DropboxFolder,~, PreProcPath,...
        ~, ~, ~,Channel1,Channel2,~,...
        Channel3, ~, movieDatabaseFolder, movieDatabase]...
        = readMovieDatabase(Prefix);

    anaphaseFile = [DropboxFolder,filesep,Prefix,filesep, 'anaphaseFrames.mat'];
    if exist(anaphaseFile, 'file')
        load(anaphaseFile, 'anaphaseFrames')
    else
        try
        [   ~, ~, ~, ~, ~, ~, ~,~,~, ~,  ~, ~, ~,...
            nc9, nc10, nc11, nc12, nc13, nc14,]...
            = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);
        anaphaseFrames = [nc9; nc10; nc11; nc12; nc13; nc14];
        catch anaphaseFrames = [0; 0; 0; 0; 0; 0]; end
            
    end

end