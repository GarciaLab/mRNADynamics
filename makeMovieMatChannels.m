function [movieMatCh1, movieMatCh2, movieMatCh3] = makeMovieMatChannels(Prefix, varargin)

[~, ~, ~, ~, PreProcPath] = DetermineLocalFolders(Prefix);


load([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMat.Mat'], 'movieMat')

%make sure the last dimension is the channel
if size(movieMat, 5) <= 3
    
    movieMatCh1 = squeeze(movieMat(:, :, :, :, 1));
    movieMatCh2 = squeeze(movieMat(:, :, :, :, 2));
    movieMatCh3 = squeeze(movieMat(:, :, :, :, 3));

    save([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMatCh1.mat'], 'movieMatCh1', '-v7.3', '-nocompression');
    save([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMatCh2.mat'], 'movieMatCh2', '-v7.3', '-nocompression');
    save([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMatCh3.mat'], 'movieMatCh3', '-v7.3', '-nocompression');

else
    
    warning('MovieMat isn''t in the right format. I need help.');
    keyboard;
    
end


end