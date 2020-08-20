function [movieMatCh1, movieMatCh2, movieMatCh3] = makeMovieMatChannels(Prefix, varargin)

movieMatCh1 = [];
movieMatCh2 = [];
movieMatCh3 = [];


[~, ~, ~, ~, PreProcPath] = DetermineLocalFolders(Prefix);


load([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMat.Mat'], 'movieMat')

%
% %make sure the last dimension is the channel
% if size(movieMat, 5) > 3
%
%     warning('MovieMat isn''t in the right format. Permuting.');
%     movieMat = permute(movieMat, [5 4 2 3 1]); %ch t z x y --> y x z t ch
%
%     save([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMat.mat'], 'movieMat', '-v7.3', '-nocompression');
%
% end

movieMatCh1 = squeeze(movieMat(:, :, :, :, 1));
if max(movieMatCh1(:)) < 256
    movieMatCh1 = uint8(movieMatCh1);
end
save([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMatCh1.mat'], 'movieMatCh1', '-v6');

if size(movieMat, 5) > 1
    movieMatCh2 = squeeze(movieMat(:, :, :, :, 2));
    if max(movieMatCh2(:)) < 256
        movieMatCh2 = uint8(movieMatCh2);
    end
    save([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMatCh2.mat'], 'movieMatCh2', '-v6');
end

if size(movieMat, 5) == 3
    
    movieMatCh3 = squeeze(movieMat(:, :, :, :, 3));
    if max(movieMatCh3(:)) < 256
        movieMatCh3 = uint8(movieMatCh3);
    end
    
    save([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMatCh3.mat'], 'movieMatCh3', '-v6');
    
elseif size(movieMat,5) > 3
    error('movie mat has greater than 3 channels. maybe the channels are permuted?');
end

end