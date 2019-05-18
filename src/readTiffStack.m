function imStack = readTiffStack(imPath, varargin)
        
    info = imfinfo(imPath);
    zInitial = 1; zFinal = length(info);

    for i = 1:length(varargin)
        if contains(varargin{i}, 'pad')
            zInitial = 2;
            zFinal = length(info) - 2;
        end
    end

    imStack = zeros(info(1).Height, info(1).Width, zFinal);

    for k = zInitial:zFinal
        imStack(:,:,k) = imread(imPath, k, 'Info', info);
    end
        
end