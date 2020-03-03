function im = sliceMovieMat(movieMat, varargin)

%y x z t ch

if ~isempty(varargin)
    
    if length(varargin) == 1
        ch = varargin{1};
        im = squeeze(movieMat(:, :, :, :, ch));
    elseif length(varargin) == 2
        ch = varargin{1}; f = varargin{2};
        im = squeeze(movieMat(:, :, :, f, ch));
    elseif length(varargin) == 3
        ch = varargin{1}; f = varargin{2}; z = varargin{1};
        im = squeeze(movieMat(:, :, z, f, ch));
    end
    
else
    
    im = movieMat;
    
end