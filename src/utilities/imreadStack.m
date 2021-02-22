function imStack = imreadStack(imPath, varargin)

%% arguments
for i = 1:length(varargin)
    if contains(varargin{i}, 'pad')
        info = imfinfo(imPath);
        
        zInitial = 2;
        zFinal = length(info) - 2;
    end
end

%%
if contains(imPath, '.tif')
    if ~exist('info', 'var')
        info = imfinfo(imPath);
    end
    
    if ~exist('zInitial', 'var')
        zInitial = 1;
    end
    if ~exist('zFinal', 'var')
        zFinal = length(info);
    end
    
    imStack = zeros(info(1).Height, info(1).Width, zFinal);
    
    for k = zInitial:zFinal
        imStack(:,:,k) = imread(imPath, k, 'Info', info);
    end
    
elseif contains(imPath, '.mat')
    imStack = load(imPath);
    varName = fieldnames(imStack);
    imStack = imStack.(varName{1});
end

assert(logical(exist('imStack', 'var')));

end