function bin = wshed(b, varargin)

%wrapper function for watersheeding. main input parameter would be
%distThresh;

displayFigures = false;
distThresh = 3;
invert = false;

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end
df = displayFigures;

%invert the image polarity if necessary. not sure if this is the best
%metric. size or euler number could also work
stats =  regionprops(logical(b), 'EulerNumber');
statsInvert =  regionprops(logical(~b), 'EulerNumber');
if min([statsInvert.EulerNumber]) < min([stats.EulerNumber])
    b=~b;
end

bw = ~b;
D = bwdist(~bw);
D = -D;

%qc
mask = imextendedmin(D,distThresh);
D2 = imimposemin(D,mask);
%

L = watershed(D2);
L(~bw) = 0;

if df
    figure();
    rgb = label2rgb(L,'jet',[.5 .5 .5]);
    imshow(rgb, []); colorbar;
end

bin = ~~L;

end