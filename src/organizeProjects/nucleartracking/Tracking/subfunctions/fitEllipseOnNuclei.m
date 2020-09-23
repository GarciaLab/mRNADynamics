function [ nuclei ] = fitEllipseOnNuclei( nuclei, names, varargin )
%FITELLIPSEONNUCLEI This function takes up a nuclei structure and the names
%of the images where they are present in order to fit ellipses on each of
%them. The data regarding the ellipses is stored in the nuclei structure
%under the field 'ellipse'.
% The optional parameter is a vector of with one field per image defining
% the approximate diameter of the nuclei on this frame, used to crop the
% image when fitting the ellipses. The value should thus be slightly over
% the size of the diameter to cope with imprecisions in position or 
% variation in size.

if nargin > 2
    diameters = varargin{1};
    if numel(diameters) == 1
        diameters = repmat(diameters,numel(names),1);
    end
    if numel(diameters) ~= numel(names)
        error('Invalid diameters argument. The diameters vector has to contain as many elements as there are images.')
    end
else
    diameters = 20*ones(numel(images),1);
end
        
        
for j = 1:numel(names)
    im = double(imread(names{j}));
    fprintf([num2str(j) ' '])
    for jj = 1:numel(nuclei)
        if nuclei(jj).indXY(j) > 0 && ~isnan(nuclei(jj).position(j,1))
            center = nuclei(jj).position(j,:);
            r = round(0.5*diameters(j));
            rctngl = [center(2)-r center(1)-r diameters(j) diameters(j)];
            imc = imcrop(im,rctngl);
            if numel(imc) == 0
                1;
            end
            if j == 24
                1;
            end
            warning('off')
            [f,dummy1,dummy2,dummy3,dummy4,dummy5] = fmgaussfit(1:size(imc,2),1:size(imc,1),imc);
            warning('on')
            nuclei(jj).ellipse(j).axis = f(3:4);
            nuclei(jj).ellipse(j).center = f([6,5])+rctngl(1:2);
            nuclei(jj).ellipse(j).angle = f(2);
        else
            nuclei(jj).ellipse(j).axis = [0 0];
            nuclei(jj).ellipse(j).center = [0 0];
            nuclei(jj).ellipse(j).angle = 0;
        end
    end
end
            

end

