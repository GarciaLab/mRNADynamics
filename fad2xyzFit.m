function [xs ys zs] = fad2xyzFit(ch, fad,addMargin)
% returns the list of xyz coordinates of all spots in the fits array of channel ch
% usage: 
%   [XS YS ZS] = FAD2XYZ(ch, fad) 
%   XYZ = FAD2XYZ(ch, fad) 
%   ... = FAD2XYZ(ch, fad,'addMargin') 
% if you specify 'addMargin', returns the coordinates in the coordinate
% frame of original images, i.e. of total size 2048 by 2048 (to use when
% overlaying on DAPI etc.) 

fits = fad.channels(ch).fits;

if ~(exist('addMargin','var') && ~isempty(addMargin))
    addMargin = 'none';
end

if isfield(fits,'amp')
    % compactFAD format
    try
        %This only seems to work with the compact FAD format
        xs=double(fits.x_fit);
        ys=double(fits.y_fit);
        zs=double(fits.z);
    catch
        %I think this is for the non-compact / full format
        xs=extractfield(fits,'x_fit');
        ys=extractfield(fits,'y_fit');
        zs=extractfield(fits,'z');
    end
else
    xs=extractfield(fits,'x');
    ys=extractfield(fits,'y');
    zs=extractfield(fits,'z');
end

if ~isempty(xs)
    if isfield(fad,'originalSize')
        originalSizeX = fad.originalSize(2);
        originalSizeY = fad.originalSize(1);
    else
        originalSizeX = 2^ceil(log2(max(xs)));
        originalSizeY = 2^ceil(log2(max(ys)));
    end

    if strcmpi(addMargin,'addMargin')
        shiftX=(originalSizeX-fad.stackSize(2))/2;
        shiftY=(originalSizeY-fad.stackSize(1))/2;
    else
        shiftX=0;
        shiftY=0;    
    end
    xs = xs + shiftX;
    ys = ys + shiftY;
end

if nargout==1
    xs = [xs(:) ys(:) zs(:)]';
end

end