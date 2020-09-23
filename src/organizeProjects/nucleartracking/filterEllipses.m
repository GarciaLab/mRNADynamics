function Ellipses = filterEllipses(Ellipses, imageDims)
%clean the Ellipses structure to prevent the code erroring
%on anything weird. 

emptyFrames = find(cellfun(@isempty, Ellipses));
for f = 1:length(emptyFrames)
    Ellipses{emptyFrames(f)} = zeros(1, size(Ellipses{f}, 2));
end

%validate sizes. the ellipse masker handles
%very large objects poorly, to say the least.
for f = 1:length(Ellipses)
    Ellipses{f} = filterEllipseFrame(Ellipses{f}, imageDims);
end