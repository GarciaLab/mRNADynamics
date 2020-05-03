function xy = fillEmptyXYFrames(xy)
%subfunction for maintracking
   
    if sum(cellfun(@(x) size(x,1),xy) < 1)
        %Find the frames where we have issues
        FramesToFix=find(cellfun(@(x) size(x,1),xy) < 1 );
        for i=1:length(FramesToFix)
            if FramesToFix(i)==1
                FrameToCopy=1;
                while sum(FramesToFix==NextFrameToCopy)
                    FrameToCopy=FrameToCopy+1;
                end
            else
                FrameToCopy=FramesToFix(i)-1;
            end
            xy{FramesToFix(i)}=xy{FrameToCopy};
        end
    end

end