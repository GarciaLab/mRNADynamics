function Ellipses = fillEmptyXYFrames(Ellipses)
%subfunction for maintracking
    
    
    if sum(cellfun(@(x) size(x,1),Ellipses) < 1)
        %Find the frames where we have issues
        FramesToFix=find(cellfun(@(x) size(x,1),Ellipses) < 1 );
        for i=1:length(FramesToFix)
            
            if FramesToFix(i) == 1
                
                FrameToCopy=1;
                while any(FramesToFix == FrameToCopy)
                    FrameToCopy = FrameToCopy + 1;
                end
                
            else
                FrameToCopy=FramesToFix(i)-1;
            end
            
            Ellipses{FramesToFix(i)}=Ellipses{FrameToCopy};
            disp(['fixed frame: ', num2str(FramesToFix(i)), ' using frame: ', num2str(FrameToCopy)]);
            
        end
    end

end