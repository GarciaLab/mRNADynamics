function keyInputHandler = NuclearEllipsesEventHandler(cntState)
 
    function keyInput(cc)
        if cc == 'C'
            % add ellipse
            [ConnectPositionx,ConnectPositiony] = ginput(1);
            
            cm = [round(ConnectPositionx), round(ConnectPositiony)];
            [Rows, Cols] = size(cntState.ImageHis);
            if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
                
                %Add a circle to this location with the mean radius of the
                %ellipses found in this frame
                
                %(x, y, a, b, theta, maxcontourvalue, time,
                %particle_id)
                if ~isempty(cntState.Ellipses{cntState.CurrentFrame})
                    MeanRadius=mean((cntState.Ellipses{cntState.CurrentFrame}(:,3)+cntState.Ellipses{cntState.CurrentFrame}(:,4))/2);
                elseif ~isempty(cntState.Ellipses{cntState.CurrentFrame+1})
                    MeanRadius=mean((cntState.Ellipses{cntState.CurrentFrame+1}(:,3)+cntState.Ellipses{cntState.CurrentFrame+1}(:,4))/2);
                elseif ~isempty(cntState.Ellipses{cntState.CurrentFrame-1})
                    MeanRadius=mean((cntState.Ellipses{cntState.CurrentFrame-1}(:,3)+cntState.Ellipses{cntState.CurrentFrame-1}(:,4))/2);
                end
                
                try
                    cntState.Ellipses{cntState.CurrentFrame}(end+1,:)=...
                        [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0,0];
                catch
                    cntState.Ellipses{cntState.CurrentFrame}(end+1,:)=...
                        [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0];
                end
            end
            
            cntState.nucleiModified = true;
        
        elseif cc == 'V'
            %remove ellipse
            [ConnectPositionx,ConnectPositiony] = ginput(1);
            
            cm = [round(ConnectPositionx), round(ConnectPositiony)];
            [Rows, Cols] = size(cntState.ImageHis);
            
            if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
                %Find out which ellipses we clicked on so we can delete it
                
                %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                Distances=sqrt((cntState.Ellipses{cntState.CurrentFrame}(:,1)-cm(1,1)).^2+...
                    (cntState.Ellipses{cntState.CurrentFrame}(:,2)-cm(1,2)).^2);
                [~,MinIndex]=min(Distances);
                
                cntState.Ellipses{cntState.CurrentFrame}=[cntState.Ellipses{cntState.CurrentFrame}(1:MinIndex-1,:);...
                    cntState.Ellipses{cntState.CurrentFrame}(MinIndex+1:end,:)];
            end
            
            cntState.nucleiModified = true;
        end
    end

    keyInputHandler = @keyInput;
end
