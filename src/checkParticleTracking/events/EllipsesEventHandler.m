function keyInputHandler = EllipsesEventHandler(cptState)
 
    function keyInput(cc)
        if cc == 'C'
            % add ellipse
            [ConnectPositionx,ConnectPositiony] = ginput(1);
            
            cm = [round(ConnectPositionx), round(ConnectPositiony)];
            [Rows, Cols] = size(cptState.ImageHis);
            if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
                
                %Add a circle to this location with the mean radius of the
                %ellipses found in this frame
                
                %(x, y, a, b, theta, maxcontourvalue, time,
                %particle_id)
                if ~isempty(cptState.Ellipses{cptState.CurrentFrame})
                    MeanRadius=mean((cptState.Ellipses{cptState.CurrentFrame}(:,3)+cptState.Ellipses{cptState.CurrentFrame}(:,4))/2);
                elseif ~isempty(cptState.Ellipses{cptState.CurrentFrame+1})
                    MeanRadius=mean((cptState.Ellipses{cptState.CurrentFrame+1}(:,3)+cptState.Ellipses{cptState.CurrentFrame+1}(:,4))/2);
                elseif ~isempty(cptState.Ellipses{cptState.CurrentFrame-1})
                    MeanRadius=mean((cptState.Ellipses{cptState.CurrentFrame-1}(:,3)+cptState.Ellipses{cptState.CurrentFrame-1}(:,4))/2);
                end
                
                try
                    cptState.Ellipses{cptState.CurrentFrame}(end+1,:)=...
                        [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0,0];
                catch
                    cptState.Ellipses{cptState.CurrentFrame}(end+1,:)=...
                        [cm(1,1),cm(1,2),MeanRadius,MeanRadius,0,0,0,0];
                end
            end
            
            cptState.nucleiModified = true;
        
        elseif cc == 'V'
            %remove ellipse
            [ConnectPositionx,ConnectPositiony] = ginput(1);
            
            cm = [round(ConnectPositionx), round(ConnectPositiony)];
            [Rows, Cols] = size(cptState.ImageHis);
            
            if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Cols)
                %Find out which ellipses we clicked on so we can delete it
                
                %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                Distances=sqrt((cptState.Ellipses{cptState.CurrentFrame}(:,1)-cm(1,1)).^2+...
                    (cptState.Ellipses{cptState.CurrentFrame}(:,2)-cm(1,2)).^2);
                [~,MinIndex]=min(Distances);
                
                cptState.Ellipses{cptState.CurrentFrame}=[cptState.Ellipses{cptState.CurrentFrame}(1:MinIndex-1,:);...
                    cptState.Ellipses{cptState.CurrentFrame}(MinIndex+1:end,:)];
            end
            
            cptState.nucleiModified = true;
        end
    end

    keyInputHandler = @keyInput;
end
