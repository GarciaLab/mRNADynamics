function [MeanOnRatio,SEOnRatio]=NuclearOnRatio(Data,MinParticles,MinEmbryos,PixelsPerLine,LinesPerFrame)

%What's the fraction of active schnitzs in a frame 10 minutes into the nuclear cycle?
OnRatio=nan(length(Data),length(Data(1).APbinID),2);
for nc=13:14
    for i=1:length(Data)

        if eval(['Data(i).nc',num2str(nc)])>0
        
            %Find the frame 10 minutes into the nuclear cycle
            [Dummy,MiddleFrame]=min(((Data(i).ElapsedTime-...
                Data(i).ElapsedTime(eval(['Data(i).nc',num2str(nc)]))-10).^2));

            %Make a list of which schnitz each particle is assigned to
            ParticleNuclei=[Data(i).CompiledParticles.Nucleus];
            if length(ParticleNuclei)~=length(Data(i).CompiledParticles)
                error('Some particles are missing a nucleus')
            end


            %Make a list of all the schnitz present in this frame
            SchnitzInFrame=[];
            for k=1:length(Data(i).schnitzcells)
                if sum(Data(i).schnitzcells(k).frames==MiddleFrame)
                    SchnitzInFrame=[SchnitzInFrame,k];
                end
            end

            %Only consider nuclei that are within a certain margin of the edge of the
            %window. We'll put this information in EllipsesFilter

            EdgeWidth=10;

            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Radii=(Data(i).Ellipses{MiddleFrame}(:,3)+...
                Data(i).Ellipses{MiddleFrame}(:,4))/2;
            EllipsesFilter=((Data(i).Ellipses{MiddleFrame}(:,1)-Radii-EdgeWidth)>0)&...
                ((Data(i).Ellipses{MiddleFrame}(:,2)-Radii-EdgeWidth)>0)&...
                ((Data(i).Ellipses{MiddleFrame}(:,1)+Radii+EdgeWidth)<PixelsPerLine)&...
                ((Data(i).Ellipses{MiddleFrame}(:,2)+Radii+EdgeWidth)<LinesPerFrame);

            for j=1:length(Data(i).APbinID)

                %Which Ellipses are in this AP position?
                if j<length(Data(i).APbinID)
                    EllipsesAPFIlter=...
                        (Data(i).EllipsePos{MiddleFrame}>=Data(i).APbinID(j))&...
                        (Data(i).EllipsePos{MiddleFrame}<Data(i).APbinID(j+1));
                else
                    EllipsesAPFIlter=...
                        (Data(i).EllipsePos{MiddleFrame}>=Data(i).APbinID(j));
                end

                %Now determine which Ellipses we want to find the corresponding schnitz for
                EllipsesToCheck=find(EllipsesFilter&EllipsesAPFIlter')';

                %We will only go ahead if we have at least a MinParticles
                %number of ellipses to check

                if length(EllipsesToCheck)>=MinParticles

                    SchnitzToCheck=[];
                    for m=EllipsesToCheck
                        for k=SchnitzInFrame
                           if Data(i).schnitzcells(k).cellno...
                                    (find(Data(i).schnitzcells(k).frames==MiddleFrame))...
                                    ==m
                                SchnitzToCheck=[SchnitzToCheck,k];
                            end
                        end
                    end

                    %Do we have a discrepancy between Schnitz and Ellipses?
                    if ~(length(EllipsesToCheck)==length(SchnitzToCheck))
                        error('The number of ellipses and schnitz found in this frame and AP bin does not match!');
                    end

                    %Finally, figure out which one of these is active
                    OnRatio(i,j,nc-12)=...
                        sum(ismember(ParticleNuclei,SchnitzToCheck))/length(SchnitzToCheck);
                end
            end
        end    
    end
end


%Average all the embryos
MeanOnRatio=nan(length(Data(1).APbinID),2);
SEOnRatio=nan(length(Data(1).APbinID),2);
for nc=13:14
    for j=1:length(Data(1).APbinID)
        Temp=OnRatio(:,j,nc-12);
        Temp=Temp(~isnan(Temp));
        if length(Temp)>=MinEmbryos
            MeanOnRatio(j,nc-12)=mean(Temp);
            SEOnRatio(j,nc-12)=std(Temp)/sqrt(length(Temp));
        end
    end
end

