function [schnitzcells,Particles]=CheckSchnitzConsistency(schnitzcells,Particles,Ellipses,varargin)

%Makes sure there are no weird inconsistencies in the schnitzcells
%structure

%Optional parameters:
%OnlyCellNo: Only check for missing cellno part of schnitzcells

OnlyCellNo=0;
if ~isempty(varargin)
    if strcmp(varargin{1},'OnlyCellNo')
        OnlyCellNo=1;
    end
end




display('Looking for problems in schnitzcells and Particles structures')


if ~OnlyCellNo

    %See which particles are assigned to nuclei
    for i=1:length(Particles)
        if ~isempty(Particles(i).Nucleus)
            AssignedNuclei(i)=Particles(i).Nucleus;
        else
            AssignedNuclei(i)=0;
        end
    end



    %Check if there are schnitzcells with repeated frames. Where do these come
    %from?
    for i=1:length(schnitzcells)
        while min(diff(schnitzcells(i).frames))==0
            %There are repeated frames!

            %Split the schnitz at the last repeated frame. We keep going in
            %case there are multiple repeated frames.
            schnitzcells=SplitSchnitzRepeatedFrames(schnitzcells,i);

            %If there is an associated particle, change its approved flag
            %to 0
            Particles(find(AssignedNuclei==i)).Approved=0;
        end
    end


    %Check if there are schnitz that are self-referential. For example, does a
    %schnitz think it's its own parent or daughter nucleus?
    for i=1:length(schnitzcells)
        if schnitzcells(i).P==i
            schnitzcells(i).P=0;
        end
        if schnitzcells(i).E==i
            schnitzcells(i).E=0;
        end 
        if schnitzcells(i).D==i
            schnitzcells(i).D=0;
        end
    end



    %Check if there are particles that have an approved field that is smaller
    %than the number of frames
    for i=1:length(Particles)
        if length(Particles(i).Frame)~=length(Particles(i).FrameApproved)
            Particles(i).FrameApproved=logical(ones(size(Particles(i).Frame)));
        end
    end


    %Check if there are particles with multiple frames. If so we'll split them
    %and reset their nuclear information

    %First sort within the particles according to their frames
    for i=1:length(Particles)
        [Particles(i).Frame,Permutations]=sort(Particles(i).Frame);
        Particles(i).Index=Particles(i).Index(Permutations);
        if isfield(Particles,'FrameApproved')
            Particles(i).FrameApproved=Particles(i).FrameApproved(Permutations);
        end
    end

    %Now look for repeated instances like we did for schnitzcells above
    for i=1:length(Particles)
        while min(diff(Particles(i).Frame))==0
            %There are repeated frames!

            %Find the last location where the particle skips a frame and create a
            %filter.
            FrameFilter=logical(zeros(size(Particles(i).Frame)));
            FrameFilter(max(find(diff(Particles(i).Frame)==0))+1:end)=1;

            NParticles=length(Particles);

            Particles(NParticles+1).Frame=Particles(i).Frame(FrameFilter);
            Particles(NParticles+1).Index=Particles(i).Index(FrameFilter);
            Particles(NParticles+1).Nucleus=[];
            if isfield(Particles,'FrameApproved')
                Particles(NParticles+1).FrameApproved=Particles(i).FrameApproved(FrameFilter);
            end
            Particles(NParticles+1).Approved=0;

            Particles(i).Frame=Particles(i).Frame(~FrameFilter);
            Particles(i).Index=Particles(i).Index(~FrameFilter);
            if isfield(Particles,'FrameApproved')
                Particles(i).FrameApproved=Particles(i).FrameApproved(~FrameFilter);
            end
            Particles(i).Approved=0;
        end
    end
end



%Sometimes cellno gets screwed up. In those cases we need to regenerate it
%from the ellipses
for k=1:length(schnitzcells)

    %Check this schnitz for consistency with cellno.
    %Otherwise fix it. I obtained the fixing code from
    %TrackmRNADynamicsV2.m
    if length(schnitzcells(k).frames)~=length(schnitzcells(k).cellno)
       %If there number of frames is different from the number of
       %cellno then use the cenx and ceny to find the cellno in
       %Ellipses an repopulate this schnitz
       if (length(schnitzcells(k).frames)==length(schnitzcells(k).cenx))&...
               (length(schnitzcells(k).frames)==length(schnitzcells(k).ceny))
           for m=1:length(schnitzcells(k).frames)
                %The information in Ellipses is
                %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                MaxDistance=2;  %Maximum pixel distance to identify an
                                %ellipse with a schnitz
                Distances=sqrt((Ellipses{schnitzcells(k).frames(m)-1}(:,1)-...
                   schnitzcells(k).cenx(m)).^2+...
                   (Ellipses{schnitzcells(k).frames(m)-1}(:,2)-...
                   schnitzcells(k).ceny(m)).^2);
                [MinValue,MinIndex]=min(Distances);

               %Make sure no other schnitz is associated to this
               %ellipse
               EllipseFoundElsewhere=0;
               for n=[1:k-1,k+1:length(schnitzcells)]
                   %Only consider it if the schnitzcell is also valid!
                   if (length(schnitzcells(n).frames)==length(schnitzcells(n).cellno))
                       if sum(schnitzcells(n).frames==schnitzcells(k).frames(m))
                          IndexToCheck=find(schnitzcells(n).frames==schnitzcells(k).frames(m));

                          %The schnitz I'm comparing to
                          %might also be screwed up.
                          %I'd have to compare its cenx
                          %and ceny to be sure
                          try
                              if schnitzcells(k).cellno(IndexToCheck)==MinIndex
                                    error('duplicated schnitz?')
                                    DistancesK=sqrt((Ellipses{schnitzcells(k).frames(m)-1}(:,1)-...
                                       schnitzcells(n).cenx(IndexToCheck)).^2+...
                                       (Ellipses{schnitzcells(k).frames(m)-1}(:,2)-...
                                       schnitzcells(n).ceny(IndexToCheck)).^2);          

                                    [MinValueK,MinIndexK]=min(DistancesK);
                                    if MinValue<MinValueK
                                        schnitzcells(n).cellno(IndexToCheck)=[];                                        
                                    end
                                end
                          end
                       end
                   end
                end

                if ~EllipseFoundElsewhere
                   schnitzcells(k).cellno(m)=MinIndex;
                else
                   MinValue
                   1+1; error('What to do here?')
                end
           end

       else
           error('Cannnot rescue schnitz')
       end
    end
end
    
