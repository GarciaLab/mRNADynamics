function [Particles,schnitzcells]=SplitSchnitzDaughters(Particles,schnitzcells,...
    SplitFrame,SourceSchnitz,DaughterSchnitz1,DaughterSchnitz2)
% function [Particles,schnitzcells]=SplitSchnitzDaughters(Particles,schnitzcells,...
%   SplitFrame,SourceSchnitz,DaughterSchnitz1,DaughterSchnitz2)
%
% DESCRIPTION
% This function splits a schnitzcell and associates it with daughter
% nuclei. Typically, this is a subfunction for CheckParticleTracking
%
% ARGUMENTS
% Particles: The particles structure created by TrackmRNADynamics and
% located in Particles.mat
% schnitzcells: The list of segmented nuclei with lineages created by
% TrackNuclei and stored in _lin.mat
% SplitFrame: The frame at which the new nuclei have already started.
% SourceSchnitz: The mother schnitz to be split into two daughters
% DaughterSchnitz1: One of the chosen daughter schnitzs
% DaughterSchnitz2: If DaughterSchnitz2 is 0 then we only associate it with one.
%
% OUTPUT
% Particles: A modified Particles
% schnitzcells: A modified schnitzcells
%
% Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
% Created: Unknown
% Last Updated: V2: Updated to support Laurent's schnitzcells
%

%Check if we are splitting a schnitz and creating daughter schnitzes
if (SourceSchnitz==DaughterSchnitz1)|(SourceSchnitz==DaughterSchnitz2)
    
    %Create a filter for the corresponding frames
    FrameFilter=schnitzcells(SourceSchnitz).frames>=SplitFrame;

    %Create the new schnitz
    Nschnitz=length(schnitzcells);
    
    %Assign the daughter information
    if (schnitzcells(SourceSchnitz).E~=0)|(schnitzcells(SourceSchnitz).D~=0)
        if SourceSchnitz==DaughterSchnitz1
            schnitzcells(Nschnitz+1).P=SourceSchnitz;
            
            schnitzcells(SourceSchnitz).E=DaughterSchnitz2;
            schnitzcells(SourceSchnitz).D=Nschnitz+1;
            
            if DaughterSchnitz2~=0
                schnitzcells(DaughterSchnitz2).P=SourceSchnitz;  

                %Check that there isn't another schnitz who thinks it's the
                %parent of DaughterSchnitz2
                schnitzcells=CheckMultipleDaughters(schnitzcells,DaughterSchnitz2,SourceSchnitz);
            end
                
                
        elseif SourceSchnitz==DaughterSchnitz2
            schnitzcells(Nschnitz+1).P=SourceSchnitz;
            
            schnitzcells(SourceSchnitz).E=DaughterSchnitz1;
            schnitzcells(SourceSchnitz).D=Nschnitz+1;
            
            schnitzcells(DaughterSchnitz1).P=SourceSchnitz; 
            
            %Check that there isn't another schnitz who thinks it's the
            %parent of DaughterSchnitz1
            schnitzcells=CheckMultipleDaughters(schnitzcells,DaughterSchnitz1,SourceSchnitz);
        end
        
        %error('I need to deal with this case')
    else
        schnitzcells(Nschnitz+1).P=SourceSchnitz;
        schnitzcells(SourceSchnitz).D=Nschnitz+1;
        schnitzcells(Nschnitz+1).D=0;
        schnitzcells(Nschnitz+1).E=0;
        if SourceSchnitz==DaughterSchnitz1
            schnitzcells(SourceSchnitz).E=DaughterSchnitz2;
            
            if DaughterSchnitz2~=0
                schnitzcells(DaughterSchnitz2).P=SourceSchnitz;

                %Check that there isn't another schnitz who thinks it's the
                %parent of DaughterSchnitz2
                schnitzcells=CheckMultipleDaughters(schnitzcells,DaughterSchnitz2,SourceSchnitz);
            end
            
        elseif SourceSchnitz==DaughterSchnitz2
            schnitzcells(SourceSchnitz).E=DaughterSchnitz1;
            
            if DaughterSchnitz1~=0
                schnitzcells(DaughterSchnitz1).P=SourceSchnitz;

                %Check that there isn't another schnitz who thinks it's the
                %parent of DaughterSchnitz1
                schnitzcells=CheckMultipleDaughters(schnitzcells,DaughterSchnitz1,SourceSchnitz);
            end
        end
    end
    
    
    Fields=fieldnames(schnitzcells);
    for i=1:length(Fields)
        if ~(strcmp(Fields{i},'P')|strcmp(Fields{i},'E')|strcmp(Fields{i},'D'))
            InfoAll=getfield(schnitzcells(SourceSchnitz),Fields{i});
            if length(InfoAll)==length(FrameFilter)
                InfoMother=InfoAll(~FrameFilter);
                InfoDaughter=InfoAll(FrameFilter);

                schnitzcells(Nschnitz+1)=setfield(schnitzcells(Nschnitz+1),Fields{i},InfoDaughter);
                schnitzcells(SourceSchnitz)=setfield(schnitzcells(SourceSchnitz),Fields{i},InfoMother);
            elseif (strcmp(Fields{i},'len')|strcmp(Fields{i},'cellno')) %I had to add this because
                                                                        %some schnitzcells had their ang and length
                                                                        %incomplete. In we leave the field empty
                                                                        %that case
                schnitzcells(Nschnitz+1)=setfield(schnitzcells(Nschnitz+1),Fields{i},[]);
            end
        end
    end
    
    %Now, check if there are any particles that were associated with the
    %new nucleus
    for i=1:length(Particles)
        if Particles(i).Nucleus==SourceSchnitz
            AssociatedParticle=i;
        end
    end
    
    if sum(Particles(AssociatedParticle).Frame==SplitFrame)
        warning('Particle associated with this nucleus. Check that things make sense.')
    end
    
else
    %Here we join schnitzes
    schnitzcells(SourceSchnitz).E=DaughterSchnitz1;
    schnitzcells(SourceSchnitz).D=DaughterSchnitz2;
    
    if DaughterSchnitz1~=0
        schnitzcells(DaughterSchnitz1).P=SourceSchnitz;
        
        %Check that there isn't another schnitz who thinks it's the
        %parent of DaughterSchnitz1
        schnitzcells=CheckMultipleDaughters(schnitzcells,DaughterSchnitz1,SourceSchnitz);
    end
    if DaughterSchnitz2~=0
        schnitzcells(DaughterSchnitz2).P=SourceSchnitz;
        
        %Check that there isn't another schnitz who thinks it's the
        %parent of DaughterSchnitz2
        schnitzcells=CheckMultipleDaughters(schnitzcells,DaughterSchnitz2,SourceSchnitz);
    end
end



