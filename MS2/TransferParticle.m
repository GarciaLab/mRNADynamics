function [fad,fad2,Particles]=...
    TransferParticle(fad,Framefad2,Particlefad,fad2,Framefad,...
    Indexfad2,Particles)

%This function moves a particle from the structure fad2 to fad. It
%also updates the Particles structure.

%If Particlefad is empty it adds the particle from fad2 to fad and at
%the end of Particles.

Fields=fieldnames(fad2.channels(Framefad2).fits);
NParticles=length(fad2.channels(Framefad2).fits.dog);

%Create a filter to get the particle out
FilterParticle=logical(zeros(size(fad2.channels(Framefad2).fits.dog)));
FilterParticle(Indexfad2)=1;


for j=1:length(Fields)
    %In the end this seems to work, the fields that have higher
    %dimensions are still devided well.
    
    if strcmp(Fields{j},'maskUsedForTotalInt')
        %Do nothing with this one
            
    elseif strcmp(Fields{j},'snippets')
        
        Temp=getfield(fad2.channels(Framefad2).fits,Fields{j});
        Temp1=Temp(:,:,FilterParticle);
        Temp2=Temp(:,:,~FilterParticle);
           
        fad.channels(Framefad2).fits=setfield(fad.channels(Framefad2).fits,Fields{j},...
               cat(3,getfield(fad.channels(Framefad2).fits,Fields{j}),Temp1));
    
        
        fad2.channels(Framefad2).fits=setfield(fad2.channels(Framefad2).fits,Fields{j},...
            Temp2);
        
    elseif length(getfield(fad2.channels(Framefad2).fits,Fields{j}))==NParticles
        
        Temp=getfield(fad2.channels(Framefad2).fits,Fields{j});
        Temp1=Temp(FilterParticle);
        Temp2=Temp(~FilterParticle);
        
        [NRows,Dummy]=size(getfield(fad.channels(Framefad2).fits,Fields{j}));
        
        %I had to add this because some entries are row vectors and some
        %other ones are column vectors.
        if NRows>1
            fad.channels(Framefad2).fits=setfield(fad.channels(Framefad2).fits,Fields{j},...
                [getfield(fad.channels(Framefad2).fits,Fields{j});Temp1]);
        else
            fad.channels(Framefad2).fits=setfield(fad.channels(Framefad2).fits,Fields{j},...
                [getfield(fad.channels(Framefad2).fits,Fields{j}),Temp1]);
        end
        
        fad2.channels(Framefad2).fits=setfield(fad2.channels(Framefad2).fits,Fields{j},...
            Temp2);
        
    else
        Fields{j}
        error('Field not processed in TransferParticle.m')
    end
end


%Now I need to insert the information about this new frame of the particle
%in the Particles structure
if ~isempty(Particlefad)
    if Particles(Particlefad).Frame(end)<Framefad2
        %If the new frame is at the beginning
        Particles(Particlefad).Frame(end+1)=Framefad2;
        Particles(Particlefad).Index(end+1)=length(fad.channels(Framefad2).fits.dog);
    elseif Particles(Particlefad).Frame(1)>Framefad2
        %If the new frame is at the end
        Particles(Particlefad).Frame=[Framefad2,Particles(Particlefad).Frame];
        Particles(Particlefad).Index=[length(fad.channels(Framefad2).fits.dog),...
            Particles(Particlefad).Index];
    else
        %This is in case the new frame corresponds to a gap
        [MinValue,MinIndex]=min((Particles(Particlefad).Frame-Framefad2).^2);
        if Particles(Particlefad).Frame(MinIndex)>Framefad2
            Particles(Particlefad).Frame=[Particles(Particlefad).Frame(1:MinIndex-1),...
                Framefad2,Particles(Particlefad).Frame(MinIndex:end)];
            Particles(Particlefad).Index=[Particles(Particlefad).Index(1:MinIndex-1),...
                length(fad.channels(Framefad2).fits.dog),...
                Particles(Particlefad).Index(MinIndex:end)];
        elseif Particles(Particlefad).Frame(MinIndex)<Framefad2
            Particles(Particlefad).Frame=[Particles(Particlefad).Frame(1:MinIndex),...
                Framefad2,Particles(Particlefad).Frame((MinIndex+1):end)];
            Particles(Particlefad).Index=[Particles(Particlefad).Index(1:MinIndex),...
                length(fad.channels(Framefad2).fits.dog),...
                Particles(Particlefad).Index((MinIndex+1):end)];            
        else
            error('Something wrong here')
        end
        
    end
else
    NParticles=length(Particles);
    Particles(NParticles+1).Frame=Framefad2;
    Particles(NParticles+1).Index=length(fad.channels(Framefad2).fits.dog);
    Particles(NParticles+1).Approved=logical(0);
    Particles(NParticles+1).FrameApproved=logical(1);
end

    

    


