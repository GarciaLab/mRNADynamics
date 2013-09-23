function [fad,fad2]=...
    TransferParticleBack(CurrentFrame,ParticleIndex,fad,fad2)

%This function is inspired on TransferParticle.m. The idea is that we will
%move a particle from fad to fad2

%This function moves a particle from the structure fad2 to fad. It
%also updates the Particles structure.

%If Particlefad is empty it adds the particle from fad2 to fad and at
%the end of Particles.

Fields=fieldnames(fad.channels(CurrentFrame).fits);
NParticles=length(fad.channels(CurrentFrame).fits.dog);

%Create a filter to get the particle out
FilterParticle=logical(zeros(size(fad.channels(CurrentFrame).fits.dog)));
FilterParticle(ParticleIndex)=1;


for j=1:length(Fields)
    %In the end this seems to work, the fields that have higher
    %dimensions are still devided well.
    
    if strcmp(Fields{j},'maskUsedForTotalInt')
        %Do nothing with this one
            
    elseif strcmp(Fields{j},'snippets')
        
        Temp=getfield(fad.channels(CurrentFrame).fits,Fields{j});
        Temp1=Temp(:,:,FilterParticle);
        Temp2=Temp(:,:,~FilterParticle);
           
        fad2.channels(CurrentFrame).fits=setfield(fad2.channels(CurrentFrame).fits,Fields{j},...
               cat(3,getfield(fad2.channels(CurrentFrame).fits,Fields{j}),Temp1));
    
        
        fad.channels(CurrentFrame).fits=setfield(fad.channels(CurrentFrame).fits,Fields{j},...
            Temp2);
        
    elseif length(getfield(fad.channels(CurrentFrame).fits,Fields{j}))==NParticles
        
        Temp=getfield(fad.channels(CurrentFrame).fits,Fields{j});
        Temp1=Temp(FilterParticle);
        Temp2=Temp(~FilterParticle);
        
        [NRows,Dummy]=size(getfield(fad2.channels(CurrentFrame).fits,Fields{j}));
        
        %I had to add this because some entries are row vectors and some
        %other ones are column vectors.
        [NRows2,NCols2]=size(Temp1);
        
        if NRows>1
            if (NCols2>1)
                Temp1=Temp1';
            end
            fad2.channels(CurrentFrame).fits=setfield(fad2.channels(CurrentFrame).fits,Fields{j},...
                [getfield(fad2.channels(CurrentFrame).fits,Fields{j});Temp1]);
        else
            if (NRows2>1)
                Temp1=Temp1';
            end
            
            fad2.channels(CurrentFrame).fits=setfield(fad2.channels(CurrentFrame).fits,Fields{j},...
                [getfield(fad2.channels(CurrentFrame).fits,Fields{j}),Temp1]);
        end
        
        fad.channels(CurrentFrame).fits=setfield(fad.channels(CurrentFrame).fits,Fields{j},...
            Temp2);
        
    else
        Fields{j}
        error('Field not processed in TransferParticle.m')
    end
end


    


