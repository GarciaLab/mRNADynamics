function [Frames,Amp]=PlotParticleTrace(CurrentParticle,Particles,fad,DataFolder,FilePrefix)

%This displays the particle trace as well as a montage of the images

%V2: Removed the ImageSnippet for now. I was running into the borders of
%the image.

NRows=3;            %Number of rows used to display the images
ImageSize=20;       %Size of the image around the found particle. The total
                    %image size will be 2*ImageSize+1.
MarginSize=1;       %Size of the white margin around the snippets


[Frame,Amp]=GetParticleTrace(CurrentParticle,Particles,fad);
Frames=Particles(CurrentParticle).Frame;
Indexes=Particles(CurrentParticle).Index;



for i=1:length(Frames)
    CurrentFrame=Frames(i);

    %Get the coordinates taking the margins into account
    [x,y]=fad2xyzFit(CurrentFrame,fad, 'addMargin'); 
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex=Particles(CurrentParticle).Index(Particles(CurrentParticle).Frame==CurrentFrame);
    
    %This is the position of the current particle
    xTrace=round(x(CurrentParticleIndex));
    yTrace=round(y(CurrentParticleIndex));
    zTrace=fad.channels(CurrentFrame).fits.z(CurrentParticleIndex);

   
    %Image=imread([DataFolder,'\Data\',FilePrefix(1:end-1),filesep,FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(zTrace,2),'.tif']);
   
    %ImageSnippet(:,:,1,i)=ones(2*ImageSize+1+2*MarginSize,2*ImageSize+1+2*MarginSize);
    
    %ImageSnippet(MarginSize+1:end-MarginSize,MarginSize+1:end-MarginSize,1,i)=mat2gray(Image(yTrace-ImageSize:yTrace+ImageSize,...
    %    xTrace-ImageSize:xTrace+ImageSize));
end



% figure(1)
% subplot(1,2,1)
% plot(Frames,Raw,'.k')
% subplot(1,2,2)
% montage(ImageSnippet,'DisplayRange',[])