function schnitzcells=ExtractNuclearFluorescence(schnitzcells,CurrentFrame,...
    Image,LinesPerFrame,PixelsPerLine,NumberSlices,Circle,IntegrationRadius)


%Find the schnitzs in the current image and extract the corresponding
%fluorescence values. This function used to exist inside TrackNuclei, but I
%had to make it independent so that I could use parfor loops.


%Create a blank image we'll use to generate the mask
Mask=logical(zeros(LinesPerFrame,PixelsPerLine));

%Get the fluorescence over the mask
if sum(schnitzcells.frames==CurrentFrame)
    CurrentIndex=find(schnitzcells.frames==CurrentFrame);
    cenx=schnitzcells.cenx(CurrentIndex);
    ceny=schnitzcells.ceny(CurrentIndex);
    Radius=schnitzcells.len(CurrentIndex);

    %Check that the nucleus and the nuclear mask fit within the image
    if ((cenx-Radius)>0&(cenx+Radius)<PixelsPerLine&...
            (ceny-Radius)>0&(ceny+Radius)<LinesPerFrame)&...
            ((round(cenx)-(3*IntegrationRadius-1)/2)>0&...
            (round(cenx)+(3*IntegrationRadius-1)/2)<PixelsPerLine&...
            (round(ceny)-(3*IntegrationRadius-1)/2)>0&...
            (round(ceny)+(3*IntegrationRadius-1)/2)<LinesPerFrame)

        %Now, add the circle
        Mask(round(ceny)-(3*IntegrationRadius-1)/2:...
            round(ceny)+(3*IntegrationRadius-1)/2,...
            round(cenx)-(3*IntegrationRadius-1)/2:...
            round(cenx)+(3*IntegrationRadius-1)/2)=Circle;

        %Save the mask we're going to use
        schnitzcells.Mask=Circle;


        for CurrentZ=1:(NumberSlices+2)
            schnitzcells.Fluo(CurrentIndex,CurrentZ)=sum(sum(immultiply(Image(:,:,CurrentZ),Mask)));
        end

    else  %If not assign NaN to the fluroescence
        schnitzcells.Fluo(CurrentIndex,1:(NumberSlices+2))=nan;
    end
end