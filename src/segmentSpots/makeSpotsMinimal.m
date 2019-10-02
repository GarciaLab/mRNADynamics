function [SpotsMinimal, SpotsExtra, Spots3D] = makeSpotsMinimal(Spots)

%Makes a smaller version of Spots. Necessary fields are in SpotsMinimal. Less
%important fields are in SpotsExtra. 3D information is stored in Spots3D.

SpotsExtra = Spots;
Spots3D = Spots;
nFrames = length(Spots);

for t = 1:nFrames
    
    %generate the structure with the most important fields: Spots
    SpotsMinimal(t).Fits = rmfield([Spots(t).Fits], 'ConfidenceIntervals');
    
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'gaussParams');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'xFitWidth');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'yFitWidth');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'FixedAreaIntensity3');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'GaussianIntensity');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'CentralIntensity');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'dogFixedAreaIntensity');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'intArea');
    SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'DOGIntensity');
    
    if isfield(SpotsMinimal(t).Fits, 'Spot1Fits3D')
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot1Fits3D');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot1CI3D');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot1Int3D');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot1Pos');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot2Fits3D');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot2CI3D');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot2Int3D');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Spot2Pos');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'Offset3D');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'OffsetCI');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'gauss3DIntensity');
        SpotsMinimal(t).Fits = rmfield([SpotsMinimal(t).Fits], 'GaussPos');

    end
    
    %generate structure with the less important fields: Spots2
    SpotsExtra(t).Fits = rmfield([Spots(t).Fits], 'FixedAreaIntensity');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'xFit');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'yFit');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'z');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Offset');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'brightestZ');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'xDoG');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'yDoG');
    SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'snippet_size');
        
   if isfield(SpotsExtra(t).Fits, 'Spot1Fits3D')
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot1Fits3D');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot1CI3D');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot1Int3D');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot1Pos');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot2Fits3D');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot2CI3D');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot2Int3D');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Spot2Pos');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'Offset3D');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'OffsetCI');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'gauss3DIntensity');
        SpotsExtra(t).Fits = rmfield([SpotsExtra(t).Fits], 'GaussPos');
   end


    
    %generate the structure with the 3D gaussian fields: Spots3D
    Spots3D(t).Fits = rmfield([Spots(t).Fits], 'ConfidenceIntervals');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'gaussParams');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'xFitWidth');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'yFitWidth');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'FixedAreaIntensity3');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'GaussianIntensity');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'CentralIntensity');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'dogFixedAreaIntensity');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'intArea');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'DOGIntensity');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'GaussPos');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'FixedAreaIntensity');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'xFit');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'yFit');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'z');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'Offset');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'brightestZ');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'xDoG');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'yDoG');
    Spots3D(t).Fits = rmfield([Spots3D(t).Fits], 'snippet_size');
    
    
end

end