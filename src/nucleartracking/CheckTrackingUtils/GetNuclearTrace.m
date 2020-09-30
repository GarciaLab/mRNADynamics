function [PlotNuclearTraceSettings] =...
    GetNuclearTrace(cntState, PlotNuclearTraceSettings)
% This function uses the total intensity mask to calculate the particles
% intensity and subtracts the background obtained from the fit.

% First, get the different intensity values corresponding to this particle.
CurrentNucleus = cntState.CurrentNucleus;
schnitzcells = cntState.schnitzcells;
NSlices = size(schnitzcells(CurrentNucleus).Fluo, 2);

for i=1:length(schnitzcells(CurrentNucleus).frames)
    
    Frame(i)=schnitzcells(CurrentNucleus).frames(i);
    MaxFluo(i) = max(schnitzcells(CurrentNucleus).Fluo(i,:));
    MaxZ(i) = find(schnitzcells(CurrentNucleus).Fluo(i,:) == max(schnitzcells(CurrentNucleus).Fluo(i,:)), 1);
    MedFluo(i) = median(schnitzcells(CurrentNucleus).Fluo(i,2:NSlices-1));
    MedZ(i) = find(schnitzcells(CurrentNucleus).Fluo(i,:) == median(schnitzcells(CurrentNucleus).Fluo(i,:)), 1);
    MidMedFluo(i) = median(schnitzcells(CurrentNucleus).Fluo(i,max(2, MaxZ(i)-5):min(size(schnitzcells(CurrentNucleus).Fluo, 2)-1, MaxZ(i)+5)));
    if ~isempty(find(schnitzcells(CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1))
        MidMedZ(i) = find(schnitzcells(CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1);
    else
        SubFluos = schnitzcells(CurrentNucleus).Fluo(i,max(2, MaxZ(i)-5):min(size(schnitzcells(CurrentNucleus).Fluo, 2)-1, MaxZ(i)+5));
        SubFluos = sort(SubFluos(SubFluos > MidMedFluo(i)));
        MidMedFluo(i) = SubFluos(1);
        MidMedZ(i) = find(schnitzcells(CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1);
        
    end

end
PlotNuclearTraceSettings.Frames = Frame;
PlotNuclearTraceSettings.MaxFluo = MaxFluo;
PlotNuclearTraceSettings.MedFluo = MedFluo;
PlotNuclearTraceSettings.MidMedFluo = MidMedFluo;
PlotNuclearTraceSettings.MaxZ = MaxZ;
PlotNuclearTraceSettings.MedZ = MedZ;
PlotNuclearTraceSettings.MidMedZ = MidMedZ;








