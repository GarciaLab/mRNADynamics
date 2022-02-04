function Ellipses = UpdateEllipsesOriginalCoordinates(Ellipses,Prefix, CompiledEmbryos,liveExperiment)
%%
nEmbryos = length(Ellipses);
for i = 1:nEmbryos
    if ~isempty(Ellipses{i})
        Ellipses{i}(:,5:6) = TransformToOriginalImageCoordinates( Ellipses{i}(:,1:2), Prefix, CompiledEmbryos,liveExperiment, i);
        Ellipses{i}(:,3:4) = max(Ellipses{i}(:,3));
    end
end