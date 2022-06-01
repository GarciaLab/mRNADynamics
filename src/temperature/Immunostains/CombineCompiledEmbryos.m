function CT_CompiledEmbryos = CombineCompiledEmbryos(SetPrefixes)
for i = 1:length(SetPrefixes)
    Prefix = SetPrefixes{i};
    liveExperiment = LiveExperiment(Prefix);

    CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
    load(CompiledEmbryoPath);
    if i == 1
        CT_CompiledEmbryos = CompiledEmbryos(1);
        CT_CompiledEmbryos.SetIDs = ones(1, size(CompiledEmbryos.CoordAs, 1));
    else
        CT_CompiledEmbryos.SetIDs = [CT_CompiledEmbryos.SetIDs ones(1, size(CompiledEmbryos.CoordAs, 1))];
        CT_CompiledEmbryos.CoordAs = [CT_CompiledEmbryos.CoordAs ;  CompiledEmbryos.CoordAs];
        CT_CompiledEmbryos.CoordPs = [CT_CompiledEmbryos.CoordPs ; CompiledEmbryos.CoordPs];
        CT_CompiledEmbryos.CoordDs = [CT_CompiledEmbryos.CoordDs ; CompiledEmbryos.CoordDs];
        CT_CompiledEmbryos.CoordVs = [CT_CompiledEmbryos.CoordVs;  CompiledEmbryos.CoordVs];
        CT_CompiledEmbryos.Flags = [CT_CompiledEmbryos.Flags  CompiledEmbryos.Flags];
        CT_CompiledEmbryos.Approved = [CT_CompiledEmbryos.Approved  CompiledEmbryos.Approved];
        CT_CompiledEmbryos.APLengths = [CT_CompiledEmbryos.APLengths  CompiledEmbryos.APLengths];
        CT_CompiledEmbryos.DVLengths = [CT_CompiledEmbryos.DVLengths  CompiledEmbryos.DVLengths];
        
        CT_CompiledEmbryos.APRotationAngles = [CT_CompiledEmbryos.APRotationAngles  CompiledEmbryos.APRotationAngles];
        CT_CompiledEmbryos.FlippedOrientation = [CT_CompiledEmbryos.FlippedOrientation  CompiledEmbryos.FlippedOrientation];
        CT_CompiledEmbryos.APSlopes = [CT_CompiledEmbryos.APSlopes  CompiledEmbryos.APSlopes];
        CT_CompiledEmbryos.APIntercepts = [CT_CompiledEmbryos.APIntercepts  CompiledEmbryos.APIntercepts];
        CT_CompiledEmbryos.Midpoints = [CT_CompiledEmbryos.Midpoints ; CompiledEmbryos.Midpoints];
        CT_CompiledEmbryos.DVSlopes = [CT_CompiledEmbryos.DVSlopes  CompiledEmbryos.DVSlopes];
        CT_CompiledEmbryos.DVIntercepts = [CT_CompiledEmbryos.DVIntercepts  CompiledEmbryos.DVIntercepts];
        
        CT_CompiledEmbryos.Checked = [CT_CompiledEmbryos.Checked  CompiledEmbryos.Checked];
        CT_CompiledEmbryos.xShift = [CT_CompiledEmbryos.xShift  CompiledEmbryos.xShift];
        CT_CompiledEmbryos.yShift = [CT_CompiledEmbryos.yShift  CompiledEmbryos.yShift];
        
        CT_CompiledEmbryos.RotatedCoordAs = [CT_CompiledEmbryos.RotatedCoordAs ; CompiledEmbryos.RotatedCoordAs];
        CT_CompiledEmbryos.RotatedCoordPs = [CT_CompiledEmbryos.RotatedCoordPs ; CompiledEmbryos.RotatedCoordPs];
        CT_CompiledEmbryos.RotatedCoordDs = [CT_CompiledEmbryos.RotatedCoordDs ; CompiledEmbryos.RotatedCoordDs];
        CT_CompiledEmbryos.RotatedCoordVs = [CT_CompiledEmbryos.RotatedCoordVs ; CompiledEmbryos.RotatedCoordVs];
        
        CT_CompiledEmbryos.nc = [CT_CompiledEmbryos.nc  CompiledEmbryos.nc];
        
        CT_CompiledEmbryos.FurrowMeasurements = [CT_CompiledEmbryos.FurrowMeasurements  CompiledEmbryos.FurrowMeasurements];
        CT_CompiledEmbryos.CellWidthMeasurements = [CT_CompiledEmbryos.CellWidthMeasurements  CompiledEmbryos.CellWidthMeasurements];
        CT_CompiledEmbryos.ZoomFurrowMeasurements = [CT_CompiledEmbryos.ZoomFurrowMeasurements  CompiledEmbryos.ZoomFurrowMeasurements];
        CT_CompiledEmbryos.ZoomCellWidthMeasurements = [CT_CompiledEmbryos.ZoomCellWidthMeasurements  CompiledEmbryos.ZoomCellWidthMeasurements];
    end
end