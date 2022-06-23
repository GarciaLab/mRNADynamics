function CT_CompiledEmbryos = CombineCompiledEmbryos(SetPrefixes)
%%

for i = 1:length(SetPrefixes)
    Prefix = SetPrefixes{i};
    liveExperiment = LiveExperiment(Prefix);

    CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
    load(CompiledEmbryoPath);

    CompiledEmbryos = AddDorsalProfInfoToCompiledEmbryos(Prefix);

    if isfield(CompiledEmbryos, 'ZoomFurrowMeasurements')
        NEmbryos = length(CompiledEmbryos.Approved);
        CompiledEmbryos.DeltaFC_um = {};
        CompiledEmbryos.DeltaFC_um.mean = NaN(1, NEmbryos);
        CompiledEmbryos.DeltaFC_um.std = NaN(1, NEmbryos);
        CompiledEmbryos.DeltaFC_um.count = NaN(1, NEmbryos);
        CompiledEmbryos.DeltaFC_um.se = NaN(1, NEmbryos);
        CompiledEmbryos.UseCustomFurrowMeasurements = zeros(1, NEmbryos, 'logical');
        for j = 1:NEmbryos
            if size(CompiledEmbryos.ZoomFurrowMeasurements{j}, 1) >= 3
                CompiledEmbryos.DeltaFC_um.mean(j) = mean(CompiledEmbryos.ZoomFurrowMeasurements{j}(:,10));
                CompiledEmbryos.DeltaFC_um.std(j) = std(CompiledEmbryos.ZoomFurrowMeasurements{j}(:,10));
                CompiledEmbryos.DeltaFC_um.count(j) = size(CompiledEmbryos.ZoomFurrowMeasurements{j}, 1);
                CompiledEmbryos.DeltaFC_um.se(j) = CompiledEmbryos.DeltaFC_um.std(j)/sqrt(CompiledEmbryos.DeltaFC_um.count(j) );
            elseif size(CompiledEmbryos.FurrowMeasurements{j}, 1) >= 3
                CompiledEmbryos.DeltaFC_um.mean(j) = mean(CompiledEmbryos.FurrowMeasurements{j}(:,10));
                CompiledEmbryos.DeltaFC_um.std(j) = std(CompiledEmbryos.FurrowMeasurements{j}(:,10));
                CompiledEmbryos.DeltaFC_um.count(j) = size(CompiledEmbryos.FurrowMeasurements{j}, 1);
                CompiledEmbryos.DeltaFC_um.se(j) = CompiledEmbryos.DeltaFC_um.std(j)/sqrt(CompiledEmbryos.DeltaFC_um.count(j) );
                CompiledEmbryos.UseCustomFurrowMeasurements(j) = true;
            end
        end
        save(CompiledEmbryoPath, 'CompiledEmbryos');
    end
    
    if i == 1
        CT_CompiledEmbryos = CompiledEmbryos(1);
        CT_CompiledEmbryos.SlideIDs = ones(1, size(CompiledEmbryos.CoordAs, 1));
    else
        CT_CompiledEmbryos.SlideIDs = [CT_CompiledEmbryos.SlideIDs i*ones(1, size(CompiledEmbryos.CoordAs, 1))];
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
        
        CT_CompiledEmbryos.MemCoordAs = [CT_CompiledEmbryos.MemCoordAs ; CompiledEmbryos.MemCoordAs];
        CT_CompiledEmbryos.MemCoordPs = [CT_CompiledEmbryos.MemCoordPs ; CompiledEmbryos.MemCoordPs];
        CT_CompiledEmbryos.MemCoordDs = [CT_CompiledEmbryos.MemCoordDs ; CompiledEmbryos.MemCoordDs];
        CT_CompiledEmbryos.MemCoordVs = [CT_CompiledEmbryos.MemCoordVs ; CompiledEmbryos.MemCoordVs];
        
        CT_CompiledEmbryos.MemAPLengths = [CT_CompiledEmbryos.MemAPLengths  CompiledEmbryos.MemAPLengths];
        CT_CompiledEmbryos.MemDVLengths = [CT_CompiledEmbryos.MemDVLengths  CompiledEmbryos.MemDVLengths];
        
        CT_CompiledEmbryos.MemAPRotationAngles = [CT_CompiledEmbryos.MemAPRotationAngles  CompiledEmbryos.MemAPRotationAngles];
        CT_CompiledEmbryos.MemAPSlopes = [CT_CompiledEmbryos.MemAPSlopes  CompiledEmbryos.MemAPSlopes];
        CT_CompiledEmbryos.MemAPIntercepts = [CT_CompiledEmbryos.MemAPIntercepts  CompiledEmbryos.MemAPIntercepts];
        CT_CompiledEmbryos.MemMidpoints = [CT_CompiledEmbryos.MemMidpoints ; CompiledEmbryos.MemMidpoints];
        CT_CompiledEmbryos.MemDVSlopes = [CT_CompiledEmbryos.MemDVSlopes  CompiledEmbryos.MemDVSlopes];
        CT_CompiledEmbryos.MemDVIntercepts = [CT_CompiledEmbryos.MemDVIntercepts  CompiledEmbryos.MemDVIntercepts];
        
        CT_CompiledEmbryos.MemxShift = [CT_CompiledEmbryos.MemxShift  CompiledEmbryos.MemxShift];
        CT_CompiledEmbryos.MemyShift = [CT_CompiledEmbryos.MemyShift  CompiledEmbryos.MemyShift];
        
        CT_CompiledEmbryos.MemRotatedCoordAs = [CT_CompiledEmbryos.MemRotatedCoordAs ; CompiledEmbryos.MemRotatedCoordAs];
        CT_CompiledEmbryos.MemRotatedCoordPs = [CT_CompiledEmbryos.MemRotatedCoordPs ; CompiledEmbryos.MemRotatedCoordPs];
        CT_CompiledEmbryos.MemRotatedCoordDs = [CT_CompiledEmbryos.MemRotatedCoordDs ; CompiledEmbryos.MemRotatedCoordDs];
        CT_CompiledEmbryos.MemRotatedCoordVs = [CT_CompiledEmbryos.MemRotatedCoordVs ; CompiledEmbryos.MemRotatedCoordVs];
        
        CT_CompiledEmbryos.FurrowMeasurements = [CT_CompiledEmbryos.FurrowMeasurements  CompiledEmbryos.FurrowMeasurements];
        CT_CompiledEmbryos.CellWidthMeasurements = [CT_CompiledEmbryos.CellWidthMeasurements  CompiledEmbryos.CellWidthMeasurements];
        CT_CompiledEmbryos.ZoomFurrowMeasurements = [CT_CompiledEmbryos.ZoomFurrowMeasurements  CompiledEmbryos.ZoomFurrowMeasurements];
        CT_CompiledEmbryos.ZoomCellWidthMeasurements = [CT_CompiledEmbryos.ZoomCellWidthMeasurements  CompiledEmbryos.ZoomCellWidthMeasurements];
        
        CT_CompiledEmbryos.SingleRowDorsalNucleiCounts = [CT_CompiledEmbryos.SingleRowDorsalNucleiCounts  CompiledEmbryos.SingleRowDorsalNucleiCounts];
        CT_CompiledEmbryos.SingleRowDorsalMedianNucleiSpacing = [CT_CompiledEmbryos.SingleRowDorsalMedianNucleiSpacing  CompiledEmbryos.SingleRowDorsalMedianNucleiSpacing];
        
        CT_CompiledEmbryos.DorsalAPProfiles = [CT_CompiledEmbryos.DorsalAPProfiles ; CompiledEmbryos.DorsalAPProfiles]; 
        CT_CompiledEmbryos.DorsalAvgAPProfiles = [CT_CompiledEmbryos.DorsalAvgAPProfiles ; CompiledEmbryos.DorsalAvgAPProfiles]; 
        CT_CompiledEmbryos.NarrowProfileNucleiFluoInfo = [CT_CompiledEmbryos.NarrowProfileNucleiFluoInfo ; CompiledEmbryos.NarrowProfileNucleiFluoInfo]; 
        CT_CompiledEmbryos.DorsalNarrowAPProfiles = [CT_CompiledEmbryos.DorsalNarrowAPProfiles ; CompiledEmbryos.DorsalNarrowAPProfiles]; 
        CT_CompiledEmbryos.DorsalAvgNarrowAPProfiles = [CT_CompiledEmbryos.DorsalAvgNarrowAPProfiles ; CompiledEmbryos.DorsalAvgNarrowAPProfiles]; 
        CT_CompiledEmbryos.ProfileNucleiFluoInfo = [CT_CompiledEmbryos.ProfileNucleiFluoInfo ; CompiledEmbryos.ProfileNucleiFluoInfo]; 
        CT_CompiledEmbryos.ProfileNarrowNucleiFluoInfo = [CT_CompiledEmbryos.ProfileNarrowNucleiFluoInfo  CompiledEmbryos.ProfileNarrowNucleiFluoInfo]; 
        CT_CompiledEmbryos.AllDorsalNucleiFluoInfo = [CT_CompiledEmbryos.AllDorsalNucleiFluoInfo ; CompiledEmbryos.AllDorsalNucleiFluoInfo]; 
        CT_CompiledEmbryos.DorsalStdAPProfiles = [CT_CompiledEmbryos.DorsalStdAPProfiles ; CompiledEmbryos.DorsalStdAPProfiles]; 
        CT_CompiledEmbryos.DorsalStdNarrowAPProfiles = [CT_CompiledEmbryos.DorsalStdNarrowAPProfiles ; CompiledEmbryos.DorsalStdNarrowAPProfiles]; 
        CT_CompiledEmbryos.DorsalCountAPProfiles = [CT_CompiledEmbryos.DorsalCountAPProfiles ; CompiledEmbryos.DorsalCountAPProfiles]; 
        CT_CompiledEmbryos.DorsalCountNarrowAPProfiles = [CT_CompiledEmbryos.DorsalCountNarrowAPProfiles ; CompiledEmbryos.DorsalCountNarrowAPProfiles]; 
        
        CT_CompiledEmbryos.DeltaFC_um.mean = [CT_CompiledEmbryos.DeltaFC_um.mean  CompiledEmbryos.DeltaFC_um.mean]; 
        CT_CompiledEmbryos.DeltaFC_um.std = [CT_CompiledEmbryos.DeltaFC_um.std  CompiledEmbryos.DeltaFC_um.std]; 
        CT_CompiledEmbryos.DeltaFC_um.count = [CT_CompiledEmbryos.DeltaFC_um.count  CompiledEmbryos.DeltaFC_um.count]; 
        CT_CompiledEmbryos.DeltaFC_um.se = [CT_CompiledEmbryos.DeltaFC_um.se  CompiledEmbryos.DeltaFC_um.se]; 
        CT_CompiledEmbryos.UseCustomFurrowMeasurements = [CT_CompiledEmbryos.UseCustomFurrowMeasurements CompiledEmbryos.UseCustomFurrowMeasurements];
        
        
    end
end