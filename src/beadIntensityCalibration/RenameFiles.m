% ExtractBeadIntensities.m
% author: Gabriella Martini
% date created: 1/4/20


rootFolder = 'E:\Gabriella\LivemRNA\Data\RawDynamicsData\2019-12-09\GreenBeadIntensityCalibrationMeasurements';
rawDataDir = dir([rootFolder, filesep, 'GreenCal*']);
newdirnames = {};
expr1A = 'GreenCal(?<size1>\d+)_(?<size2>\d)um';
expr1B = 'GreenCal(?<size1>\d+)um';
expr2A = '^(?<BeadInt>\d+)Set(?<SetID>[A-Za-z])';
expr2B = '^_(?<BeadInt>\d+)Set(?<SetID>[A-Za-z])';
expr2C = 'ControlSet(?<SetID>[A-Za-z])';
expr3A = '^(?<LaserInt>\d+)Percent';
expr3B = '^(?<LaserInt1>\d+)_(?<LaserInt2>\d+)Percent';
for i=1:length(rawDataDir)
    tempname = rawDataDir(i).name;
    tempnamesplit = strsplit(tempname, '-');
    split1 = tempnamesplit{1};
    split2 = tempnamesplit{2};
    split3 = tempnamesplit{3};
    beadsize = regexp(split1, expr1A, 'names');
    if length(beadsize)>0
        if length(beadsize.size1)== 1
            newsize = ['0', beadsize.size1, beadsize.size2];
        elseif length(beadsize.size1) == 2
            newsize = [beadsize.size1, beadsize.size2];
        end
    else
        beadsize = regexp(split1, expr1B, 'names');
        if length(beadsize)>0
            newsize = ['0', beadsize.size1, '0'];
        else
            warning(['Filename for i = ', num2str(i),' does not have expected format'])
            continue
        end
    end

    % Extract Bead Intensities and SetID
    info2A = regexp(split2, expr2A, 'names');
    info2B = regexp(split2, expr2B, 'names');
    info2C = regexp(split2, expr2C, 'names');
    if length(info2A) > 0
        if length(info2A.BeadInt) == 1
            newBeadInt = ['00', info2A.BeadInt, '0'];
        elseif length(info2A.BeadInt) == 2
            newBeadInt = ['0', info2A.BeadInt, '0'];
        elseif length(info2A.BeadInt) == 3
            newBeadInt = [info2A.BeadInt, '0'];
        else
            warning(['Unexpected Bead Intensity for i = ', num2str(i)])
            continue
        end
        SetID = info2A.SetID;
    elseif length(info2B) > 0
        if length(info2B.BeadInt) == 1
            newBeadInt = ['000', info2B.BeadInt];
        else 
            warning(['Unexpected Bead Intensity for i = ', num2str(i)])
            continue
        end
        SetID = info2B.SetID;
    elseif length(info2C) > 0
        newBeadInt = ['0000'];
        setID = info2C.SetID;
    end

    info3A = regexp(split3, expr3A, 'names');
    info3B = regexp(split3, expr3B, 'names');
    if length(info3A) > 0
        tempLaser = info3A.LaserInt;
        if length(tempLaser) == 1
            newLaser = ['00', tempLaser, '0'];
        elseif length(tempLaser) == 2
            newLaser = ['0', tempLaser, '0'];
        elseif length(tempLaser) == 3
            newLaser = [tempLaser, '0'];
        else
            warning(['Unexpected Laser Power for i = ', num2str(i)])
            continue
        end
    elseif length(info3B) > 0
        tempLaser = info3B.LaserInt1;
        if length(tempLaser) == 1
            newLaser = ['00', tempLaser, info3B.LaserInt2];
        elseif length(tempLaser) == 2
            newLaser = ['0', tempLaser, info3B.LaserInt2];
        elseif length(tempLaser) == 3
            newLaser = [tempLaser, info3B.LaserInt2];
        else
            warning(['Unexpected Laser Power for i = ', num2str(i)])
            continue
        end

    else
        warning(['Unexpected Laser Power for i = ', num2str(i)])
        continue
    end

    newdirname = ['488nmLaser-',newLaser, 'Power_Green', newsize, '00nm-', newBeadInt, 'Bead_Set', SetID];
    newfilename = ['2019-12-09-', newdirname, '.lif'];
    
    newdirnames{i} = newdirname;
    APath = fullfile(rootFolder, tempname);
    BPath = fullfile(rootFolder, newdirname);

    subdir = dir([APath, filesep, '*.lif']);
    if length(subdir) == 1
        movefile(APath, fullfile(rootFolder, newdirname));
        CPath = fullfile(BPath, subdir(1).name);
        DPath = fullfile(BPath, newfilename);
        try
            movefile(CPath, DPath)
        catch 
            warning(['i = ', num2str(i)]);
        end
    else 
        warning(['i = ', num2str(i)]);
    end
end
