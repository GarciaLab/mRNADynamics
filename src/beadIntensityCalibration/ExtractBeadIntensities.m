% ExtractBeadIntensities.m
% author: Gabriella Martini
% date created: 1/4/20

clear all, close all
rootFolder = 'E:\Gabriella\LivemRNA\Data\RawDynamicsData\2019-12-09\GreenBeadIntensityCalibrationMeasurements';
rawDataDir = dir([rootFolder, filesep, '*Green*']);
rawDataDir = rawDataDir([rawDataDir.isdir]);
expr4 = '(?<lambda>\d+)nmLaser';
expr5 = '(?<power>\d+)Power_Green(?<size>\d+)nm';
expr6 = '(?<beadintensity>\d+)Bead_Set(?<SetID>[A-Za-z]).lif';
varnames = {'BeadIntensity','LaserPower', 'SetID', 'Replicate', 'z',...
    'BeadID','CenterRow', 'CenterCol','Radius','Intensity'};
vartypes = {'double', 'double', 'string','int64', 'int64', 'int64', 'double',...
    'double', 'double', 'double'};
nvars = length(varnames);
BeadData = table('Size', [0 nvars],...
    'VariableTypes', vartypes,'VariableNames', varnames);
varnames2 = varnames;
varnames2{5} = 'zcount';
MeanBeadData = table('Size', [0 nvars],...
    'VariableTypes', vartypes,'VariableNames', varnames2);


varnames3 = {'BeadIntensity','LaserPower', 'SetID', 'Replicate','MeanIntensity',...
    'StdIntensity', 'BeadCount'};
vartypes3 = {'double', 'double', 'string', 'int64', 'double', 'double', 'int64'};
nvars3 = length(varnames3);
SummaryData = table('Size', [0 nvars3],...
    'VariableTypes', vartypes3,'VariableNames', varnames3);

for i=30:30%1:length(rawDataDir)
display(['i = ', num2str(i)]);
DataDir = dir([rawDataDir(i).folder, filesep, rawDataDir(i).name, filesep, '*.lif']);
LIFPath = fullfile(DataDir(1).folder, DataDir(1).name);
filename = DataDir(1).name;
fnsplit = strsplit(filename, '-');
year = str2num(fnsplit{1});
month = str2num(fnsplit{2});
day = str2num(fnsplit{3});
fn4 = regexp(fnsplit{4}, expr4, 'names');
fn5 = regexp(fnsplit{5}, expr5, 'names');
fn6 = regexp(fnsplit{6}, expr6, 'names');
power = str2num(fn5.power)/10;
beadsize = str2num(fn5.size)/1000; % units microns
lambda = str2num(fn4.lambda);
beadintensity = str2num(fn6.beadintensity)/10;
setID = fn6.SetID;
r = bfGetReader();
% Decorate the reader with the Memoizer wrapper
r = loci.formats.Memoizer(r);
r.setId(LIFPath);
LIFImages = bfopen(LIFPath);
LIFMeta = LIFImages{:, 4};
r.close();
PixelSize = double(LIFMeta.getPixelsPhysicalSizeX(0).value);% units: microns
radius = beadsize/PixelSize/2;
sigma = radius/5;
radmin = round(radius*1/2, 0);
radmax = round(radius, 0);
rep = 1;
for j =1:size(LIFImages, 1)

     NSlices = length(LIFImages{j,1});

    %images = zeros([size(LIFImages{1}{1, 1}, 1), size(LIFImages{1}{1, 1}, 2), NSlices]);
    rn = 1;
    TempBeadData = table('Size', [0 nvars],...
    'VariableTypes', vartypes,'VariableNames', varnames);

    TempMeanBeadData = table('Size', [0 nvars],...
    'VariableTypes', vartypes,'VariableNames', varnames2);

    MaxBeadID = 1;
    for k=1:NSlices
        %images(:,:,k) = LIFImages{j}{k,1};
        img = LIFImages{j}{k,1};
        [centers, radii] = imfindcircles(img, [radmin, radmax],...
            'Sensitivity', 0.9, 'Method','twostage');
        if length(radii) > 1
            keepcenters = ones(length(radii), 1);
            for ic =1:(length(radii)-1)
                for jc=(ic+1):length(radii)
                    d = sqrt((centers(ic, 1)-centers(jc, 1))^2+(centers(ic, 2)-centers(jc,2))^2);
                    if d <= radii(ic)+radii(jc)
                        keepcenters(jc) = 0;
                    end
                end
            end
            centers = centers(keepcenters == 1, :);
            radii = radii(keepcenters == 1);
        end

        if length(radii) == 0
            continue
        else
            I_blurred = imgaussfilt(img,sigma, 'FilterSize', ceil(2*sigma));
            for l=1:length(radii)
                [xgrid, ygrid] = meshgrid(1:size(I_blurred,1), 1:size(I_blurred,2));
                mask = ((xgrid-centers(l,1)).^2 + (ygrid-centers(l,2)).^2) <= (sigma).^2;
                values = I_blurred(mask);
                intensity = mean(values);
                %intensity = I_blurred(round(centers(l,1)), round(centers(l,2)));
                if size(TempBeadData, 1) ~= 0 && k ~= 1
                    TempBeadData.Distance = ((TempBeadData.CenterRow-centers(l,1)).^2+...
                        (TempBeadData.CenterCol-centers(l,2)).^2).^(1/2);
                    MinD = min(TempBeadData.Distance(TempBeadData.z == k-1));
                    if MinD < radii(l)
                        index = find((TempBeadData.Distance == min(TempBeadData.Distance)));
                        BeadID = TempBeadData.BeadID(index);
                    else
                        BeadID = MaxBeadID;
                        MaxBeadID = MaxBeadID + 1;
                    end
                    TempBeadData.Distance = [];
                    TempBeadData(rn,:) = {beadintensity, power,setID,rep, k, BeadID,...
                        centers(l, 1),centers(l, 2), radii(l), intensity};
                else
                    BeadID = MaxBeadID;
                    MaxBeadID = MaxBeadID + 1;
                    TempBeadData(rn,:) = {beadintensity, power,setID,rep, k,BeadID,...
                        centers(l, 1),centers(l, 2), radii(l), intensity};
                end
                rn = rn + 1;
            end
        end
        if mod(k, 21) == 0
            for bid =1:max(TempBeadData.BeadID)
                Is = TempBeadData.Intensity(TempBeadData.BeadID == bid);
                Zcount = length(Is);
                if Zcount >= 12

                    MeanI = mean(maxk(Is, 5));
                    MeanRow = mean(TempBeadData.CenterRow(TempBeadData.BeadID == bid));
                    MeanCol = mean(TempBeadData.CenterCol(TempBeadData.BeadID == bid));
                    MeanRad = mean(TempBeadData.Radius(TempBeadData.BeadID == bid));

                    TempMeanBeadData(height(TempMeanBeadData)+1,:) = {beadintensity,...
                        power,setID,rep, Zcount,bid,...
                        MeanRow,MeanCol, MeanRad, MeanI};
                end
            end
            % use middle 50% of beads to calculate intensity 
            testset = TempMeanBeadData.Intensity;
            CondMean = mean(TempMeanBeadData.Intensity);
            CondStd = std(TempMeanBeadData.Intensity);
            CondCount = height(TempMeanBeadData);
            SummaryData(height(SummaryData)+1,:) = {beadintensity,power,...
                setID, rep, CondMean,CondStd,CondCount};
            SummaryData(height(SummaryData),:)
            MeanBeadData = [MeanBeadData;TempMeanBeadData];

            BeadData = [BeadData; TempBeadData]; 
            MaxBeadID = 1;
            rep = rep + 1;
            rn = 1;
            TempBeadData = table('Size', [0 nvars],...
                'VariableTypes', vartypes,'VariableNames', varnames);
        end
    end

    if mod(k, 21) ~= 0
        for bid =1:max(TempBeadData.BeadID)
            Is = TempBeadData.Intensity(TempBeadData.BeadID == bid);
            Zcount = length(Is);
            if Zcount >= 12
                MeanI = mean(maxk(Is, 5));
                MeanRow = mean(TempBeadData.CenterRow(TempBeadData.BeadID == bid));
                MeanCol = mean(TempBeadData.CenterCol(TempBeadData.BeadID == bid));
                MeanRad = mean(TempBeadData.Radius(TempBeadData.BeadID == bid));
                TempMeanBeadData(height(TempMeanBeadData)+1,:) = {beadintensity,...
                    power,setID,rep, Zcount,bid,...
                    MeanRow,MeanCol, MeanRad, MeanI};
            end
        end
        CondMean = mean(TempMeanBeadData.Intensity);
        CondVar = var(TempMeanBeadData.Intensity);
        CondCount = height(TempMeanBeadData);
        SummaryData(height(SummaryData)+1,:) = {beadintensity,power,...
            setID, rep, CondMean,CondVar,CondCount};
        SummaryData(height(SummaryData),:)
        BeadData = [BeadData; TempBeadData]; 
        MeanBeadData = [MeanBeadData;TempMeanBeadData];
        rep = rep +1;
    end


end
end
%% 


outfile = 'E:\Gabriella\LivemRNA\Data\DynamicsResults\2019-12-09-GreenBeadIntensityCalibrationMeasurements\AllBeadIntensityValues.mat';
save(outfile, 'MeanBeadData', 'SummaryData', 'BeadData');

%% 
%scatter(SummaryData.LaserPower, SummaryData.MeanIntensity, 10, SummaryData.BeadIntensity)
figure(1)
SubData = SummaryData(SummaryData.BeadIntensity == .3, :);
scatter(SubData.BeadIntensity.*SubData.LaserPower,...
    SubData.MeanIntensity, 50, 'filled')
hold on
SubData = SummaryData(SummaryData.BeadIntensity == 1, :);
scatter(SubData.BeadIntensity.*SubData.LaserPower,...
    SubData.MeanIntensity, 50,  'filled')
SubData = SummaryData(SummaryData.BeadIntensity == 3, :);
scatter(SubData.BeadIntensity.*SubData.LaserPower,...
    SubData.MeanIntensity, 50,  'filled')
SubData = SummaryData(SummaryData.BeadIntensity == 10, :);
scatter(SubData.BeadIntensity.*SubData.LaserPower,...
    SubData.MeanIntensity, 50,  'filled')
xlabel('Bead Intensity * Laser power')
ylabel('Measured Intensity')
legend('.3', '1', '3', '10')
hold off

%scatter(SummaryData.LaserPower, SummaryData.MeanIntensity, 10, SummaryData.BeadIntensity)
figure(2)
SubData = SummaryData(SummaryData.BeadIntensity == .3, :);
scatter(SubData.LaserPower,...
    SubData.MeanIntensity, 50, 'filled')
hold on
SubData = SummaryData(SummaryData.BeadIntensity == 1, :);
scatter(SubData.LaserPower,...
    SubData.MeanIntensity, 50,  'filled')
SubData = SummaryData(SummaryData.BeadIntensity == 3, :);
scatter(SubData.LaserPower,...
    SubData.MeanIntensity, 50,  'filled')
SubData = SummaryData(SummaryData.BeadIntensity == 10, :);
scatter(SubData.LaserPower,...
    SubData.MeanIntensity, 50,  'filled')
xlabel('Laser power')
ylabel('Measured Intensity')
legend('.3', '1', '3', '10')
hold off





