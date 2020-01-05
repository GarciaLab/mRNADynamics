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

for i=10:15%length(rawDataDir)
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
radmin = round(radius*1/2, 0);
radmax = round(radius, 0);
rep = 1;
for j =1:size(LIFImages, 1)
    NSlices = length(LIFImages{j,1});
    %images = zeros([size(LIFImages{1}{1, 1}, 1), size(LIFImages{1}{1, 1}, 2), NSlices]);
    rn = 1;
    TempBeadData = table('Size', [0 nvars],...
    'VariableTypes', vartypes,'VariableNames', varnames);
    MaxBeadID = 1;
    for k=1:NSlices
        %images(:,:,k) = LIFImages{j}{k,1};
        img = LIFImages{j}{k,1};
        [centers, radii] = imfindcircles(img, [radmin, radmax],...
            'Sensitivity', 0.9, 'Method','twostage');
        keepcenters = ones(length(centers), 1);
        for ic =1:(length(centers)-1)
            for jc=(ic+1):length(centers)
                d = sqrt((centers(ic, 1)-centers(jc, 1))^2+(centers(ic, 2)-centers(jc,2))^2);
                if d < radii(ic)+radii(jc)
                    keepcenters(jc) = 0;
                end
            end
        end
        centers = centers(keepcenters == 1, :);
        radii = radii(keepcenters == 1);
        if length(centers) == 0
            continue
        else
        for l=1:length(centers)
            [xgrid, ygrid] = meshgrid(1:size(img,1), 1:size(img,2));
            mask = ((xgrid-centers(l,1)).^2 + (ygrid-centers(l,2)).^2) <= (radii(l)/2).^2;
            values = img(mask);
            intensity = mean(values);
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
        if k == 21
            for bid =1:max(TempBeadData.BeadID)
                Is = TempBeadData.Intensity(TempBeadData.BeadID == bid);
                Zcount = length(Is);
                if Zcount > 10
                    MeanI = mean(maxk(Is, 5));
                    MeanRow = mean(TempBeadData.CenterRow(TempBeadData.BeadID == bid));
                    MeanCol = mean(TempBeadData.CenterCol(TempBeadData.BeadID == bid));
                    MeanRad = mean(TempBeadData.Radius(TempBeadData.BeadID == bid));
                    MeanBeadData(height(MeanBeadData)+1,:) = {beadintensity,...
                        power,setID,rep, Zcount,bid,...
                        MeanRow,MeanCol, MeanRad, MeanI};
                end
            end
            BeadData = [BeadData; TempBeadData]; 
            MaxBeadID = 1;
            rep = rep + 1;
            rn = 1;
            TempBeadData = table('Size', [0 nvars],...
                'VariableTypes', vartypes,'VariableNames', varnames);
        end
    end
    if k < 21
        rep = rep +1;
    end
    for bid =1:max(TempBeadData.BeadID)
        Is = TempBeadData.Intensity(TempBeadData.BeadID == bid);
        Zcount = length(Is);
        if Zcount > 10
            MeanI = mean(maxk(Is, 5));
            MeanRow = mean(TempBeadData.CenterRow(TempBeadData.BeadID == bid));
            MeanCol = mean(TempBeadData.CenterCol(TempBeadData.BeadID == bid));
            MeanRad = mean(TempBeadData.Radius(TempBeadData.BeadID == bid));
            MeanBeadData(height(MeanBeadData)+1,:) = {beadintensity,...
                power,setID,rep, Zcount,bid,...
                MeanRow,MeanCol, MeanRad, MeanI};
        end
    end
    BeadData = [BeadData; TempBeadData]; 
end
end
%% 




%% 


img = LIFImages{1}{8,1};

[centers, radii] = imfindcircles(img, [radmin, radmax],...
    'Sensitivity', 0.9, 'Method','twostage');

figure(1)
imshow(uint8(img))
h = viscircles(centers, radii)
[xgrid, ygrid] = meshgrid(1:size(img,1), 1:size(img,2));
mask = ((xgrid-centers(1,1)).^2 + (ygrid-centers(1,2)).^2) <= radii(1).^2;
values = img(mask);
intensity = mean(values);
%% 

img2 = images(:,:,2);

[centers2, radii2] = imfindcircles(img2, [radmin, radmax],...
    'Sensitivity', 0.95);
figure(2)
imshow(images(:,:,2))
h = viscircles(centers2, radii2)
[xgrid, ygrid] = meshgrid(1:size(img2,1), 1:size(img2,2));
mask = ((xgrid-centers2(1,1)).^2 + (ygrid-centers2(1,2)).^2) <= radii2(1).^2;
values2 = img(mask);
intensity2 = mean(values2);


