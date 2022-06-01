function ComputeFixedNucleiFluo(Prefix, varargin)


ShowPlots = false;

liveExperiment= LiveExperiment(Prefix);

CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

Ellipses = getEllipses(liveExperiment);
FrameInfo = getFrameInfo(liveExperiment);
try
    clrmp = single(hsv(20));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
    %in case the user doesn't have this colormap, just keep going.
end
nSlices=FrameInfo(1).NumberSlices;
%%
ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
Channel4 = liveExperiment.Channel4;
Channel5 = liveExperiment.Channel5;


xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;

%Get the nuclei segmentation data
Ellipses = getEllipses(liveExperiment);
FrameInfo = getFrameInfo(liveExperiment);

try
    clrmp = single(hsv(20));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
    %in case the user doesn't have this colormap, just keep going.
end
Channels = {Channel1, Channel2, Channel3, Channel4, Channel5};





InputChannelIndexes = find(contains(Channels, 'input', 'IgnoreCase', true));
HisChannelIndexes = find(contains(Channels, 'his', 'IgnoreCase', true));
ChannelsToIntegrate = unique([HisChannelIndexes InputChannelIndexes]);
if isempty(InputChannelIndexes)
    warning(['No input channel found. Check correct definition in MovieDatabase.',...
        ' Input channels should use the :input notation.'])
    return;
end
%%

close all
nEmbryos = length(FrameInfo);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);
GoodEmbryos = 1:nEmbryos;
GoodEmbryos = GoodEmbryos(CompiledEmbryos.Approved);
EllipsesFluoInfo = cell(size(Ellipses));
%%
close all
if sum(InputChannelIndexes)
    for CurrentEmbryo =1:nEmbryos
        if ismember(CurrentEmbryo, GoodEmbryos) & ~isempty(Ellipses{CurrentEmbryo})
            
            ellipseFrame = double(Ellipses{CurrentEmbryo});
            ellipseFrame = ellipseFrame(ellipseFrame(:,1) > 0, :);
            ellipseFrame = ellipseFrame(ellipseFrame(:,2) > 0, :);
            Ellipses{CurrentEmbryo} = ellipseFrame;
            IntegrationRadius=ellipseFrame(1,3)/2;       %Radius of the integration region in um
            x = 1:floor(ceil(2*IntegrationRadius)/2)*2+1;
            y = 1:floor(ceil(2*IntegrationRadius)/2)*2+1;
            [X, Y] = meshgrid(x, y);
            FluoFilter = double(sqrt((X-median(x)).^2 + (Y-median(y)).^2) < ceil(2*IntegrationRadius)/2);
            %Initialize fields
            EmbryoSchnitzFluo = NaN(size(ellipseFrame, 1),8+nSlices, length(Channels));
            
            CoordA = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,:);
            CoordP = CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,:);
            CoordD = CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,:);
            CoordV = CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,:);
            %bIndex = boundary(ellipseFrame(:,1), ellipseFrame(:,2), 1);
            APLength = CoordP(1)-CoordA(1);
            DVLength = CoordV(2)-CoordD(2);
            
            for channelIndex = 1:length(Channels)
                
                EmbryoSchnitzFluo(:,1:3, channelIndex) = ellipseFrame(:,1:3);
                EmbryoSchnitzFluo(:,4:5, channelIndex) = ellipseFrame(:,5:6);
                EmbryoSchnitzFluo(:,6, channelIndex) = (EmbryoSchnitzFluo(:,1, channelIndex)-CoordA(1))/APLength;
                EmbryoSchnitzFluo(:,7, channelIndex) = (EmbryoSchnitzFluo(:,2, channelIndex)-CoordD(2))/DVLength;
                if ismember(channelIndex, ChannelsToIntegrate)
                EmbryoChStackFile = [liveExperiment.preFolder, filesep, Prefix,...
                    '_Position',iIndex(CurrentEmbryo,3),'_001_ch', iIndex(channelIndex,2),'.tif'];
                
                
                EmbryoChStack = imreadStack2(EmbryoChStackFile, liveExperiment.yDim, liveExperiment.xDim,...
                    liveExperiment.nFrames);
                zIndex = 1;
                if ShowPlots
                    close all
                    EmbryoFig = figure(1);
                    fullAxes = axes(EmbryoFig);
                    EmbryoImage = imshow(EmbryoChStack(:,:,zIndex), [min(min(EmbryoChStack(:,:,zIndex))), ...
                        max(max(EmbryoChStack(:,:, zIndex)))], 'Border','Tight','Parent',fullAxes);
                    FilterEmbryoImage = imgaussfilt(EmbryoChStack(:,:,zIndex), 2/PixelSize_um);
                    hold(fullAxes, 'on')
                    
                    for k=1:size(ellipseFrame, 1)
                        n =k;
                        %         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
                        %             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
                        %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
                        %             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
                        colorhash = uint8(mod(round(ellipseFrame(n, 1)+ellipseFrame(n, 2)),20)+1);
                        PlotHandle{k} = ellipse(ellipseFrame(n, 3), ellipseFrame(n, 4),...
                            ellipseFrame(n, 5) * (360/(2*pi)), ellipseFrame(n, 5),...
                            ellipseFrame(n, 6), clrmp(colorhash,:), 10, fullAxes, 0.5);
                        
                    end
                    hold(fullAxes, 'off')
                end
                
                
                
                %%
                
                
                %Create the circle that we'll use as the mask
                
                
                
                
                %Extract the fluroescence of each schnitz, for each channel,
                
                
                
                %         nameSuffix=['_ch',iIndex(ChN,2)];
                
                %
                convImage = imfilter(EmbryoChStack, FluoFilter, 'same');
                
                for j=1:length(EmbryoSchnitzFluo)
                    for zIndex = 1:nSlices
                        x_coord = uint16(round(EmbryoSchnitzFluo(j, 4, channelIndex)));
                        y_coord = uint16(round(EmbryoSchnitzFluo(j, 5, channelIndex)));
                        if x_coord < 1
                            x_coord = 1;
                        end
                        if x_coord > size(convImage, 2)
                            x_coord = size(convImage, 2);
                        end
                        
                        if y_coord < 1
                            y_coord = 1;
                        end
                        if y_coord > size(convImage, 1)
                            y_coord = size(convImage, 1);
                        end
               
                        EmbryoSchnitzFluo(j, 8+zIndex, channelIndex) = convImage(y_coord, x_coord, zIndex)/sum(FluoFilter(:));
                     
                      
                    end
                    
                end %loop of nuclei in a frame
                EmbryoSchnitzFluo(:, 8, channelIndex) = max(EmbryoSchnitzFluo(:, 9:9+nSlices-1, channelIndex), [], 2);
                end
            end
            EllipsesFluoInfo{CurrentEmbryo} = EmbryoSchnitzFluo;
        else
            EllipsesFluoInfo{CurrentEmbryo} = NaN(size(0,8+nSlices, length(Channels)));
        end
    end
    
    
    
else
    
    error('Input channel not recognized. Check correct definition in MovieDatabase.Input channels should use the :input notation.');
    
end

%%
EmbryoSchnitzPath = [liveExperiment.resultsFolder, filesep, 'EllipsesFluoInfo.mat'];
save(EmbryoSchnitzPath, 'EllipsesFluoInfo')
