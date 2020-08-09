function [ParticleStitchInfo, SpotFilter, Particles, Spots,...
    PreviousParticle, CurrentParticle, ZoomMode, GlobalZoomMode] =...
    ...
    addSpot(...
    ...
    ZoomMode, GlobalZoomMode, Particles, ParticleStitchInfo, CurrentChannel, ...
    CurrentParticle, CurrentFrame, CurrentZ, Spots,...
    SpotFilter, cc, Prefix,UseHistoneOverlay,...
    schnitzcells, nWorkers, plot3DGauss, imStack)

%ADDSPOT


liveExperiment = LiveExperiment(Prefix);
globalMotionModel = getGlobalMotionModel(liveExperiment);
% 
% movieMat = getMovieMat(thisExperiment);
% imStack = movieMat(:, :, :, CurrentFrame, CurrentChannel);

FrameInfo = getFrameInfo(liveExperiment);
LinesPerFrame = liveExperiment.yDim;
PixelsPerLine = liveExperiment.xDim;
pixelSize_nm = liveExperiment.pixelSize_nm;
PreProcPath = liveExperiment.preFolder;
snippetSize_px = liveExperiment.snippetSize_px;
nSlices = liveExperiment.zDim;


numParticles = length(Particles{CurrentChannel});
startParallelPool(nWorkers, 0, 1);

%Check that we're in zoom mode. If not, set it up.
PreviousParticle = 0; % resets particle so trace will refresh
if ~(ZoomMode || GlobalZoomMode)
    disp('You need to be in Zoom Mode to do this. You can switch using ''o'' or ''+''. Run the ''['' command again.')
else
    %Click on the region we're going to fit in order to find a new
    %spot and add that spot to the current particle
    
    %Check that this particle doesn't already have a spot assigned
    %in this frame
    if sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame) &&  ~GlobalZoomMode
        warning('There is a spot assigned to this particle in this frame already.')
    else
        
        
        [ConnectPositionx,ConnectPositiony] = ginput(1);
        
        ConnectPositionx = round(ConnectPositionx);
        ConnectPositiony = round(ConnectPositiony);
        
        % check that the clicked particle isn't too close to the
        % edge of the frame
        if (ConnectPositionx > snippetSize_px/2) && (ConnectPositionx + snippetSize_px/2 < PixelsPerLine)...
                && (ConnectPositiony > snippetSize_px/2) && (ConnectPositiony + snippetSize_px/2 < LinesPerFrame)
            SpotsIndex = length(Spots{CurrentChannel}(CurrentFrame).Fits)+1;
            breakflag = 0; %this catches when the spot addition was unsuccessful and allows checkparticletracking to keep running and not error out
            use_integral_center = 1;
            
            FitCell = cell(1, nSlices);
            
            for z = 1:nSlices
                spotsIm = imStack(:, :, z);
                try
                    imAbove= imStack(:, :, z+1);
                   imBelow= imStack(:, :, z-1);
                catch
                    imAbove = nan(size(spotsIm,1),size(spotsIm,2));
                    imBelow = nan(size(spotsIm,1),size(spotsIm,2));
                end
                
                Threshold = min(min(spotsIm));
                dog = spotsIm;
                im_thresh = dog >= Threshold;
                [im_label, ~] = bwlabel(im_thresh);
                microscope = FrameInfo(1).FileMode;
                show_status = 0;
                fig = [];
                k = 1; %This is supposed to be the index for the particles in an image.
                %However, this image only contains one particle
                neighborhood_px = round(1300 / pixelSize_nm); %nm
                %Get the information about the spot on this z-slice
                if cc == '['
                    Fit = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, im_label, dog, neighborhood_px, snippetSize_px, ...
                        pixelSize_nm, show_status, fig, microscope,...
                        [1, ConnectPositionx, ConnectPositiony], [], '', CurrentFrame, [], z);
                elseif cc == '{'
                     Fit = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, im_label, dog, neighborhood_px, snippetSize_px, ...
                        pixelSize_nm, show_status, fig, microscope,...
                        [1, ConnectPositionx, ConnectPositiony], [ConnectPositionx, ConnectPositiony], '', CurrentFrame, [], z);
                end
                
                FitCell{z} = Fit;
                
                
            end
            Fits = [];
            
           for z = 1:nSlices
                if ~isempty(FitCell{z})
                    fieldnames = fields(FitCell{z});
                    if isempty(Fits)
                        Fits = FitCell{z};
                    else
                        for field = 1:length(fieldnames)
                            Fits.(fieldnames{field}) =  [Fits.(fieldnames{field}), FitCell{z}.(fieldnames{field})];
                        end
                    end
                end
            end
            
            if isempty(Fits)
                breakflag = true;
            end
            
            if ~breakflag
                if cc == '['
                    force_z = 0;
                elseif cc == '{'
                    force_z = CurrentZ;
                end
                [Fits,~, ~] = findBrightestZ(Fits, -1, use_integral_center, force_z, []);
                
                
                
                if SpotsIndex ~= 1
                    if ~isempty(setdiff(fields(Spots{CurrentChannel}(CurrentFrame).Fits), fields(Fits)))...
                            | ~isempty(setdiff(fields(Fits),fields(Spots{CurrentChannel}(CurrentFrame).Fits)))
                        [Spots{CurrentChannel}(CurrentFrame).Fits, Fits] = addFields(Spots{CurrentChannel}(CurrentFrame).Fits, Fits);
                    end
                    Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex) = Fits;
                else
                    Spots{CurrentChannel}(CurrentFrame).Fits = Fits;
                end
                %%
                if plot3DGauss
                    nSpots = 1;
                    Spots{CurrentChannel}(CurrentFrame) =...
                        ...
                        fitSnip3D(...
                        ...
                        Spots{CurrentChannel}(CurrentFrame), CurrentChannel, SpotsIndex, CurrentFrame,...
                        Prefix, PreProcPath, FrameInfo, nSpots, imStack);
                end
                %%
                %Add this to SpotFilter, which tells the code that this spot is
                %above the threshold. First, check whether the
                %dimensions of SpotFilter need to be altered. If so, pad it with NaNs
                if size(SpotFilter{CurrentChannel},2)>SpotsIndex
                    SpotFilter{CurrentChannel}(CurrentFrame,SpotsIndex)=1;
                else
                    %Pad with NaNs
                    SpotFilter{CurrentChannel}(:,end:SpotsIndex)=NaN;
                    SpotFilter{CurrentChannel}(CurrentFrame,SpotsIndex)=1;
                end
                
                %Turn this spot into a new particle. This is the equivalent of
                %the 'u' command.
%                 try
                    [SpotFilter{CurrentChannel},Particles{CurrentChannel}]=...
                        TransferParticle(Spots{CurrentChannel},...
                        SpotFilter{CurrentChannel},Particles{CurrentChannel},...
                        CurrentFrame,SpotsIndex,globalMotionModel);
%                 catch
%                     warning('failed to add spot for unknown reason.')
%                     return;
%                 end
                numParticles = numParticles + 1;
                
                %Connect this particle to the CurrentParticle. This is
                %the equivalent of running the 'c' command.
                if ~GlobalZoomMode
                    [Particles{CurrentChannel}, ParticleStitchInfo{CurrentChannel}]=...
                        JoinParticleTraces(CurrentParticle,...
                        numParticles,Particles{CurrentChannel},ParticleStitchInfo{CurrentChannel});
                else %NL: is this "else" option ever realized?
                    CurrentParticle = length(Particles{CurrentChannel});
                    Particles = addNucleusToParticle(Particles, CurrentFrame, ...
                        CurrentChannel, UseHistoneOverlay, schnitzcells, CurrentParticle);
                    GlobalZoomMode = false;
                    ZoomMode = true;
                    disp('Creating new particle trace...');
                end                                
                
                disp('Spot addded to the current particle.')
            else
                warning('You clicked too close to the edge. A spot can''t be added here.');
                msgbox('You clicked too close to the edge. A spot can''t be added here.');
            end
        end
    end
end
clear breakflag;

end

