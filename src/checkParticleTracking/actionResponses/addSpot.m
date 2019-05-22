function [numParticles, SpotFilter, Particles, Spots, PreviousParticle] =...
    addSpot(ZoomMode, GlobalZoomMode, Particles, CurrentChannel, ...
    CurrentParticle, CurrentFrame, CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
    LinesPerFrame, Spots, ZSlices, PathPart1, PathPart2, Path3, FrameInfo, pixelSize, ...
    SpotFilter, numParticles, cc, xSize, ySize, NDigits, intScale, Prefix, PreProcPath, ProcPath)
%ADDSPOT Summary of this function goes here
%   Detailed explanation goes here

zStep = FrameInfo(1).ZStep;
saveType = '.tif';

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

        [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'r', 'linewidth',1, 'FigHandle', Overlay);
        if ConnectPositionx < 1 || ConnectPositiony < 1
            %sometimes ginputc returns the wrong coordinates for an
            %unknown reason. if that happens, we'll resort to a
            %black crosshair from ginput.
            [ConnectPositionx,ConnectPositiony] = ginput(1);
        end

        ConnectPositionx = round(ConnectPositionx);
        ConnectPositiony = round(ConnectPositiony);

        % check that the clicked particle isn't too close to the
        % edge of the frame
        if (ConnectPositionx > snippet_size/2) && (ConnectPositionx + snippet_size/2 < PixelsPerLine)...
                && (ConnectPositiony > snippet_size/2) && (ConnectPositiony + snippet_size/2 < LinesPerFrame)
            SpotsIndex = length(Spots{CurrentChannel}(CurrentFrame).Fits)+1;
            breakflag = 0; %this catches when the spot addition was unsuccessful and allows checkparticletracking to keep running and not error out
            maxWorkers = 8;
            use_integral_center = 1;
            
%             startParallelPool(maxWorkers, false, true);
            Fits = [];
            for z = 1:ZSlices %#ok<PFUIX>
                imAbove = [];
                imBelow = [];
                spotsIm = [];
                spotsIm=imread([PathPart1,iIndex(CurrentFrame,NDigits),...
                    '_z',iIndex(z,2),PathPart2]);
                try
                    imAbove = double(imread([PathPart1,iIndex(CurrentFrame,NDigits),...
                        '_z',iIndex(z-1,2),PathPart2]));
                    imBelow = double(imread([PathPart1,iIndex(CurrentFrame,NDigits),...
                        '_z',iIndex(z+1,2),PathPart2]));
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
                k = 1; %This is supposed to be the index for the partiles in an image.
                %However, this image only contains one particle
                neighborhood = round(1300 / pixelSize); %nm
                %Get the information about the spot on this z-slice
                if cc == '['
                    [temp_particles{z}, Fit] = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, im_label, dog, neighborhood, snippet_size, ...
                        pixelSize, show_status, fig, microscope, [1, ConnectPositionx, ConnectPositiony], [], '', intScale, CurrentFrame, [], z, use_integral_center);
                elseif cc == '{'
                    [temp_particles{z}, Fit] = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, im_label, dog, neighborhood, snippet_size, ...
                        pixelSize, show_status, fig, microscope, [1, ConnectPositionx, ConnectPositiony], [ConnectPositionx, ConnectPositiony], '', intScale, CurrentFrame, [], z, use_integral_center);
                end
                
                if ~isempty(Fit)
                    fieldnames = fields(Fit);
                    if isempty(Fits)
                        Fits = Fit;
                    else
                        for i = 1:length(fieldnames)
                            Fits.(fieldnames{i}) =  [Fits.(fieldnames{i}), Fit.(fieldnames{i})];
                        end
                    end
                end
            end
            
%             if ~isempty(Fits)
%                 if ~isempty(Spots{CurrentChannel}(CurrentFrame).Fits)
%                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex) = Fits;
%                 else
%                      Spots{CurrentChannel}(CurrentFrame).Fits = Fits;
%                 end
%             else
%                 disp('No spot added. Did you click too close to the image boundary?');
%                 %                     breakflag = 1;
%                 %                     break
%             end
%             Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex) = Fits;

% %             for z = 1:zSlices
%                 if ~isempty(Fits)
%                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex) = Fits;
%                 else
%                     disp('No spot added. Did you click too close to the image boundary?');
% %                     breakflag = 1;
% %                     break
%                 end
%             end
            
%             for zIndex = 1:ZSlices
%                 if ~isempty(temp_particles{zIndex})
%                     %Copy the information stored in temp_particles into the
%                     %Spots structure
%                     if ~isempty(temp_particles{zIndex}{1})
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity(zIndex)=...
%                             temp_particles{zIndex}{1}{1};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xFit(zIndex)=...
%                             temp_particles{zIndex}{1}{2};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yFit(zIndex)=...
%                             temp_particles{zIndex}{1}{3};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Offset(zIndex)=...
%                             temp_particles{zIndex}{1}{4};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Area(zIndex)=...
%                             temp_particles{zIndex}{1}{6};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yDoG(zIndex)=...
%                             temp_particles{zIndex}{1}{9};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xDoG(zIndex)=...
%                             temp_particles{zIndex}{1}{10};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).GaussianIntensity(zIndex)=...
%                             temp_particles{zIndex}{1}{11};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).CentralIntensity(zIndex)=...
%                             temp_particles{zIndex}{1}{12};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).DOGIntensity(zIndex)=...
%                             temp_particles{zIndex}{1}{13};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).ConfidenceIntervals{zIndex}=...
%                             temp_particles{zIndex}{1}{19};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).z(zIndex)=...
%                             zIndex;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).gaussParams{zIndex}=...
%                             temp_particles{zIndex}{1}{22};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).intArea=...
%                             temp_particles{zIndex}{1}{23};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).discardThis=...
%                             0;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).frame=...
%                             CurrentFrame;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).r = 0;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity3 = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity5 = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).cylIntensity(zIndex) = temp_particles{zIndex}{1}{24};
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).brightestZ = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).IntegralZ = use_integral_center;
% 
%                     else
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity(zIndex)=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xFit(zIndex)=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yFit(zIndex)=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Offset(zIndex)=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Area(zIndex)= nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yDoG(zIndex)= nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xDoG(zIndex)= nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).GaussianIntensity(zIndex)=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).CentralIntensity(zIndex)=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).DOGIntensity(zIndex)=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).ConfidenceIntervals{zIndex}=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).z(zIndex)=zIndex;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).gaussParams{zIndex}=nan;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).discardThis=0;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).frame=CurrentFrame;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).r=0;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity3 = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).cylIntensity(zIndex) = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity5 = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).brightestZ = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).IntegralZ = use_integral_center;
%                     end
%                 else
%                     disp('No spot added. Did you click too close to the image boundary?')
%                     breakflag = 1;
%                     break
%                 end
%             end

%             allNaNs = sum(isnan(Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity)) ==...
%                 length(Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity);
%             if allNaNs
%                 Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex) = [];
%                 if isempty(Spots{CurrentChannel}(CurrentFrame).Fits)
%                     Spots{CurrentChannel}(CurrentFrame).Fits = [];
%                 end
%                 breakflag = 1;
%             end

            if isempty(Fits)
                breakflag = true;
            end
            
            if ~breakflag
                if cc == '['
                    force_z = 0;
                elseif cc == '{'
                    force_z = CurrentZ;
                end
%                 [Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex),~, ~] = findBrightestZ(Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex), -1, use_integral_center, force_z, []);
%                   [a,~, ~] = findBrightestZ(Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex), -1, use_integral_center, force_z, []);
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
                   Spots{CurrentChannel} =...
                   ...
                   fitSnip3D(...
                   ...
                   Spots{CurrentChannel}, CurrentChannel, SpotsIndex, CurrentFrame,...
                   Prefix, PreProcPath, ProcPath, FrameInfo, [], false, saveType); 
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
                [SpotFilter{CurrentChannel},Particles{CurrentChannel}]=...
                    TransferParticle(Spots{CurrentChannel},...
                    SpotFilter{CurrentChannel},Particles{CurrentChannel},...
                    CurrentFrame,SpotsIndex);
                numParticles = numParticles + 1;

                %Connect this particle to the CurrentParticle. This is
                %the equivalent of running the 'c' command.
                if ~GlobalZoomMode
                    Particles{CurrentChannel}=...
                        JoinParticleTraces(CurrentParticle,...
                        numParticles,Particles{CurrentChannel});
                else
                    disp('Re-run TrackmRNADynamics to associate this particle with a nucleus.')
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

