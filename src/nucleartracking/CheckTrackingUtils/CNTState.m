% CNTState.m
% author: Gabriella Martini
% date created: 9/7/20
% date last modified: 9/7/20
classdef CNTState < handle
    
    properties
        liveExperiment
        
        
        schnitzcells
        FrameInfo
        ImageMat
        MaxImageMat
        MedImageMat
        storedTimeProjection
        multiImage
        maxTimeCell
        
        Ellipses
        nucleiModified
        
        Frames
        CurrentFrame {mustBeEmptyOrScalar(CurrentFrame)}
        PreviousFrame {mustBeEmptyOrScalar(PreviousFrame)}
        
        CurrentX {mustBeEmptyOrScalar(CurrentX)}
        CurrentZ {mustBeEmptyOrScalar(CurrentZ)}
        ManualZFlag
        ZSlices
        MaxFluo
        MedFluo
        MidMedFluo
        MaxZ
        MedZ
        MidMedZ
        
        %CurrentZ {mustBeEmptyOrScalar(CurrentZ)}
        %CurrentZIndex {mustBeEmptyOrScalar(CurrentZIndex)}
        
        CurrentNucleusCellNo {mustBeEmptyOrScalar(CurrentNucleusCellNo)}
        CurrentNucleus {mustBeEmptyOrScalar(CurrentNucleus)}
        PreviousNucleus {mustBeEmptyOrScalar(PreviousNucleus)}
        lastNucleus {mustBeEmptyOrScalar(lastNucleus)}
       
        
        CurrentChannel {mustBeEmptyOrScalar(CurrentChannel)}
        CurrentChannelIndex {mustBeEmptyOrScalar(CurrentChannelIndex)}
        PreviousChannel {mustBeEmptyOrScalar(PreviousChannel)}
        PreviousChannelIndex {mustBeEmptyOrScalar(PreviousChannelIndex)}
        inputChannel % replaces coatChannel in CPTState
        



        
        DisplayRange
        DisplayRangeSpot
        UseHistoneOverlay
        ImageHis
        HideApprovedFlag
        
        nameSuffix
        
        nWorkers
        
        projectionMode
    end
    
    methods
        function this = CNTState(liveExperiment, schnitzcells, Ellipses,...
                FrameInfo, UseHistoneOverlay, nWorkers, projectionMode)
            
            this.liveExperiment = liveExperiment;
            this.schnitzcells = schnitzcells;
            this.FrameInfo = FrameInfo;
            this.ImageMat = [];
            this.storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called
            this.multiImage = {};
            this.maxTimeCell = [];
            
            this.Ellipses = Ellipses;
            this.nucleiModified = false;
            
            this.Frames = [];
            this.CurrentFrame = 0;
            this.PreviousFrame = this.CurrentFrame;
            
            this.ManualZFlag = 0;
            numberZSlices = this.FrameInfo(1).NumberSlices;
            this.ZSlices = numberZSlices + 2; %Note that the blank slices are included
            
% First, get the different intensity values corresponding to this particle.
            this.CurrentX = [];
            this.CurrentNucleus = 1;
            this.PreviousNucleus = 1;
            this.CurrentNucleusCellNo = 1;
            this.lastNucleus = 0;
            
            Frame = [];
            MaxFluo = [];
            MaxZ = [];
            MedFluo = [];
            MedZ = [];
            MidMedFluo = [];
            MidMedZ = [];
            for i=1:length(schnitzcells(this.CurrentNucleus).frames)

                Frame(i)=schnitzcells(this.CurrentNucleus).frames(i);
                MaxFluo(i) = max(schnitzcells(this.CurrentNucleus).Fluo(i,:));
                MaxZ(i) = find(schnitzcells(this.CurrentNucleus).Fluo(i,:) == max(schnitzcells(this.CurrentNucleus).Fluo(i,:)), 1);
                MedFluo(i) = median(schnitzcells(this.CurrentNucleus).Fluo(i,2:this.ZSlices-1));
                MedZ(i) = find(schnitzcells(this.CurrentNucleus).Fluo(i,:) == median(schnitzcells(this.CurrentNucleus).Fluo(i,:)), 1);
                MidMedFluo(i) = median(schnitzcells(this.CurrentNucleus).Fluo(i,max(2, MaxZ(i)-5):min(this.ZSlices-1, MaxZ(i)+5)));
                if ~isempty(find(schnitzcells(this.CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1))
                    MidMedZ(i) = find(schnitzcells(this.CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1);
                else
                    SubFluos = schnitzcells(this.CurrentNucleus).Fluo(i,max(2, MaxZ(i)-5):min(this.ZSlices-1, MaxZ(i)+5));
                    SubFluos = sort(SubFluos(SubFluos > MidMedFluo(i)));
                    MidMedFluo(i) = SubFluos(1);
                    MidMedZ(i) = find(schnitzcells(this.CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1);

                end

            end
            this.CurrentZ = MidMedZ(find(Frame == this.CurrentFrame, 1));
            this.Frames = Frame;
            this.MaxFluo = MaxFluo;
            this.MedFluo = MedFluo;
            this.MidMedFluo = MidMedFluo;
            this.MaxZ = MaxZ;
            this.MedZ = MedZ;
            this.MidMedZ = MidMedZ;
            %this.CurrentZ = round(this.ZSlices / 2);
            
            this.CurrentNucleus = 1;
            this.PreviousNucleus = 1;
            this.CurrentNucleusCellNo = 1;
            this.lastNucleus = 0; %this gets flagged if there's a drop to one particle within the Particles structure.
            
            this.CurrentChannel = 1;
            this.CurrentChannelIndex = 1;
            this.PreviousChannel = this.CurrentChannel;
            this.PreviousChannelIndex = this.CurrentChannelIndex;
            

          
            
            this.DisplayRange = [];
            this.DisplayRangeSpot = [];
            this.UseHistoneOverlay = UseHistoneOverlay;
            this.HideApprovedFlag = 0;
            
            this.nameSuffix = '';
            
            this.nWorkers = nWorkers;
            
            this.projectionMode = projectionMode;
        end
        
        function numNuclei = numNuclei(this)
            numNuclei = length(this.schnitzcells);
        end
        
        function numValidFrames = numValidFrames(this)
            numValidFrames = length(this.Ellipses);
        end
        
        % 9/7/20 (GM): not sure what the nucleus equivalent would be here 
        % I don't think I need this functionality without having multiple
        % channels
%         function currentEllipses = getCurrentEllipses(this)
%             currentEllipses = this.Ellipses;
%         end
        
        function currentFrameEllipses = getCurrentFrameEllipses(this)
            currentEllipses = this.Ellipses();
            currentFrameEllipses = currentEllipses{this.CurrentFrame};
        end
        % I don't think I need this functionality without having multiple
        % channels
%         function currentSchnitzCells = getCurrentChannelParticles(this)
%             currentParticles = this.Particles{this.CurrentChannelIndex};
%         end
        
        function currentSchnitzCell = getCurrentSchnitzCell(this)
            currentSchnitzCells = this.schnitzcells;
            currentSchnitzCell = currentSchnitzCells(this.CurrentNucleus);
        end
        
%         function currentParticleFit = getCurrentParticleFit(this)
%             currentFrameSpots = this.getCurrentFrameSpots();
%             currentParticleFit = currentFrameSpots.Fits(this.CurrentParticleIndex);
%         end

        
%         function currentXDoG = getCurrentXDoG(this)
%             currentFit = this.getCurrentParticleFit();
%             currentXDoG = double(currentFit.xDoG(this.CurrentZIndex));
%         end

          function currentNucleusX = getCurrentX(this)
              currentFrameIndex = find(this.schnitzcells(this.CurrentNucleus).frames == this.CurrentFrame, 1);
              currentNucleusX = this.schnitzcells(this.CurrentNucleus).cenx(currentFrameIndex);
          end
%         
%         function currentYDoG = getCurrentYDoG(this)
%             currentFit = this.getCurrentParticleFit();
%             currentYDoG = double(currentFit.yDoG(this.CurrentZIndex));
%         end
          function currentNucleusY = getCurrentY(this)
              currentFrameIndex = find(this.schnitzcells(this.CurrentNucleus).frames == this.CurrentFrame, 1);
              currentNucleusY = this.schnitzcells(this.CurrentNucleus).ceny(currentFrameIndex);
          end
        
%         function currentXFit = getCurrentXFit(this)
%             currentFit = this.getCurrentParticleFit();
%             currentXFit = double(currentFit.xFit(this.CurrentZIndex));
%         end
%         
%         function currentYFit = getCurrentYFit(this)
%             currentFit = this.getCurrentParticleFit();
%             currentYFit = double(currentFit.yFit(this.CurrentZIndex));
%         end
        
        % WHAT SHOULD THIS BE?
%         function updateCurrentZIndex(this)
%             this.CurrentZIndex = this.getMaxZIndex();
%         end
        
%         function maxZIndex = getMaxZIndex(this)
%             currentFrameIndex = find(this.getCurrentSchnitzCell().frames == this.CurrentFrame);
%             maxZIndex = find(this.getCurrentSchnitzCell().Fluo(currentFrameIndex,:) == max(this.getCurrentSchnitzCell().Fluo(currentFrameIndex,:)), 1);
%         end
%         
        function currentNucleusCellNo = getCurrentNucleusCellNo(this)
            currentNucleusCellNo = this.getCurrentSchnitzCell().cellno(...
                this.getCurrentSchnitzCell().frames == this.CurrentFrame);
        end
        
        function [xApproved, yApproved] = getApprovedSchnitzCells(this, x, y)
            IndexApprovedSchnitzCells = [];
            numNuclei = this.numNuclei();
            currentSchnitzCells = this.schnitzcells;
            
            for i = 1:numNuclei
                if sum(currentSchnitzCells(i).frames == this.CurrentFrame) &&...
                        sum(currentSchnitzCells(i).Approved == 1)
                    IndexApprovedSchnitzCells = [IndexApprovedSchnitzCells,...
                        currentSchnitzCells(i).cellno(currentSchnitzCells(i).frames == this.CurrentFrame)];
                end
            end
            
            xApproved = x(IndexApprovedSchnitzCells);
            yApproved = y(IndexApprovedSchnitzCells);
        end
        
        function [xDisapproved, yDisapproved] = getDisapprovedSchnitzCells(this, x, y)
            IndexDisapprovedSchnitzCells=[];
            numNuclei = this.numNuclei();
            currentSchnitzCells = this.schnitzcells;
            
            for i = 1:numNuclei
                if sum(currentSchnitzCells(i).frames == this.CurrentFrame) &&...
                        sum(currentSchnitzCells(i).Approved == -1)
                    IndexDisapprovedSchnitzCells=[IndexDisapprovedSchnitzCells,...
                        currentSchnitzCells(i).cellno(currentSchnitzCells(i).frames == this.CurrentFrame)];
                end
            end
            
            xDisapproved = x(IndexDisapprovedSchnitzCells);
            yDisapproved = y(IndexDisapprovedSchnitzCells);
        end
        
        function [xNonFlagged, yNonFlagged] = getNonFlaggedSchnitzCells(this, x, y)
            IndexNonFlaggedSchnitzCells = [];
            numNuclei = this.numNuclei();
            currentSchnitzCells = this.schnitzcells;
            
            
            for i = 1:numNuclei
                if sum(currentSchnitzCells(i).frames == this.CurrentFrame) &&...
                        ~(sum(currentSchnitzCells(i).Approved == -1) || sum(currentSchnitzCells(i).Approved == 1))
                    IndexNonFlaggedSchnitzCells=[IndexNonFlaggedSchnitzCells,...
                        currentSchnitzCells(i).cellno(currentSchnitzCells(i).frames == this.CurrentFrame)];
                end
            end
            
            xNonFlagged = x(IndexNonFlaggedSchnitzCells);
            yNonFlagged = y(IndexNonFlaggedSchnitzCells);
        end
        
       
      
        function processImageMatrices(this, movieMat)
            
            fr_idx = find(this.Frames == this.CurrentFrame);
            maxz = this.MaxZ(fr_idx);
            medz = this.MedZ(fr_idx);
            midmedz = this.MidMedZ(fr_idx);
            

            if ~isempty(movieMat)
                if ~isempty(midmedz)
                    this.ImageMat = movieMat(:, :, midmedz,...
                        this.CurrentFrame, this.CurrentChannel);
                else
                    this.ImageMat = zeros(size(movieMat, 1), size(movieMat, 2), 'uint8');
                end
                if ~isempty(maxz)
                    this.MaxImageMat = movieMat(:, :, maxz,...
                        this.CurrentFrame, this.CurrentChannel);
                else
                    this.MaxImageMat = zeros(size(movieMat, 1), size(movieMat, 2), 'uint8');
                end
                if ~isempty(medz)
                    this.MedImageMat = movieMat(:, :, medz,...
                        this.CurrentFrame, this.CurrentChannel);
                else
                    this.MedImageMat = zeros(size(movieMat, 1), size(movieMat, 2), 'uint8');
                end
            else
                if ~isempty(midmedz)
                    this.ImageMat = getMovieSlice(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel, midmedz);
                else
                    dummyMat = getMovieSlice(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel, 1);
                    this.ImageMat = zeros(size(dummyMat, 1), size(dummyMat, 2), 'uint8');
                end
                if ~isempty(maxz)
                    this.MaxImageMat = getMovieSlice(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel, maxz);
                else
                    dummyMat = getMovieSlice(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel, 1);
                    this.MaxImageMat = zeros(size(dummyMat, 1), size(dummyMat, 2), 'uint8');
                end
                if ~isempty(medz)
                    this.MedImageMat = getMovieSlice(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel, medz);
                else
                    dummyMat = getMovieSlice(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel, 1);
                    this.MedImageMat = zeros(size(dummyMat, 1), size(dummyMat, 2), 'uint8');
                end
            end
            
        disp('stop here')
        end
        % NOT SURE WHAT THESE & and y INPUTS ARE
        function [xTrace, yTrace] = getXYTraces(this, x, y)
            xTrace = x(this.CurrentNucleusCellNo);
            yTrace = y(this.CurrentNucleusCellNo);
        end
        
        % NOT SURE WHAT THESE & and y INPUTS ARE
        function updateTraceInfo(this)
            Frame = [];
            MaxFluo = [];
            MaxZ = [];
            MedFluo = [];
            MedZ = [];
            MidMedFluo = [];
            MidMedZ = [];
            for i=1:length(this.schnitzcells(this.CurrentNucleus).frames)

                Frame(i)=this.schnitzcells(this.CurrentNucleus).frames(i);
                MaxFluo(i) = max(this.schnitzcells(this.CurrentNucleus).Fluo(i,:));
                MaxZ(i) = find(this.schnitzcells(this.CurrentNucleus).Fluo(i,:) == max(this.schnitzcells(this.CurrentNucleus).Fluo(i,:)), 1);
                MedFluo(i) = median(this.schnitzcells(this.CurrentNucleus).Fluo(i,2:this.ZSlices-1));
                MedZ(i) = find(this.schnitzcells(this.CurrentNucleus).Fluo(i,:) == median(this.schnitzcells(this.CurrentNucleus).Fluo(i,:)), 1);
                MidMedFluo(i) = median(this.schnitzcells(this.CurrentNucleus).Fluo(i,max(2, MaxZ(i)-5):min(this.ZSlices-1, MaxZ(i)+5)));
                if ~isempty(find(this.schnitzcells(this.CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1))
                    MidMedZ(i) = find(this.schnitzcells(this.CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1);
                else
                    SubFluos = this.schnitzcells(this.CurrentNucleus).Fluo(i,max(2, MaxZ(i)-5):min(this.ZSlices-1, MaxZ(i)+5));
                    SubFluos = sort(SubFluos(SubFluos > MidMedFluo(i)));
                    MidMedFluo(i) = SubFluos(1);
                    MidMedZ(i) = find(this.schnitzcells(this.CurrentNucleus).Fluo(i,:) == MidMedFluo(i), 1);

                end

            end
            this.CurrentZ = MidMedZ(find(Frame == this.CurrentFrame, 1));
            this.Frames = Frame;
            this.MaxFluo = MaxFluo;
            this.MedFluo = MedFluo;
            this.MidMedFluo = MidMedFluo;
            this.MaxZ = MaxZ;
            this.MedZ = MedZ;
            this.MidMedZ = MidMedZ;
            
        end
        
        function updateCurrentNucleusCellNo(this)
            this.CurrentNucleusCellNo = this.getCurrentNucleusCellNo();
        end
    end
end

function mustBeEmptyOrScalar(in)

if numel(in) > 1
    error(['Value assigned to property is not '...
        'scalar or empty']);
end

end