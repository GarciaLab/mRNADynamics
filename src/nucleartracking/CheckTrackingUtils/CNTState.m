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
        storedTimeProjection
        multiImage
        maxTimeCell
        
        Ellipses
        nucleiModified
        
        Frames
        CurrentFrame {mustBeEmptyOrScalar(CurrentFrame)}
        PreviousFrame {mustBeEmptyOrScalar(PreviousFrame)}
        
        ManualZFlag
        ZSlices
        CurrentZ {mustBeEmptyOrScalar(CurrentZ)}
        CurrentZIndex {mustBeEmptyOrScalar(CurrentZIndex)}
        
        CurrentNucleusCellNo {mustBeEmptyOrScalar(CurrentNucleusCellNo)}
        CurrentNucleus {mustBeEmptyOrScalar(CurrentNucleus)}
        PreviousNucleus {mustBeEmptyOrScalar(PreviousNucleus)}
        lastNucleus {mustBeEmptyOrScalar(lastNucleus)}
       
        
        CurrentChannel {mustBeEmptyOrScalar(CurrentChannel)}
        CurrentChannelIndex {mustBeEmptyOrScalar(CurrentChannelIndex)}
        PreviousChannel {mustBeEmptyOrScalar(PreviousChannel)}
        PreviousChannelIndex {mustBeEmptyOrScalar(PreviousChannelIndex)}
        inputChannel % replaces coatChannel in CPTState
        
        FrameIndicesToFit
        Coefficients
        fitApproved
        lineFitted
        
        ZoomMode
        GlobalZoomMode
        xForZoom
        yForZoom
        
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
            this.CurrentZ = round(this.ZSlices / 2);
            
            this.CurrentNucleus = 1;
            this.PreviousNucleus = 1;
            this.CurrentNucleusCellNo = 1;
            this.lastNucleus = 0; %this gets flagged if there's a drop to one particle within the Particles structure.
            
            this.CurrentChannel = 1;
            this.CurrentChannelIndex = 1;
            this.PreviousChannel = this.CurrentChannel;
            this.PreviousChannelIndex = this.CurrentChannelIndex;
            
            
            this.FrameIndicesToFit = 0; % index of the current particle that were used for fitting
            this.Coefficients = []; % coefficients of the fitted line
            this.fitApproved = 0; %JP: I think functions should use this instead of calculating fitApproved on their own
            this.lineFitted = 0; % equals 1 if a line has been fitted
            
            this.ZoomMode = 0;
            this.GlobalZoomMode = 0;
            this.xForZoom = 0;
            this.yForZoom = 0;
            
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
            numValidFrames = length({this.Ellipses});
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
              currentFrameIndex = find(this.getCurrentSchnitzCell().frames == this.CurrentFrame);
              currentNucleusX = getCurrentSchnitzCell().cenx(currentFrameIndex);
          end
%         
%         function currentYDoG = getCurrentYDoG(this)
%             currentFit = this.getCurrentParticleFit();
%             currentYDoG = double(currentFit.yDoG(this.CurrentZIndex));
%         end
          function currentNucleusY = getCurrentY(this)
              currentFrameIndex = find(this.getCurrentSchnitzCell().frames == this.CurrentFrame);
              currentNucleusY = getCurrentSchnitzCell().ceny(currentFrameIndex);
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
%             this.CurrentZIndex = find(this.getCurrentParticleFit().z == this.CurrentZ);
%         end
        
        function maxZIndex = getMaxZIndex(this)
            currentFrameIndex = find(this.getCurrentSchnitzCell().frames == this.CurrentFrame);
            maxZIndex = find(this.getCurrentSchnitzCell().Fluo(currentFrameIndex,:) == max(this.getCurrentSchnitzCell().Fluo(currentFrameIndex,:)));
        end
        
        function currentNucleusCellNo = getCurrentNucleusCellNo(this)
            currentNucleusCellNo = this.getCurrentSchnitzCell().cellno(...
                this.getCurrentSchnitzCell().frames == this.CurrentFrame);
        end
        
        function [xApproved, yApproved] = getApprovedSchnitzCells(this, x, y)
            IndexApprovedSchnitzCells = [];
            numNuclei = this.numNuclei();
            currentSchnitzCells = this.schnitzcells();
            
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
            currentSchnitzCells = this.schnitzcells();
            
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
            currentSchnitzCells = this.schnitzcells();
            
            
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
        
       
        % not exactly sure what this does (GM 9/7/20)
        % HASN'T BEEN CHANGED AT ALL and neither has anything that comes
        % after it. 
        function processImageMatrices(this, nFrames,...
                nSlices, blankImage, currentNC,...
                ncFramesFull, movieMat, maxMat)
            
            if strcmpi(this.projectionMode, 'None')
                
                if ~isempty(movieMat)
                    this.ImageMat = movieMat(:, :, this.CurrentZ,...
                        this.CurrentFrame, this.CurrentChannel);
                else
                    this.ImageMat = getMovieSlice(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel, this.CurrentZ);
                end
                
                %to have a 3x3 square of time and z images on the screen. 
%                 if multiView
%                     for z = 1:-1:-1
%                         for f = -1:1:1
%                             if any( 1:nSlices == this.CurrentZ + z) &&...
%                                     any( 1:nFrames == this.CurrentFrame + f)
%                                 if ~isempty(movieMat)
%                                     this.multiImage{z+2, f+2} =...
%                                         movieMat(:, :, this.CurrentZ+z,...
%                                         this.CurrentFrame+f, this.CurrentChannel);
%                                 else
%                                     this.multiImage{z+2, f+2} =...
%                                         getMovieSlice(this.liveExperiment, this.CurrentFrame+f,...
%                                         this.CurrentChannel,...
%                                         this.CurrentZ+z);
%                                 end
%                             else
%                                 this.multiImage{z+2, f+2} = blankImage;
%                             end
%                         end % loop over frames
%                     end % loop over z slices
%                 end
                
            elseif strcmpi(this.projectionMode, 'Max Z')
                
                if ~isempty(maxMat)
                    if nFrames > 1
                        this.ImageMat = maxMat(:, :, this.CurrentFrame, this.CurrentChannel);
                    else
                        this.ImageMat = maxMat;
                    end
                else
                    
                    imStack = getMovieFrame(this.liveExperiment,...
                        this.CurrentFrame, this.CurrentChannel);
                    
                    this.ImageMat = squeeze(max(imStack, [], 3));
                    
                end
                
                %not currently supported when loading single stacks
            elseif strcmpi(this.projectionMode, 'Max Z and Time')
                if isempty(this.maxTimeCell)
                    this.ImageMat = max(max(movieMat(...
                        :,:,:, ncFramesFull(currentNC):ncFramesFull(currentNC+1),...
                        this.CurrentChannel), [], 2), [], 3); % ch z t x y
                end
            end
        end
        
        % NOT SURE WHAT THESE & and y INPUTS ARE
        function [xTrace, yTrace] = getXYTraces(this, x, y)
            xTrace = x(this.CurrentNucleusCellNo);
            yTrace = y(this.CurrentNucleusCellNo);
        end
        
        % NOT SURE WHAT THESE & and y INPUTS ARE
        function updateZIndex(this, x, y, z)
            [xTrace, yTrace] = this.getXYTraces(x, y);
            
            if (~isempty(xTrace)) && (~this.ManualZFlag)
                this.CurrentZ = z(this.CurrentNucleusCellNo);
                this.CurrentZIndex = find(this.schnitzcells(this.CurrentNucleus).Fluo(this.CurrentFrame,:)...
                    == max(this.schnitzcells(this.CurrentNucleus).Fluo(this.CurrentFrame,:)));
                this.ManualZFlag = 0;
            end
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