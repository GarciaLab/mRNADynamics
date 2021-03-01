classdef CPTState < handle
    
    properties
        liveExperiment
        
        Spots
        Particles
        SpotFilter
        schnitzcells
        FrameInfo
        ImageMat
        storedTimeProjection
        multiImage
        maxTimeCell
        
        UseTwinTraces 
        HideSingleSliceTrace
        
        PlotInputChannel
        
        
        Ellipses
        nucleiModified
        
        Frames
        CurrentFrame {mustBeEmptyOrScalar(CurrentFrame)}
        PreviousFrame {mustBeEmptyOrScalar(PreviousFrame)}
        
        ManualZFlag
        ZSlices
        CurrentZ {mustBeEmptyOrScalar(CurrentZ)}
        CurrentZIndex {mustBeEmptyOrScalar(CurrentZIndex)}
        
        CurrentParticleIndex {mustBeEmptyOrScalar(CurrentParticleIndex)}
        CurrentParticle {mustBeEmptyOrScalar(CurrentParticle)}
        PreviousParticle {mustBeEmptyOrScalar(PreviousParticle)}
        lastParticle {mustBeEmptyOrScalar(lastParticle)}
        
        TwinParticleIndex {mustBeEmptyOrScalar(TwinParticleIndex)}
        TwinParticle {mustBeEmptyOrScalar(TwinParticle)}
        PreviousTwinParticle {mustBeEmptyOrScalar(PreviousTwinParticle)}
        
        CurrentChannel {mustBeEmptyOrScalar(CurrentChannel)}
        CurrentChannelIndex {mustBeEmptyOrScalar(CurrentChannelIndex)}
        PreviousChannel {mustBeEmptyOrScalar(PreviousChannel)}
        PreviousChannelIndex {mustBeEmptyOrScalar(PreviousChannelIndex)}
        coatChannel
        
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
        plot3DGauss
        
        projectionMode
        
        UseCompiledParticles
    end
    
    methods
        function this = CPTState(liveExperiment, Spots, Particles, SpotFilter, schnitzcells, Ellipses,...
                FrameInfo, UseHistoneOverlay, nWorkers, plot3DGauss, projectionMode, UseTwinTraces, PlotInputChannel, UseCompiledParticles)
            
            this.liveExperiment = liveExperiment;
            this.Spots = Spots;
            this.Particles = Particles;
            this.SpotFilter = SpotFilter;
            this.schnitzcells = schnitzcells;
            this.FrameInfo = FrameInfo;
            this.ImageMat = [];
            this.storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called
            this.multiImage = {};
            this.maxTimeCell = [];
            
            if exist('UseTwinTraces', 'var')
                this.UseTwinTraces = UseTwinTraces;
            else
                this.UseTwinTraces = false;
            end
            
            if exist('PlotInputChannel', 'var')
                this.PlotInputChannel = PlotInputChannel;
            else
                this.PlotInputChannel = false;
            end
            
            this.Ellipses = Ellipses;
            this.nucleiModified = false;
            
            this.Frames = [];
            this.CurrentFrame = 0;
            this.PreviousFrame = this.CurrentFrame;
            
            this.ManualZFlag = 0;
            numberZSlices = this.FrameInfo(1).NumberSlices;
            this.ZSlices = numberZSlices + 2; %Note that the blank slices are included
            this.CurrentZ = round(this.ZSlices / 2);
            
            this.CurrentParticle = 1;
            this.PreviousParticle = 1;
            this.lastParticle = 0; %this gets flagged if there's a drop to one particle within the Particles structure.
            
            
            this.CurrentChannel = 1;
            this.CurrentChannelIndex = 1;
            this.PreviousChannel = this.CurrentChannel;
            this.PreviousChannelIndex = this.CurrentChannelIndex;
            
            % Added by G. Martini to support 2 spot data analysis
            this.TwinParticle = findTwinParticle(this.CurrentParticle, this.Particles, this.CurrentChannelIndex);
            this.PreviousTwinParticle = 1;
            
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
            this.plot3DGauss = plot3DGauss;
            
            this.projectionMode = projectionMode;
            
            if exist('UseCompiledParticles', 'var')
                this.UseCompiledParticles = UseCompiledParticles;
            else
                this.UseCompiledParticles = false;
            end
            
            if this.UseTwinTraces
                this.HideSingleSliceTrace = true;
            else
                this.HideSingleSliceTrace = false;
            end
        end
        
        function numParticles = numParticles(this)
            numParticles = length(this.Particles{this.CurrentChannelIndex});
        end
        
        function numValidFrames = numValidFrames(this)
            numValidFrames = length({this.Spots{1}.Fits});
        end
        
        function currentSpots = getCurrentChannelSpots(this)
            currentSpots = this.Spots{this.CurrentChannelIndex};
        end
        
        function currentFrameSpots = getCurrentFrameSpots(this)
            currentSpots = this.getCurrentChannelSpots();
            currentFrameSpots = currentSpots(this.CurrentFrame);
        end
        
        function currentParticles = getCurrentChannelParticles(this)
            currentParticles = this.Particles{this.CurrentChannelIndex};
        end
        
        function currentParticle = getCurrentParticle(this)
            currentParticles = this.getCurrentChannelParticles();
            currentParticle = currentParticles(this.CurrentParticle);
        end
        
        function twinParticle = getTwinParticle(this)
            currentParticles = this.getCurrentChannelParticles();
            twinParticle = currentParticles(this.TwinParticle);
        end
        
     
        
        function currentParticleFit = getCurrentParticleFit(this)
            currentFrameSpots = this.getCurrentFrameSpots();
            currentParticleFit = currentFrameSpots.Fits(this.CurrentParticleIndex);
        end
        
        function currentXDoG = getCurrentXDoG(this)
            currentFit = this.getCurrentParticleFit();
            currentXDoG = double(currentFit.xDoG(this.CurrentZIndex));
        end
        
        function currentYDoG = getCurrentYDoG(this)
            currentFit = this.getCurrentParticleFit();
            currentYDoG = double(currentFit.yDoG(this.CurrentZIndex));
        end
        
        function currentXFit = getCurrentXFit(this)
            currentFit = this.getCurrentParticleFit();
            currentXFit = double(currentFit.xFit(this.CurrentZIndex));
        end
        
        function currentYFit = getCurrentYFit(this)
            currentFit = this.getCurrentParticleFit();
            currentYFit = double(currentFit.yFit(this.CurrentZIndex));
        end
        
        function updateCurrentZIndex(this)
            this.CurrentZIndex = find(this.getCurrentParticleFit().z == this.CurrentZ);
        end
        
        function maxZIndex = getMaxZIndex(this)
            maxZIndex = find(this.getCurrentParticleFit().z == this.getCurrentParticleFit().brightestZ);
        end
        
        function currentParticleIndex = getCurrentParticleIndex(this)
            currentParticleIndex = this.getCurrentParticle().Index(...
                this.getCurrentParticle().Frame == this.CurrentFrame);
        end
        
        function twinParticleIndex = getTwinParticleIndex(this)
            if this.UseTwinTraces
                TwinParticle = this.getTwinParticle();
                if ~isempty(TwinParticle)
                    twinParticleIndex = TwinParticle.Index(...
                        TwinParticle.Frame == this.CurrentFrame);
                else
                    twinParticleIndex = [];
                end
            else
                twinParticleIndex = [];
            end
        end
        
        function [xApproved, yApproved] = getApprovedParticles(this, x, y)
            IndexApprovedParticles = [];
            numParticles = this.numParticles();
            currentChannelParticles = this.getCurrentChannelParticles();
            
            for i = 1:numParticles
                if sum(currentChannelParticles(i).Frame == this.CurrentFrame) &&...
                        sum(currentChannelParticles(i).Approved == 1)
                    IndexApprovedParticles = [IndexApprovedParticles,...
                        currentChannelParticles(i).Index(currentChannelParticles(i).Frame == this.CurrentFrame)];
                end
                
            end
            
            xApproved = x(IndexApprovedParticles);
            yApproved = y(IndexApprovedParticles);
        end
        
        function [xDisapproved, yDisapproved] = getDisapprovedParticles(this, x, y)
            IndexDisapprovedParticles=[];
            numParticles = this.numParticles();
            currentChannelParticles = this.getCurrentChannelParticles();
            
            for i = 1:numParticles
                if sum(currentChannelParticles(i).Frame == this.CurrentFrame) && sum(currentChannelParticles(i).Approved == -1)
                    IndexDisapprovedParticles=[IndexDisapprovedParticles,...
                        currentChannelParticles(i).Index(currentChannelParticles(i).Frame == this.CurrentFrame)];
                end
            end
            
            xDisapproved = x(IndexDisapprovedParticles);
            yDisapproved = y(IndexDisapprovedParticles);
        end
        
        function [xNonFlagged, yNonFlagged] = getNonFlaggedParticles(this, x, y)
            IndexNonFlaggedParticles = [];
            numParticles = this.numParticles();
            currentChannelParticles = this.getCurrentChannelParticles();
            
            for i = 1:numParticles
                if sum(currentChannelParticles(i).Frame == this.CurrentFrame) &&...
                        ~(sum(currentChannelParticles(i).Approved == -1) || sum(currentChannelParticles(i).Approved == 1))
                    IndexNonFlaggedParticles=[IndexNonFlaggedParticles,...
                        currentChannelParticles(i).Index(currentChannelParticles(i).Frame == this.CurrentFrame)];
                end
            end
            
            xNonFlagged = x(IndexNonFlaggedParticles);
            yNonFlagged = y(IndexNonFlaggedParticles);
        end
        
        function [DaughterE, DaughterD, Mother] = getMotherDaughters(this)
            if isfield(this.schnitzcells, 'E')
                DaughterE = this.schnitzcells(this.getCurrentParticle().Nucleus).E;
                DaughterD = this.schnitzcells(this.getCurrentParticle().Nucleus).D;
                Mother = this.schnitzcells(this.getCurrentParticle().Nucleus).P;
            else
                DaughterE = 0;
                DaughterD = 0;
                Mother = 0;
            end
        end
        
        function processImageMatrices(this, multiView, nFrames,...
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
                if multiView
                    for z = 1:-1:-1
                        for f = -1:1:1
                            if any( 1:nSlices == this.CurrentZ + z) &&...
                                    any( 1:nFrames == this.CurrentFrame + f)
                                if ~isempty(movieMat)
                                    this.multiImage{z+2, f+2} =...
                                        movieMat(:, :, this.CurrentZ+z,...
                                        this.CurrentFrame+f, this.CurrentChannel);
                                else
                                    this.multiImage{z+2, f+2} =...
                                        getMovieSlice(this.liveExperiment, this.CurrentFrame+f,...
                                        this.CurrentChannel,...
                                        this.CurrentZ+z);
                                end
                            else
                                this.multiImage{z+2, f+2} = blankImage;
                            end
                        end % loop over frames
                    end % loop over z slices
                end
                
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
        
        function [xTrace, yTrace] = getXYTraces(this, x, y)
            xTrace = x(this.CurrentParticleIndex);
            yTrace = y(this.CurrentParticleIndex);
        end
        
        function updateZIndex(this, x, y, z)
            [xTrace, yTrace] = this.getXYTraces(x, y);
            
            if (~isempty(xTrace)) && (~this.ManualZFlag)
                this.CurrentZ = z(this.CurrentParticleIndex);
                this.CurrentZIndex = find(this.getCurrentParticleFit().z == this.CurrentZ);
                this.ManualZFlag = 0;
            end
        end
        
        function updateCurrentParticleIndex(this)
            if this.UseTwinTraces
                this.TwinParticle = findTwinParticle(this.CurrentParticle, this.Particles, this.CurrentChannelIndex);
            end
            this.CurrentParticleIndex = this.getCurrentParticleIndex();
            this.TwinParticleIndex = this.getTwinParticleIndex();
        end
    end
end

function mustBeEmptyOrScalar(in)

if numel(in) > 1
    error(['Value assigned to property is not '...
        'scalar or empty']);
end

end