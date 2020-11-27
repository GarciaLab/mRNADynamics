function Particles = addNucleusProbabilities(liveExperiment, trackingOptions, FrameInfo, Particles)
    
    % path to files 
    nucleusProbDirFinal = [liveExperiment.procFolder 'nucleusProbabilityMapsFull' filesep];
    
    % Get frame dimensions
    xDim = FrameInfo(1).PixelsPerLine;
    yDim = FrameInfo(1).LinesPerFrame;
    zDim = FrameInfo(1).NumberSlices+2;

    for Channel = 1:trackingOptions.NCh
        ParticlesCh = Particles{Channel};          
        % generate indexing vectors
        FrameVec = [];
        ParticleIndexVec = [];
        ParticleSubIndexVec = [];
        for p = 1:length(ParticlesCh)
            FrameVec = [FrameVec ParticlesCh(p).Frame];
            ParticleIndexVec = [ParticleIndexVec repelem(p,length(ParticlesCh(p).Frame))];
            ParticleSubIndexVec = [ParticleSubIndexVec 1:length(ParticlesCh(p).Frame)];
            % initialize field
            ParticlesCh(p).nucleusProbability = NaN(size(ParticlesCh(p).Frame));
        end
        % generate position vectors
        xPosVec = round(vertcat(ParticlesCh.xPos));
        xPosVec(xPosVec<1) = 1;
        xPosVec(xPosVec>xDim) = xDim;
        yPosVec = round(vertcat(ParticlesCh.yPos));
        yPosVec(yPosVec<1) = 1;
        yPosVec(yPosVec>yDim) = yDim;
        zPosVec = round(vertcat(ParticlesCh.zPos));
        zPosVec(zPosVec<1) = 1;
        zPosVec(zPosVec>zDim) = zDim;

        % convert these to linear indices
        linIndVec = sub2ind([yDim,xDim,zDim],yPosVec,xPosVec,zPosVec);
        FrameIndex = unique(FrameVec);
        % iterate through frames and determine likelihood of each
        % detection event of being inside a nucleus
        probCell = cell(size(FrameInfo));
        parfor frame = 1:max(FrameIndex)
            % load probability stack 
            if any(FrameVec==frame)
                readName = ['prob' liveExperiment.Prefix '_' sprintf('%03d',frame) '_ch00.tif'];
                probStack = double(imreadStack([nucleusProbDirFinal  readName]));
                probStack(probStack>1e4) = 1e4;
                probStack = probStack/1e4;
                % now iterate through particles and assign              
                probCell{frame} = probStack(linIndVec(FrameVec==frame));       
            end
        end
        % record
        for frame = FrameIndex
            frameIndices = find(FrameVec==frame);
            for frameInd = 1:length(frameIndices)                  
                % record
                ParticlesCh(ParticleIndexVec(frameIndices(frameInd))).nucleusProbability(ParticleSubIndexVec(frameIndices(frameInd)))...
                            = probCell{frame}(frameInd);
            end
        end
        Particles{Channel} = ParticlesCh;
    end