function RawParticles = track03PredictParticlePaths(...
                          RawParticles, FrameInfo, retrack, displayFigures)
                        
  % This script uses inferred GHMM models to (a) infer "true" motion states
  % for observed particle frames and to predict particles positions before
  % and after first and last observed frames. These hypothetical tracks can
  % then be used to stitch together particle fragments
  
  
  NCh = length(RawParticles);
  ncVec = [FrameInfo.nc];
  frameIndex = 1:length(ncVec);
  
  for Channel = 1:NCh
    wb = waitbar(0,['Simulating particle paths (channel ' num2str(Channel) ')']);
    for p = 1:length(RawParticles{Channel})
      waitbar(p/length(RawParticles{Channel}),wb);
      % calculate forward and backward frames to calculate
      frameVec = RawParticles{Channel}(p).Frame;
      nc = ncVec(frameVec(1)==frameIndex);
      ncFrames = frameIndex(ncVec==nc);
      fwdFrames = ncFrames(ncFrames>frameVec(end));
      bkdFrames = ncFrames(ncFrames<frameVec(1));
      % extract hmm model
      hmmModel = RawParticles{Channel}(p).hmmModel;
      % Loop through each variable. Generate decoded jump state vector and
      % use this to predict past and future particle positions
      for h = 1:length(hmmModel)
        varName = hmmModel(h).var; % current variable name
        varPos = RawParticles{Channel}(p).(varName); % extract 
        varDeltas = diff(varPos); % The HMM operates on jumps, not positions

        if length(varDeltas) > 1
          % get vector of data likelihoods
          B = mixgauss_prob(varDeltas, hmmModel(h).Mu, hmmModel(h).Sigma);
          [fwdProbs, bkdProbs, ~,  ~, ~] = ...
            fwdback(hmmModel(h).Prior, hmmModel(h).Transmat, B);

          % calculate jump probabilities for each frame
          ssProbArray = log(fwdProbs) + log(bkdProbs); % fwd
          ssProbArray = ssProbArray - logsumexp(ssProbArray);
          ssProbArray = exp(ssProbArray);
        else
          ssProbArray = hmmModel(h).Prior;
        end      
        
        % calculate forward and backward probabilities   
        Mu = hmmModel(h).Mu';
        Sigma = reshape(hmmModel(h).Sigma,[],1);
        A = hmmModel(h).Transmat;
        % forward
        fwdArray = NaN(length(Mu),length(fwdFrames));
        prev = ssProbArray(:,end);
        for f = 1:length(fwdFrames)
          fwdArray(:,f) = A'*prev;
          prev = fwdArray(:,f);
        end
        % backward        
        bkdArray = NaN(length(Mu),length(bkdFrames));
        post = ssProbArray(:,1);
        for f = length(bkdFrames):-1:1
          bkdArray(:,f) = sum(A'.*post',2)/sum(sum(A'.*post',2));
          post = bkdArray(:,f);
        end
        
        % calculate expected mean and variance
        pdFwdMean = cumsum(sum(fwdArray.*Mu));
        pdFwdSE1 = cumsum(sum(fwdArray.*Sigma)); % expecation of variance
        pdFwdSE2 = cumsum(sum(fwdArray.*(Mu.^2))-sum((fwdArray.*Mu).^2)); % variance of expectation
        pdFwdSE = sqrt(pdFwdSE1 + pdFwdSE2);
                
        pdBkdMean = fliplr(cumsum(fliplr(sum(bkdArray.*-Mu))));
        pdBkdSE1 = cumsum(fliplr(sum(bkdArray.*Sigma)));
        pdBkdSE2 = cumsum(fliplr(sum(bkdArray.*(Mu.^2))-sum((bkdArray.*Mu).^2))); % variance of expectation
        pdBkdSE = fliplr(sqrt(pdBkdSE1 + pdBkdSE2));

        % generate aggregate path and helper vectors
        hmmModel(h).pathVec = [pdBkdMean+varPos(1) varPos pdFwdMean+varPos(end)];
        hmmModel(h).sigmaVec = [pdBkdSE zeros(size(varPos)) pdFwdSE];
        hmmModel(h).dfVec = [abs(bkdFrames-frameVec(1)) zeros(size(varPos)) fwdFrames-frameVec(end)];
      end
      RawParticles{Channel}(p).hmmModel = hmmModel;
    end
    close(wb);
  end
  