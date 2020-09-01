function [hmmModel, pathArray, sigmaArray] = simulatePathsWrapper(Particle,hmmModel,ncFrameFilter)

  % calculate forward and backward frames to calculate
  frameVec = Particle.Frame; 
  ncFrames = find(ncFrameFilter);
  fwdFrames = ncFrames(ncFrames>frameVec(end));
  bkdFrames = ncFrames(ncFrames<frameVec(1));
  % initialize arrays to store path info 
  pathArray = NaN(length(ncFrameFilter),length(hmmModel));
  sigmaArray = NaN(length(ncFrameFilter),length(hmmModel));
  % Loop through each variable. Generate decoded jump state vector and
  % use this to predict past and future particle positions
  for h = 1:length(hmmModel)
    varName = hmmModel(h).var; % current variable name
    varPos = Particle.(varName); % extract 
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
    pathArray(ncFrameFilter,h) = [pdBkdMean+varPos(1) varPos pdFwdMean+varPos(end)];
    hmmModel(h).pathVec = pathArray(:,h); 
    sigmaArray(ncFrameFilter,h) = [pdBkdSE zeros(size(varPos)) pdFwdSE];
    hmmModel(h).sigmaVec = sigmaArray(:,h);
    
%     hmmModel(h).dfVec = [abs(bkdFrames-frameVec(1)) zeros(size(varPos)) fwdFrames-frameVec(end)];
  end