function Particles=JoinParticleTraces(OriginalParticle,ClickedParticle,Particles,liveExperiment,Spots)

%This function joins two particle traces and renumbers all particles in the
%Particles structure accordingly

%Transfer the information to the original particle
vectorFields = {'xPos','yPos','zPos','zPosDetrended','FrameApproved','Index'};
[Particles(OriginalParticle).Frame, sortOrder] = sort([Particles(OriginalParticle).Frame,Particles(ClickedParticle).Frame]);
% generate logical vector for ease of indexing
newFrameFlags = ismember(Particles(OriginalParticle).Frame,Particles(ClickedParticle).Frame);

for v = 1:length(vectorFields)
    vField = vectorFields{v};
    if isfield(Particles,vField) && v<=4
        vec = vertcat(Particles(OriginalParticle).(vField),Particles(ClickedParticle).(vField));
        Particles(OriginalParticle).(vField) = vec(sortOrder);
    elseif isfield(Particles,vField)
        vec = [Particles(OriginalParticle).(vField),Particles(ClickedParticle).(vField)];
        Particles(OriginalParticle).(vField) = vec(sortOrder);
    end
end
if isfield(Particles,'FrameApproved')
    Particles(OriginalParticle).FrameApproved=logical(Particles(OriginalParticle).FrameApproved);
end

Particles(OriginalParticle).Approved=0;
%Particles(OriginalParticle).nc=[Particles(OriginalParticle).nc,Particles(ClickedParticle).nc];

%Now, get rid of the clicked particle
Particles=Particles([1:ClickedParticle-1,ClickedParticle+1:end]);

%Deals with the indexing changing because of the removal of
%the old particle.
if ClickedParticle<OriginalParticle
   OriginalParticle=OriginalParticle-1;
end

%Perform Kalman filtering
trackingOptions = getTrackingOptions(liveExperiment);
trackingStruct.MeasurementVec = [Particles(OriginalParticle).xPos Particles(OriginalParticle).yPos...
                    Particles(OriginalParticle).zPosDetrended];
trackingStruct.Frame = Particles(OriginalParticle).Frame;
FluoVec = [];
for f = 1:length(trackingStruct.Frame)
    spot = Spots(trackingStruct.Frame(f)).Fits(Particles(OriginalParticle).Index(f));
    brightestZ = spot.brightestZ;
    brightestZIndex = spot.z == brightestZ;
    FluoVec(f) = spot.FixedAreaIntensity(brightestZIndex);
end
trackingStruct.MeasurementVec(:,end+1) = FluoVec./trackingOptions.kalmanOptions.fluoFactor;

% make particle path predictions
Particles(OriginalParticle) = pathPrediction(Particles(OriginalParticle), trackingStruct, trackingOptions, trackingOptions.kalmanOptions);

% update QC info
Particles(OriginalParticle).logL = nansum(Particles(OriginalParticle).logLDistance(:,1:3),2);
Particles(OriginalParticle).logLMean = nanmean(Particles(OriginalParticle).logL);

qcVariables = {'NucleusBoundaryFlags', 'SpotlogLFlags', 'earlyFlags'};
for q = 1:length(qcVariables)
    qcVec = false(size(newFrameFlags));
    qcVec(~newFrameFlags) = Particles(OriginalParticle).(qcVariables{q});
    Particles(OriginalParticle).(qcVariables{q}) = qcVec;
end
Particles(OriginalParticle).FractionFlagged = mean(Particles(OriginalParticle).FrameApproved);

% update qc scores
qcScoreArray = zeros(length(newFrameFlags),size(Particles(OriginalParticle).qcScoreArray,2));
qcScoreArray(~newFrameFlags,:) = Particles(OriginalParticle).qcScoreArray;
Particles(OriginalParticle).qcScoreArray = qcScoreArray;

% just set nucleus probability to Inf 
ncProbs = Inf(size(newFrameFlags));
ncProbs(~newFrameFlags) = Particles(OriginalParticle).nucleusProbability;
Particles(OriginalParticle).nucleusProbability = ncProbs;

end