function indMit = generateIndMit(anaphaseFrames, nFrames)
% Subfunction for TrackNuclei that prepares anaphase timing array
% for use in MainTracking.

if iscolumn(anaphaseFrames)
    anaphaseFrames = anaphaseFrames';
end

%This checks whether all ncs have been defined
if length(anaphaseFrames)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
end

if length( find(isnan(anaphaseFrames))) ==...
        length(anaphaseFrames) || length(anaphaseFrames) < 6
    error('Have the ncs been defined in MovieDatabase or anaphaseFrames.mat?')
end

%Checking for an edge case that happens on occassion and leads to _very_
%opaque errors down the road
for k = 2:length(anaphaseFrames)
    if ~isnan(anaphaseFrames(k)) && isnan(anaphaseFrames(k-1))
        error(['NaN entry cannot preceed non-NaN mitosis entries. Please correct ',...
            'this in chooseanaphaseframes or moviedatabase.']);
    end
end

anaphaseFrames=anaphaseFrames(anaphaseFrames~=0);

%Note that I'm adding a range of two frames frames before and after the
%determines mitosis
indMit=[anaphaseFrames'-2,anaphaseFrames'+2];

%Make sure no indices are negative. This could happen is the nuclear cycle
%started at frame 1, for example.
indMit(indMit<1)=1;

%Check whether nc14 occurred very close to the end of the movie. For those
%frames we'll move the boundary for tracking purposes
indMit(indMit>=nFrames)=indMit(indMit>=nFrames)-2;

%If indMit is empty this means that we didn't see any mitosis. Deal with
%this by assigning the first frames to it
if isempty(indMit)
    indMit=[1,2];
    %If we don't have nc14 we'll fool the code into thinking that the last
    %frame of the movie was nc14
elseif isnan(indMit(end,1))
    indMit(end,1)=nFrames-3;
    indMit(end,2)=nFrames-2;
end

end