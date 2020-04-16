function [ xy, mapping, nuclei ] = trackMitosis(FrameInfo, hisMat, firstFrameNumber, lastFrameNumber, shifts, diameter, embryoMask, xy, manualMapping, nuclei, varargin )
%BACKTRACKMITOSIS This function tracks the nuclei through mitosis. The main
% difference with respect to the tracking through interphases is that here,
% nuclei can be mapped to two nuclei on the next frame.
%   Detailed explanation goes here

if numel(varargin) > 0
    h_waitbar_tracking = varargin{1};
end

space_resolution = getDefaultParameters(FrameInfo,'space resolution');
totalNumberOfFrames = size(hisMat,3);
maxNucleusStep = 40*0.22/space_resolution;
maxShiftCorrection = getDefaultParameters(FrameInfo,'max Shift Correction', 'trackToTheNextFrame')*diameter/space_resolution;
nFrames = (lastFrameNumber-firstFrameNumber+1);
% edgeClearance = getDefaultParameters(FrameInfo,'edge clearance')*diameter/space_resolution;

mapping = {};

% If nuclei are not provided, segment the images.
if ~exist('xy','var') || isempty(xy)
    for j = 1:nFrames
        
        tmp = findNuclei(FrameInfo,hisMat,firstFrameNumber-1+j,diameter,embryoMask);
        
        xy{j} = tmp;%+repmat(sum(shifts(firstFrameNumber:firstFrameNumber-1+j,:),1),size(tmp,1),1);
        
    end
else
    if numel(xy) ~= nFrames
        error(['Provided nuclei cell array is not the right size. Nuclei for ' num2str(nFrames) ' are needed (frame ' num2str(firstFrameNumber) ' to ' num2str(lastFrameNumber) ' ) but ' num2str(numel(xy)) ' were provided.'])
    end
end
if exist('nuclei','var') && ~isempty(nuclei)
    fill_nuclei_struct = true;
    nucleiIndices = nan(size(xy{1},1),1);
    for j = 1:size(xy{1},1)
        for jj = 1:numel(nuclei)
            if nuclei(jj).indXY(firstFrameNumber) == j
                nucleiIndices(j) = jj;
                break
            end
        end
    end
    if any(isnan(nucleiIndices))
        error(' NAN ')
    end
else
    nucleiIndices = 1:size(xy{1},1);
    fill_nuclei_struct = false;
end

% Apply the threshold on the shift correction:
shifts = max(min(shifts,maxShiftCorrection),-maxShiftCorrection);

alreadyDivided_jm1 = false(size(xy{1},1),1);

for j = 2:nFrames%:-1:2
    
    %
    % 1. Initialize variables
    %
    mapping_tmp = zeros(size(xy{j-1},1),2); % temporary array where found matches are stored.
    
    tmpAlreadyDivided_jm1 = alreadyDivided_jm1;
    tmpAlreadyDivided_j = false(size(xy{j},1),1);
    mapping_tmp(tmpAlreadyDivided_jm1,2) = -1;
    
    inverse_mapping = zeros(size(xy{j},1),1);
    
    keep_looping = true; % true as long as some nuclei were mapped together during the last loop.

    % List of nuclei that could still be mapped
    remainingNuclei1 = xy{j-1};
    remainingNuclei2 = xy{j};
    % List of nuclei that could still be mapped, corrected for bulk shift.
    shiftedXY1 = xy{j-1}+repmat(shifts(firstFrameNumber+j-2,:),size(xy{j-1},1),1);
    
    % These 'indices' vectors will help to keep track of the original index
    % of each nucleus. (During the algorithm, nuclei that were successfully
    % mapped are taken out of the list of remaining nuclei, which
    % introduces ambiguity on their indices).
    indices1 = 1:size(remainingNuclei1,1);
    indices2 = 1:size(remainingNuclei2,1);
    
    mappedNuc1 = []; % List of nuclei indices for which matches were found on the next frame
    mappedNuc2 = []; % List of nuclei indices for which matches were found on the previous frame
    
    if exist('manualMapping','var') && ~isempty(manualMapping{j-1}) && any(manualMapping{j-1}(:))
        
        % Extract the manual corrections
        mapp = manualMapping{j-1};
        ind0 = any(mapp < 0,2); % Values below 0 prevent further mapping, i.e. if one of the two mapping indices is -1, the nucleus will only be mapped to one nucleus on the next frame and if both of them are -1, the nucleus will disappear on the next frame.
        ind1 = find(any(mapp > 0,2)); % When two values are provided, division is enforced.
        ind2 = mapp(any(mapp > 0,2),:); % Indices of the two daughters.
        
        % Save the manual mapping data as found matches.
        mapping_tmp(ind1,:) = ind2;
        inverse_mapping(ind2(:,1)) = ind1;
        inverse_mapping(ind2(ind2(:,2)>0,2)) = ind1(ind2(:,2)>0);
        
        % If necessary, store everything in the 'nuclei' structure
        if fill_nuclei_struct
            for k = 1:numel(ind1)
                
                
                if mapp(ind1(k),2) < 1 % no division, simply map the nucleus to the next frame.
                
                    
                    for jj = 1:numel(nuclei)
                        indx = find(nuclei(jj).indXY(firstFrameNumber+j-2) == ind1(k));
                        if isempty(indx)
                            continue
                        else
                            if numel(indx) > 1
                                error('Several nuclei have the same position.');
                            else
                                nuclei(jj).indXY(firstFrameNumber+j-1) = ind2(indx);
                            end
                        end
                    end
                    
                    
                else if mapp(ind1(k),2) > 0 % Division. This nucleus ends here and two daughters start on the next frame.
                    
                        
                            indx = nucleiIndices(ind1(k));
                                if numel(indx) > 1
                                    error('Several nuclei have the same position.');
                                else
                                    indE = numel(nuclei)+1;
                                    indD = numel(nuclei)+2;
                                    nuclei(indE).position = nan(totalNumberOfFrames,2);
                                    nuclei(indE).indXY = zeros(totalNumberOfFrames,1);
                                    nuclei(indD).position = nan(totalNumberOfFrames,2);
                                    nuclei(indD).indXY = zeros(totalNumberOfFrames,1);
                                    nuclei(nucleiIndices(ind1(k))).D = ind1;
                                    nuclei(nucleiIndices(ind1(k))).E = ind2;
                                    nuclei(indE).position(firstFrameNumber+j-1,:) = xy{j}(mapping_tmp(ind1(k),1),:);
                                    nuclei(indE).indXY(firstFrameNumber+j-1) = mapping_tmp(ind1(k),1);
                                    nuclei(indE).P = nucleiIndices(ind1(k));
                                    nuclei(indD).position(firstFrameNumber+j-1,:) = xy{j}(mapping_tmp(ind1(k),2),:);
                                    nuclei(indD).indXY(firstFrameNumber+j-1) = mapping_tmp(ind1(k),2);
                                    nuclei(indD).P = nucleiIndices(ind1(k));
                                    indToDelete = [indToDelete ind1(k)];
                                    try
                                    nucInd(mapping_tmp(ind1(k),1)) = indE;
                                    nucInd(mapping_tmp(ind1(k),2)) = indD;
                                    catch
                                        1
                                    end
                                    tmpAlreadyDivided_jm1(ind1(k)) = true;
                                    tmpAlreadyDivided_j(mapping_tmp(ind1(k),1)) = true;
                                    tmpAlreadyDivided_j(mapping_tmp(ind1(k),2)) = true;

                                end
                    end
                end
                
            end
        end
        
        

        
        
        
        mappedNuc1 = [mappedNuc1; remainingNuclei1(ind0,:)];
        mappedNuc1 = [mappedNuc1; remainingNuclei1(ind1,:)];
        mappedNuc2 = [mappedNuc2; remainingNuclei2(ind2(:,1),:)];
        mappedNuc2 = [mappedNuc2; remainingNuclei2(ind2(ind2(:,2) > 0,2),:)];
        
        indMappedTo2Nuclei = all(mapp ~= 0,2);
        remainingNuclei1(indMappedTo2Nuclei | ind0,:) = [];

        remainingNuclei2(ind2(ind2>0),:) = [];
        shiftedXY1(indMappedTo2Nuclei | ind0,:) = [];
        
        
        indices1(indMappedTo2Nuclei | ind0) = [];
        indices2(ind2(ind2>0)) = [];

        
    end
    
    %
    % 2. Loop over all nuclei and at each iteration map the ones that are
    %    the closest matches.
    %
    
    while keep_looping
        
        [dist, ix] = pdist2(remainingNuclei2,...
            shiftedXY1,'euclidean','Smallest',size(remainingNuclei2,1));
        % dist is an m x n array where m is the numer of nuclei in
        % remainingNuclei2 (frame j) and n is the number of nuclei in
        % shiftedXY1 (frame j-1).
        
        indexCorr =  repmat((1:size(dist,2))-1,size(dist,1),1)*size(dist,1);
        unsortedDist = zeros(size(dist));
        unsortedDist(ix+indexCorr) = dist;
        unsortedDist = reshape(unsortedDist,size(dist));
        
        
        % Build correspondences
        
        % First map together nuclei that are mutually the closest
        mc = dist(1:min(2,end),:);
        ic = ix(1:min(2,end),:);
        
        [mr,ir] = min(unsortedDist,[],2);
        
        indc = sub2ind(size(dist),ic,repmat(1:size(ic,2),size(ic,1),1));
        
        indr = sub2ind(size(dist),1:numel(ir),ir');
        
        
        indc(mc>maxNucleusStep) = [];
        indr(mr>maxNucleusStep) = [];
        
        
        Mc = zeros(size(dist));
        Mc(indc) = 1;
        Mr = zeros(size(dist));
        Mr(indr) = 1;
        
        attributedIndices = false(length(remainingNuclei1),1);
        ind = find(Mc & Mr);
        [rowInd,colInd] = ind2sub(size(dist),ind);
        [dummy,order] = sort(mr(rowInd),'ascend');
        try
            rowInd = rowInd(order);
        catch
            1
        end
        colInd = colInd(order);
        for jj = 1:numel(colInd)
            
            try
                if mapping_tmp(indices1(colInd(jj)),1) == 0
                    mapping_tmp(indices1(colInd(jj)),1) = indices2(rowInd(jj));
                else if mapping_tmp(indices1(colInd(jj)),2) == 0
                        mapping_tmp(indices1(colInd(jj)),2) = indices2(rowInd(jj));
                        tmpAlreadyDivided_jm1(indices1(colInd(jj))) = true;
                        tmpAlreadyDivided_j(indices2(rowInd(jj))) = true;
                        tmpAlreadyDivided_j(mapping_tmp(indices1(colInd(jj)),1)) = true;
                    else % mapping_tmp(.,2) is negative, meaning that the nucleus didn't divide in this frame and should only be mapped to one nucleus.
                        continue;
                    end
                end
            catch
                1
            end
            inverse_mapping(indices2(rowInd(jj))) = indices1(colInd(jj));
            attributedIndices(colInd(jj)) = true;
        end
        mappedNuc1 = [mappedNuc1; remainingNuclei1(colInd,:)];
        mappedNuc2 = [mappedNuc2; remainingNuclei2(unique(rowInd),:)];
        
        indMappedTo2Nuclei = all(mapping_tmp(indices1,:) ~= 0,2);
        remainingNuclei1(indMappedTo2Nuclei,:) = [];
        remainingNuclei2(rowInd,:) = [];
        shiftedXY1(indMappedTo2Nuclei,:) = [];
        
        %tmpAlreadyDivided(indMappedTo2Nuclei) = true;
        
        indices1(indMappedTo2Nuclei) = [];
        indices2(rowInd) = [];
        %st = tpaps(mappedNuc1', mappedNuc2',0.95);
        
        %shiftedXY1 = fnval(st,remainingNuclei1')';
        
        if isempty(remainingNuclei1) || isempty(remainingNuclei2) || ~any(attributedIndices)
            keep_looping = false;
        end
        
    end % while loop
    
    alreadyDivided_jm1 = tmpAlreadyDivided_j;
    tmpAlreadyDivided_jm1 = tmpAlreadyDivided_j;
    id = mapping_tmp(:,2) == -1;
    mapping_tmp(id,2) = 0;
    mapping{j-1} = mapping_tmp;
    
    %%%%%%%%
    if fill_nuclei_struct
        % Put the output in the nuclei structure
        mappingInd = j-1;%firstFrameNumber+j-2;
        newFrameInd = mappingInd+1;
        ind = inverse_mapping;
        nucInd = zeros(numel(ind),1); % nucInd is the size of the number of nuclei on the next frame and is used to temporarly store their index in the nuclei structure. E.g. 
        
            % Loop over all mapped nuclei
            
        indToDelete = [];
        for jj = 1:size(mapping{mappingInd},1)
            if all(mapping{mappingInd}(jj,:) == 0) % No match found. Leave it.
                continue;
            end
            if mapping{mappingInd}(jj,2) == 0 % No division. Assign to position on the next frame.
                indXYFrame2 = mapping{mappingInd}(jj,1);
                indNucFrame1 = nucleiIndices(jj);
                
                nuclei(nucleiIndices(jj)).indXY(firstFrameNumber+j-1) = indXYFrame2;
                nuclei(nucleiIndices(jj)).position(firstFrameNumber+j-1,:) = xy{newFrameInd}(indXYFrame2,:);
                %nucInd(jj) = mapping{mappingInd}(nucleiIndices(jj),1); % element jj in 'nuclei' becomes mapping(jj).
                nucInd(mapping{mappingInd}(jj,1)) = indNucFrame1; % element jj in 'nuclei' becomes mapping(jj).
                
            else % The nucleus divided. Assign parent/daughter indices.
                ind1 = numel(nuclei)+1;
                ind2 = numel(nuclei)+2;
                nuclei(ind1).position = nan(totalNumberOfFrames,2);
                nuclei(ind1).indXY = zeros(totalNumberOfFrames,1);
                nuclei(ind2).position = nan(totalNumberOfFrames,2);
                nuclei(ind2).indXY = zeros(totalNumberOfFrames,1);
                nuclei(nucleiIndices(jj)).D = ind1;
                nuclei(nucleiIndices(jj)).E = ind2;
                nuclei(ind1).position(firstFrameNumber+j-1,:) = xy{newFrameInd}(mapping{mappingInd}(jj,1),:);
                nuclei(ind1).indXY(firstFrameNumber+j-1) = mapping{mappingInd}(jj,1);
                nuclei(ind1).P = nucleiIndices(jj);
                nuclei(ind2).position(firstFrameNumber+j-1,:) = xy{newFrameInd}(mapping{mappingInd}(jj,2),:);
                nuclei(ind2).indXY(firstFrameNumber+j-1) = mapping{mappingInd}(jj,2);
                nuclei(ind2).P = nucleiIndices(jj);
                indToDelete = [indToDelete jj];
                nucInd(mapping{mappingInd}(jj,1)) = ind1;
                nucInd(mapping{mappingInd}(jj,2)) = ind2;
 %% =>               nucInd(ind1) = mapping{mappingInd};
            end
            %nuclei(nucleiIndices(jj)).score(startingFrame+j) = score(jj);
        end
        % Add the nuclei that were found in the next frame but that weren't
        % mapped to any current nucleus.
        orphanNuclei = find(ind == 0);
        
        indToDelete = [];
        for jj = 1:numel(orphanNuclei)
            xpos = xy{newFrameInd}(orphanNuclei(jj),1);
            ypos = xy{newFrameInd}(orphanNuclei(jj),2);
            
            %if xpos >= 1+edgeClearance && xpos <= size(img,1)-edgeClearance && ypos >= 1+edgeClearance && ypos <= size(img,2)-edgeClearance
                
                IND = numel(nuclei)+1;
                nuclei(IND).position = nan(totalNumberOfFrames,2);
                nuclei(IND).indXY = zeros(totalNumberOfFrames,1);
                nuclei(IND).P = [];
                nuclei(IND).D = [];
                nuclei(IND).E = [];
                nuclei(IND).position(firstFrameNumber+j-1,:) = xy{newFrameInd}(orphanNuclei(jj),:);
                nuclei(IND).indXY(firstFrameNumber+j-1) = orphanNuclei(jj);
                nucInd(orphanNuclei(jj)) = IND;
        end
        indices = 1:size(xy{newFrameInd},1);
        indices(orphanNuclei(indToDelete)) = [];
        xy{newFrameInd}(orphanNuclei(indToDelete),:) = [];
        nucInd(orphanNuclei(indToDelete)) = [];
        nucleiIndices = nucInd;
%         for jj = 1:numel(nuclei)
%             if nuclei(jj).indXY(firstFrameNumber+j-1) > 0
%                 nuclei(jj).indXY(firstFrameNumber+j-1) = find(nuclei(jj).indXY(firstFrameNumber+j-1) == indices);
%             end
%         end
    end
waitbar((j-1)/(totalNumberOfFrames-1), h_waitbar_tracking, ['Tracking progress : processing frames ' num2str(j-1) ' and ' num2str(j) ' out of ' num2str(totalNumberOfFrames) '...']);
end


end

