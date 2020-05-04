function indNuclei = useTargetNumber(targetNumber, nuc1, indNuclei)
%% Check the segmentation if a target number was provided.
if exist('targetNumber','var') && numel(targetNumber) == 1 &&...
    isnumeric(targetNumber) && abs(sum(indNuclei(1:numel(nuc1)))-targetNumber)/targetNumber > 0.25

    % Enforce a certain number of nuclei on the current frame, plus or
    % minus 'perc' percent
    perc = 25;
    thresh = graythresh(mat2gray(nuc1));
    indNuclei = im2bw(mat2gray(nuc1),thresh);
    % Test whether the number of nuclei found is within the given range. If
    % not, then chose the threshold to get as close as possible to the 
    % target number.
    if abs(sum(indNuclei(1:numel(nuc1)))-targetNumber)/targetNumber > perc/100
        indNucleiTmp = nan(21,numel([nuc1;nuc_prev;nuc_next]));
        T = 0:.005:1;
        % test a range of threshold and pick the one that returns the
        % number of nuclei the closest to the target number.
        for j = 1:numel(T)
            indNucleiTmp(j,:) = imbinarize(mat2gray([nuc1;nuc_prev;nuc_next]),T(j));
        end
        indNucleiTmp(:,(numel(nuc1)+1):end) = [];
        [~, IND] =  min(abs(sum(indNucleiTmp,2)-targetNumber));
        indNuclei = logical(indNucleiTmp(IND,:));
    end
end

end