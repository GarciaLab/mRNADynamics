function schnitzcells = AppendSchnitzcellsData(schnitzcells,row1,row2)

%pass fields' info by concatenating them to the
%outer loop schnitz ones

schnitzcells(row1).StitchedTo = [schnitzcells(row1).StitchedTo schnitzcells(row2).StitchedTo];
schnitzcells(row1).StitchedFrom = [schnitzcells(row1).StitchedFrom schnitzcells(row2).StitchedFrom];
schnitzcells(row1).frames = [schnitzcells(row1).frames ; schnitzcells(row2).frames];
schnitzcells(row1).cenx = [schnitzcells(row1).cenx schnitzcells(row2).cenx];
schnitzcells(row1).ceny = [schnitzcells(row1).ceny schnitzcells(row2).ceny];

if isfield(schnitzcells, 'len')
    schnitzcells(row1).len = [schnitzcells(row1).len schnitzcells(row2).len];
end

schnitzcells(row1).cellno = [schnitzcells(row1).cellno schnitzcells(row2).cellno];

if isfield(schnitzcells,'Fluo')
    schnitzcells(row1).Fluo = [schnitzcells(row1).Fluo;schnitzcells(row2).Fluo]; %
end
if isfield(schnitzcells, 'APpos')
    schnitzcells(row1).APpos = [schnitzcells(row1).APpos;schnitzcells(row2).APpos];
end
if isfield(schnitzcells, 'APpos')
    schnitzcells(row1).DVpos =  [schnitzcells(row1).DVpos;schnitzcells(row2).DVpos];
end
if isfield(schnitzcells, 'FrameApproved')
    schnitzcells(row1).FrameApproved =  [schnitzcells(row1).FrameApproved,schnitzcells(row2).FrameApproved];
end
if isfield(schnitzcells, 'Approved')
    schnitzcells(row1).Approved =  [schnitzcells(row1).Approved || schnitzcells(row2).Approved];
end
if isfield(schnitzcells, 'FluoTimeTrace')
    schnitzcells(row1).FluoTimeTrace =  [schnitzcells(row1).FluoTimeTrace;schnitzcells(row2).FluoTimeTrace];
end

if isfield(schnitzcells,'ExtendedIntoFutureWithThisThresh')
    schnitzcells(row1).ExtendedIntoFutureWithThisThresh = schnitzcells(row1).ExtendedIntoFutureWithThisThresh;
end

end













%
% F = getfield(S,FIELD) returns the contents of the specified field. For
%     example, if S.a = 1, then getfield(S,'a') returns 1. FIELD can be a
%     character vector or string scalar.
%
%         NAMES = fieldnames(S) returns a cell array of character vectors
%     containing the names of the fields in structure S.