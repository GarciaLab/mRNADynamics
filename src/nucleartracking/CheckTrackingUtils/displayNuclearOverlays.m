function [OverlayFig, overlayAxes, ImageHandle] = ...
    displayNuclearOverlays(OverlayFig,overlayAxes,...
    ImageHandle, cntState,ImageHisMat)
hasInputChannels = ~isempty(cntState.liveExperiment.inputChannels);
hisImage = ImageHisMat(:, :, cntState.CurrentFrame);
numFrames = length(cntState.Ellipses);
ApprovedNuclei = [cntState.schnitzcells.Approved];
FlaggedNuclei = [cntState.schnitzcells.Flag];
CheckedNuclei = [cntState.schnitzcells.Checked];
% Define the overlayAxes
if length(overlayAxes.Children) > 1
    delete(overlayAxes.Children(1:end-1))
end
if cntState.UseHistoneOverlay
    if hasInputChannels
        if isempty(cntState.DisplayRange)
            HisOverlayImageMat=cat(3,mat2gray(hisImage),mat2gray(cntState.MaxImageMat),zeros(size(cntState.MaxImageMat)));
        else
            HisOverlayImageMat=cat(3,mat2gray(hisImage,double(cntState.DisplayRange)),mat2gray(cntState.MaxImageMat),zeros(size(cntState.MaxImageMat)));
        end
    else
        if isempty(cntState.DisplayRange)
            HisOverlayImageMat=cat(3,mat2gray(hisImage),zeros(size(hisImage)),zeros(size(hisImage)));
        else
            HisOverlayImageMat=cat(3,mat2gray(hisImage,double(cntState.DisplayRange)),zeros(size(hisImage)),zeros(size(hisImage)));
        end
    end

    ImageHandle.CData = HisOverlayImageMat;
else
    ImageHandle.CData = cntState.MaxImageMat;
    
end

%% 


hold(overlayAxes, 'on')

EllipseHandle = [];
EllipseHandleBlue = [];
EllipseHandleRed = [];
EllipseHandleMagenta = [];
EllipseHandleOrange = [];
EllipseHandleGreen = [];
EllipseHandleWhite= [];

schnitzCellNo_ApprovedUnchecked=[];
schnitzCellNo_ApprovedUncheckedPassed=[];
schnitzCellApproved=[];
schnitzCellCheckedRejected=[];
schnitzCellUncheckedRejected=[];
schnitzCellUncheckedRejectedUnflagged = [];


for i=1:cntState.numNuclei()
    if ~isempty(cntState.CurrentNucleusCellNo) & (i == cntState.CurrentNucleus)
        EllipseHandleGreen = ellipseCellCNT(cntState, cntState.CurrentNucleusCellNo, 'g', 10, overlayAxes);
    elseif (cntState.schnitzcells(i).Checked ==0) & (cntState.schnitzcells(i).Approved ==1) & (cntState.schnitzcells(i).FirstPass == 0)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellNo_ApprovedUnchecked= [schnitzCellNo_ApprovedUnchecked,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Checked ==0) & (cntState.schnitzcells(i).Approved ==1) & (cntState.schnitzcells(i).FirstPass == 1)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellNo_ApprovedUncheckedPassed= [schnitzCellNo_ApprovedUncheckedPassed,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ==1) & (cntState.schnitzcells(i).Checked ==1) 
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellApproved = [schnitzCellApproved,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ~=1) & (cntState.schnitzcells(i).Checked ==1) & (cntState.schnitzcells(i).Flag ~= 0)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellCheckedRejected = [schnitzCellCheckedRejected,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ~=1) & (cntState.schnitzcells(i).Checked ==1)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellUncheckedRejectedUnflagged = [schnitzCellUncheckedRejectedUnflagged,cntState.schnitzcells(i).cellno(frame_idx)];
        end    
    elseif (cntState.schnitzcells(i).Approved ~=1) & (cntState.schnitzcells(i).Checked ==0)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellUncheckedRejected = [schnitzCellUncheckedRejected,cntState.schnitzcells(i).cellno(frame_idx)];
        end    
    end
end

EllipseHandle = notEllipseCellCNT(cntState, schnitzCellNo_ApprovedUnchecked, 'y', 10, overlayAxes);
EllipseHandleMagenta = notEllipseCellCNT(cntState, schnitzCellNo_ApprovedUncheckedPassed, 'm', 10, overlayAxes);
EllipseHandleBlue = notEllipseCellCNT(cntState, schnitzCellApproved, 'c', 10, overlayAxes);
EllipseHandleRed = notEllipseCellCNT(cntState, schnitzCellCheckedRejected, 'r', 10, overlayAxes);
EllipseHandleOrange = notEllipseCellCNT(cntState, schnitzCellUncheckedRejected, [1, 0.5, 0], 10, overlayAxes);
EllipseHandleWhite = notEllipseCellCNT(cntState, schnitzCellUncheckedRejectedUnflagged, [1, 0.5, 0], 10, overlayAxes);



if (ApprovedNuclei(cntState.CurrentNucleus) == 1) 
    set(OverlayFig,'Color','g')
elseif (ApprovedNuclei(cntState.CurrentNucleus) <= 0)
    set(OverlayFig,'Color','r')
elseif (ApprovedNuclei(cntState.CurrentNucleus) == 2) & (CheckedNuclei(cntState.CurrentNucleus) == 1)
    set(OverlayFig,'Color','y')
else
    set(OverlayFig,'Color','default')
end


if isfield(cntState.FrameInfo, 'nc')
    if isfield(cntState.schnitzcells, 'anaphaseFrame')
        set(OverlayFig,'Name',['Schnitz Cell: ',num2str(cntState.CurrentNucleus),'/',num2str(cntState.numNuclei()),...
            ', Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),...
            ', Ana: ',num2str(cntState.schnitzcells(cntState.CurrentNucleus).anaphaseFrame),' nc: ', num2str(cntState.schnitzcells(cntState.CurrentNucleus).cycle),...
            ', Flag: ',num2str(FlaggedNuclei(cntState.CurrentNucleus))])

    else
        set(OverlayFig,'Name',['Schnitz Cell: ',num2str(cntState.CurrentNucleus),'/',num2str(cntState.numNuclei()),...
            ', Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),...
            ', Z: ',num2str(cntState.CurrentZ),'/',num2str(cntState.ZSlices),' nc: ', num2str(cntState.schnitzcells(cntState.CurrentNucleus).cycle),...
            ', Flag: ',num2str(FlaggedNuclei(cntState.CurrentNucleus))])
    end
end






end
