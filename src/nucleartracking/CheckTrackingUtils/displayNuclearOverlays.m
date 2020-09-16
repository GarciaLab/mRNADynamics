function [OverlayFig, overlayAxes, ImageHandle] = ...
    displayNuclearOverlays(OverlayFig,overlayAxes,...
    ImageHandle, cntState,ImageHisMat)

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
    if isempty(cntState.DisplayRange)
        HisOverlayImageMat=cat(3,mat2gray(hisImage),mat2gray(cntState.MaxImageMat),zeros(size(cntState.MaxImageMat)));
    else
        HisOverlayImageMat=cat(3,mat2gray(hisImage,double(cntState.DisplayRange)),mat2gray(cntState.MaxImageMat),zeros(size(cntState.MaxImageMat)));
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
EllipseHandleGreen = [];


schnitzCellNo_Unchecked=[];

for i=1:cntState.numNuclei()
    if (cntState.schnitzcells(i).Checked ==0) 
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellNo_Unchecked= [schnitzCellNo_Unchecked,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    end
end

EllipseHandle = notEllipseCellCNT(cntState, schnitzCellNo_Unchecked, 'y', 10, overlayAxes);

% Show the ones that have been approved
schnitzCellApproved=[];

for i=1:cntState.numNuclei()
    if (cntState.schnitzcells(i).Approved ==1) & (cntState.schnitzcells(i).Checked ==1) 
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellApproved = [schnitzCellApproved,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    end
end

EllipseHandleBlue = notEllipseCellCNT(cntState, schnitzCellApproved, 'b', 10, overlayAxes);

% Show the ones that have been rejected
schnitzCellRejected=[];

for i=1:cntState.numNuclei()
    if (cntState.schnitzcells(i).Approved ~= 1)  & (cntState.schnitzcells(i).Checked ==1) 
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellRejected = [schnitzCellRejected,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    end
end

EllipseHandleRed = notEllipseCellCNT(cntState, schnitzCellRejected, 'r', 10, overlayAxes);

% Show the corresponding nucleus
if ~isempty(cntState.CurrentNucleus) &&...
    cntState.CurrentNucleus > 0
    if ~isempty(cntState.CurrentNucleusCellNo)
        EllipseHandleGreen = ellipseCellCNT(cntState, cntState.CurrentNucleusCellNo, 'g', 10, overlayAxes);
    end
end







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
    set(OverlayFig,'Name',['Schnitz Cell: ',num2str(cntState.CurrentNucleus),'/',num2str(cntState.numNuclei()),...
        ', Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),...
        ', Z: ',num2str(cntState.CurrentZ),'/',num2str(cntState.ZSlices),' nc: ', num2str(cntState.FrameInfo(cntState.CurrentFrame).nc),...
        ', Flag: ',num2str(FlaggedNuclei(cntState.CurrentNucleus))])
end






end
