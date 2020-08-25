function qcCheckBoxCallback(hObject,~,checkBoxId,cptState)

  newValue = get(hObject,'Value');
  CC = cptState.CurrentChannelIndex;
  CP = cptState.CurrentParticle;
  FF = cptState.CurrentFrame==cptState.Particles{CC}(CP).Frame;
  varName = cptState.qcFlagFields{checkBoxId};
  
  oldValue = cptState.Particles{CC}(CP).(varName)(FF);
  
  cptState.Particles{CC}(CP).(varName)(FF) = newValue;
  
  % revise approval status as appropriate
  if oldValue == 2
    FrameApproved = cptState.Particles{CC}(CP).FrameApproved;
    FrameApproved(FF) = 1;
    meanUrgentFlags = 0;
    for v = 1:length(cptState.qcFlagFields)
      FrameApproved(FF) = FrameApproved(FF)&&cptState.Particles{CC}(CP).(cptState.qcFlagFields{v})(FF)~=2;
      meanUrgentFlags = meanUrgentFlags + mean(cptState.Particles{CC}(CP).(cptState.qcFlagFields{v})==2);
    end
    cptState.Particles{CC}(CP).FrameApproved = FrameApproved;
    if meanUrgentFlags<0.5
      cptState.Particles{CC}(CP).Approved = 0;
    end
  end
end