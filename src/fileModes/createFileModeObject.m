function FileModeObj = createFileModeObj(fileModeIdentifier)
	if fileModeIdentifier == 'DSPIN'
		FileModeObj = DSPINFileMode();
	elseif fileModeIdentifier == 'TIF'
		FileModeObj = TIFFileMode();
	else
		error('FileMode not supported');
	end
end
