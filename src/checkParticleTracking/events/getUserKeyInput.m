function cc = getUserKeyInput()
	% Wait for user input to select command to execute
    ct = waitforbuttonpress; % ct == 0 for click and ct == 1 for keypress
    cc = get(Overlay, 'CurrentCharacter');
    
    if strcmpi(cc, '') || ct == 0
        cc = 'donothing';
    end
end