function x = keyPressedCallback(src, event)
	fprintf('Key %s was pressed\n', event.Key); % event.Character also works, review difference later
	x = event.Key;
end
