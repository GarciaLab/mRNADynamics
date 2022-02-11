function sizesListInPixels = convertSizeToPixels(sizesCellArray, ratio)
	if ratio == 1
		sizesListInPixels = sizesCellArray;
	else
		sizesListInPixels = cellfun(@(x) {round(x / ratio)}, sizesCellArray);
	end
end