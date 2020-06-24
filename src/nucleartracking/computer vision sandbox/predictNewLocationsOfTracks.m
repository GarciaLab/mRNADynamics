function tracks = predictNewLocationsOfTracks(tracks)

  for k = 1:length(tracks)

      bbox = tracks(k).bbox;

      % Predict the current location of the track.           
      predictedState = predict(tracks(k).kalmanFilter);

      predictedCentroid = predictedState(1:2);

      % Shift the bounding box so that its center is at
      % the predicted location.
      predictedCentroid = int32(predictedCentroid) - bbox(3:4) / 2;
      tracks(k).bbox = [predictedCentroid, bbox(3:4)];

  end
end
