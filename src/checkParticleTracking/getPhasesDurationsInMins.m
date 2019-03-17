function [anaphaseInMins, prophaseInMins, metaphaseInMins] = getPhasesDurationsInMins(anaphase, prophase, metaphase, ElapsedTime)

  anaphaseInMins = anaphase;

  for i = 1:length(anaphase)
    if anaphase(i) > 0
      anaphaseInMins(i) = ElapsedTime(anaphase(i)); % in units of minutes
    end
  end

  %prophase and metaphase
  prophaseInMins = [];
  metaphaseInMins = [];

  try
    prophaseInMins = prophase;

    for i = 1:length(prophase)
      if prophase(i) > 0
        prophaseInMins(i) = ElapsedTime(prophase(i)); %mins
      end
    end

    metaphaseInMins = metaphase;

    for i = 1:length(metaphase)
      if metaphase(i) > 0
        metaphaseInMins(i) = ElapsedTime(metaphase(i)); %mins
      end
    end

  end
end