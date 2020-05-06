function deltaT_s = getDeltaTFromLasX(time1, time2)

%may have to be adjusted for AM/PM differences. fine for now.

% testingStringF = '04:20:31 PM.872';
% testingString0 = '03:45:07 PM.838';

t0 = datenum(time1, 'HH:MM:SS PM.FFF' );
tf = datenum(time2, 'HH:MM:SS PM.FFF' );

days2seconds = 86400;
deltaT_s = (tf - t0) * days2seconds; 

end