function [Time, Profile, SE] = getMembraneFurrowProfiles(ProfStr)
%%
if  strcmp(lower(ProfStr), 'dubuis')
    
    DataPath = 'S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements/DubuisFigure2BMembraneData.csv';
    DubuisFurrowData = csvread(DataPath,1, 0);
    DubuisTimes = DubuisFurrowData(:,1).';
    DubuisMeanProfile = DubuisFurrowData(:,10).';
    DubuisFixedMeanProfile = DubuisFurrowData(:,11).';
    TF = islocalmin(DubuisMeanProfile);
    DubuisTimes = DubuisTimes(find(TF, 1)+1:end);
    DubuisMeanProfile = DubuisMeanProfile(find(TF, 1)+1:end);
    TF = islocalmin(DubuisMeanProfile);
    DubuisTimes = DubuisTimes(~TF);
    DubuisMeanProfile= DubuisMeanProfile(~TF);
    TF = islocalmax(DubuisMeanProfile);
    DubuisTimes = DubuisTimes(1:find(TF, 1, 'last'));
    DubuisMeanProfile= DubuisMeanProfile(1:find(TF, 1, 'last'));
    DubuisTimes = [0 DubuisTimes];
    DubuisMeanProfile = [0 DubuisMeanProfile];
    Time = DubuisTimes;
    Profile = DubuisMeanProfile;
    SE = NaN(size(Profile));
    
else
    if strcmp(lower(ProfStr), 'yw25csquished_nopv')
        DataPath = 'S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_SquishedywNoPV\MembraneFurrowProfile.mat';
        
    elseif strcmp(lower(ProfStr), 'yw25csquished')
        DataPath ='S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_Squishedyw\MembraneFurrowProfile.mat';
    elseif strcmp(lower(ProfStr), 'hisrfp25c_nopv')
        DataPath = 'S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_AllHisRFPNoPV\MembraneFurrowProfile.mat';
    elseif strcmp(lower(ProfStr), 'hisrfp25c')
        DataPath ='S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_AllHisRFP\MembraneFurrowProfile.mat';
    elseif strcmp(lower(ProfStr), 'hisrfp25csquished_nopv')
        DataPath = 'S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_SquishedHisRFPNoPV\MembraneFurrowProfile.mat';
    elseif strcmp(lower(ProfStr), 'hisrfp25csquished')
        DataPath ='S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_SquishedHisRFP\MembraneFurrowProfile.mat';
    elseif strcmp(lower(ProfStr), 'hisrfp25cunsquished_nopv')
        DataPath ='S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_UnsquishedHisRFPNoPV\MembraneFurrowProfile.mat';
    elseif strcmp(lower(ProfStr), 'hisrfp25cunsquished')
        DataPath = 'S:/Gabriella/Dropbox/FurrowCanalMovieMeasurements\T25C_UnsquishedHisRFP\MembraneFurrowProfile.mat';
    end
    load(DataPath);
    RawProfile = WildTypeProfile.Profile;
    RawTimes = WildTypeProfile.Times;
    RawSE = WildTypeProfile.SE;
    if strcmp(lower(ProfStr), 'hisrfp25c_nopv')
        RawProfile = RawProfile(WildTypeProfile.Counts > 2);
        RawTimes = RawTimes(WildTypeProfile.Counts > 2);
        RawSE = RawSE(WildTypeProfile.Counts > 2);
    elseif strcmp(lower(ProfStr), 'yw25csquished_nopv')
        RawProfile = RawProfile(1:end-2);
        RawTimes = RawTimes(1:end-2);
        RawSE = RawSE(1:end-2);
    elseif strcmp(lower(ProfStr), 'hisrfp25cunsquished_nopv')
        RawProfile = RawProfile(WildTypeProfile.Times < 65 | WildTypeProfile.Counts > 3 );
        RawTimes = RawTimes(WildTypeProfile.Times < 65 | WildTypeProfile.Counts > 3);
        RawSE = RawSE(WildTypeProfile.Times < 65 | WildTypeProfile.Counts > 3);
    end
    
    DiffTimes = diff(RawTimes);
    deltaT = DiffTimes(1);
    sigma = floor(1/deltaT)+1;
    SmoothedProfile = imgaussfilt(RawProfile, sigma);
    SmoothedSE = sqrt(imgaussfilt(RawSE.^2, sigma));
    
    
    
    MaxT = floor(max(round(RawTimes,6)));
    MinT= ceil(min(round(RawTimes,6)));
    InterpolatedTimes = MinT:MaxT;
    
    
    InterpolatedProfile = interp1(RawTimes, SmoothedProfile, InterpolatedTimes);
    InterpolatedSEs = NaN(1, length(InterpolatedProfile));
    
    
    for k = 1:length(InterpolatedTimes)
        t = InterpolatedTimes(k);
        Tindex = find(round(RawTimes,5) == t);
        if isempty(Tindex)
            Tlow = find(RawTimes < t, 1, 'last');
            Thigh = find(RawTimes > t, 1);
            if ~isempty(Tlow) & ~isempty(Thigh)
                t1 = RawTimes(Tlow);
                s1  = SmoothedSE(Tlow);
                t2 = RawTimes(Thigh);
                s2  = SmoothedSE(Thigh);
                dfdt1 = (t2-t)/(t2-t1);
                dfdt2 = (t-t1)/(t2-t1);
                InterpolatedSEs(k) = sqrt(dfdt1^2*s1^2+dfdt2^2*s2^2);
            end
        else
            InterpolatedSEs(k) =  SmoothedSE(Tindex);
        end
    end
    
    TF = ~isnan(InterpolatedProfile);
    InterpolatedProfile  = InterpolatedProfile(TF);
    InterpolatedSEs = InterpolatedSEs(TF);
    InterpolatedTimes = InterpolatedTimes(TF);
    
    TF = islocalmax(InterpolatedProfile) & InterpolatedTimes > max(InterpolatedTimes)*.8;
    if sum(TF) > 0
        InterpolatedTimes = InterpolatedTimes(1:find(TF, 1, 'last'));
        InterpolatedProfile = InterpolatedProfile(1:find(TF, 1, 'last'));
        InterpolatedSEs = InterpolatedSEs(1:find(TF, 1, 'last'));
    end
    
    TF = islocalmin(InterpolatedProfile) & InterpolatedTimes < max(InterpolatedTimes)*.2;
     if sum(TF) > 0 
        InterpolatedTimes = InterpolatedTimes(find(TF, 1):end);
        InterpolatedProfile = InterpolatedProfile(find(TF, 1):end); 
        InterpolatedSEs = InterpolatedSEs(find(TF, 1):end);
    end
    
    TF = islocalmin(InterpolatedProfile) | islocalmax(InterpolatedProfile);
    
    if sum(TF) > 0
        InterpolatedTimes = InterpolatedTimes(~TF);
        InterpolatedProfile = InterpolatedProfile(~TF); 
        InterpolatedSEs = InterpolatedSEs(~TF);
    end
    
    if min(InterpolatedTimes) > MinT 
        P = polyfit(InterpolatedTimes(1:8), InterpolatedProfile(1:8),1);
        
        InterpolatedTimes = [MinT InterpolatedTimes];
        InterpolatedProfile = [P(1)*MinT+P(2) InterpolatedProfile];
        InterpolatedSEs = [NaN InterpolatedSEs];
    end
    
    Time = InterpolatedTimes;
    Profile = InterpolatedProfile;
    SE = InterpolatedSEs;
    
      
    MaxT = floor(max(round(Time,6)));
    MinT= ceil(min(round(Time,6)));
    InterpolatedTimes = MinT:MaxT;
    
    
    InterpolatedProfile = interp1(Time, Profile, InterpolatedTimes);
    InterpolatedSEs = NaN(1, length(InterpolatedProfile));
    
    
    for k = 1:length(InterpolatedTimes)
        t = InterpolatedTimes(k);
        Tindex = find(round(Time,5) == t);
        if isempty(Tindex)
            Tlow = find(Time < t, 1, 'last');
            Thigh = find(Time > t, 1);
            if ~isempty(Tlow) & ~isempty(Thigh)
                t1 = Time(Tlow);
                s1  = SE(Tlow);
                t2 = Time(Thigh);
                s2  = SE(Thigh);
                dfdt1 = (t2-t)/(t2-t1);
                dfdt2 = (t-t1)/(t2-t1);
                InterpolatedSEs(k) = sqrt(dfdt1^2*s1^2+dfdt2^2*s2^2);
            end
        else
            InterpolatedSEs(k) =  SE(Tindex);
        end
    end
    
    
    Time = InterpolatedTimes;
    Profile = InterpolatedProfile;
    SE = InterpolatedSEs;
    
     if min(Time) > MinT 
        P = polyfit(Time(1:8), Time(1:8),1);
        
        Time = [MinT Time];
        Profile = [P(1)*MinT+P(2) Profile];
        SE = [NaN SE];
    end
    
    
end
end
