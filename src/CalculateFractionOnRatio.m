function [MeanFracOn,SEFracOn,FracOn]=CalculateFractionOnRatio(Data,MinParticles,MinEmbryos)

%Calculate the fraction of ON nuclei given input CompiledParticles
%datasets. This calculates the fraction of ON nuclei given the variables
%EllipsesOnAP and TotalEllipsesAP, so the fraction calculated is using any
%nucleus that turns over over the course of the nuclear cycle.

%Inputs:
%   Data: array of compiled data (e.g. from LoadMS2Sets)
%   MinParticles: minimum number of particles in each AP bin per dataset
%   MinEmbryos: minimum number of embryos to average in each AP bin

%Last updated 7/9/19 by Jonathan Liu.

FracOn=nan(length(Data),length(Data(1).APbinID),3);
for i=1:length(Data)
    nc12 = Data(i).nc12;
    nc13 = Data(i).nc13;
    nc14 = Data(i).nc14;

    %If any nuclear cycle indices are zero, set them to one.
    if nc12 == 0
        nc12 = 1;
    end
    if nc13 == 0
        nc13 = 1;
    end
    if nc14 == 0
        nc14 = 1;
    end
    %We will only go ahead if we have at least a MinParticles
    %number of ellipses to check

    %Max # of particles per AP bin for each nuclear cycle
    NParticlesAPMax = [max(Data(i).NParticlesAP(nc12:nc13,:),[],1);...
        max(Data(i).NParticlesAP(nc13:nc14,:),[],1);...
        max(Data(i).NParticlesAP(nc14:end,:),[],1)];
    NParticlesAPMax = NParticlesAPMax'; %Transpose to keep same dimensions as EllipsesOnAP

    %Calculate fraction on for each AP bin
    FracOnTemp = Data(i).EllipsesOnAP./Data(i).TotalEllipsesAP;


    %Replace AP bins that don't have enough particles with nan
    FracOnTemp(NParticlesAPMax < MinParticles) = nan;

    %Save fraction on
    FracOn(i,:,:) = FracOnTemp;
end

%Average all the embryos
MeanFracOn=nan(length(Data(1).APbinID),3);
SEFracOn=nan(length(Data(1).APbinID),3);
for nc=1:3
    for j=1:length(Data(1).APbinID)
        Temp=FracOn(:,j,nc);
        Temp=Temp(~isnan(Temp)); %Remove datasets with nan in this AP bin and nuclear cycle
        if length(Temp)>=MinEmbryos
            MeanFracOn(j,nc)=mean(Temp);
            SEFracOn(j,nc)=std(Temp)/sqrt(length(Temp));
        end
    end
end

