function [TotalProd,TotalProdError,TotalProdN,...
    MeanTotalProd,SDTotalProd,SETotalProd]=IntegratemRNA(Data,MinParticles,MinEmbryos,varargin)

%Calculates the total amount of mRNA produced per nucleus in each AP bin
%If an extra parameter 'IntegrateAll' is given then it does not distinguish
%between AP positions
%Edited by Yang Joon Kim on 09/2018
%If an extra parameter 'InputOutput' is given, then it will extract the
%Data.Particles field and reconstruct them into structure for input. This
%is because for inputoutput ExperimentType, it makes datastructure
%different to the others (since it has protein channel, which is Data.Nuclei)

IntegrateAll=0;
InputOutput=0;
if ~isempty(varargin)
    if strcmp(varargin{1},'IntegrateAll')
        IntegrateAll=1;
    elseif strcmp(varargin{1},'InputOutput')
        InputOutput = 1;
    end
end

% (optional)Re-construct the data structure for inputoutput datatype
if InputOutput
    for i=1:length(Data)
        DataTemp(i) = Data(i).Particles;
    end
    Data = DataTemp;
end

%Figure out the minimum starting NC
StartNC=[];
for i=1:length(Data)
    StartNC=min([StartNC,Data(i).ncFilterID]);
end


if ~IntegrateAll

    %Calculate the total production for each embryo
    TotalProd=nan(length(Data),length(Data(1).APbinID),14);
    TotalProdError=nan(length(Data),length(Data(1).APbinID),14);
    TotalProdN=nan(length(Data),length(Data(1).APbinID),14);
    for nc=StartNC:14
        for i=1:length(Data)
            %If using old data that doesn't use cell arrays for channels, convert to
            %the new format
            if ~iscell(Data(i).APFilter)
                tempAPFilter = Data(i).APFilter;
                tempCompiledParticles = Data(i).CompiledParticles;
                Data(i).APFilter = {};
                Data(i).CompiledParticles = {};
                Data(i).APFilter{1,1} = tempAPFilter;
                Data(i).CompiledParticles{1,1} = tempCompiledParticles;
            end
            
            for j=1:length(Data(1).APbinID)
                if Data(i).APDivision(nc,j)
                    %Use only AP bins with enough area (those that do not have a
                    %NaN)
                    jMin=min(find(~isnan(Data(i).APbinArea)));
                    jMax=max(find(~isnan(Data(i).APbinArea)));

                    if (j>=jMin)&(j<=jMax)

                        if nc==14
                           FrameRange=Data(i).APDivision(nc,j):length(Data(i).ElapsedTime);
                        else
                           FrameRange=Data(i).APDivision(nc,j):Data(i).APDivision(nc+1,j);
                        end
                        
                        ParticleFilter=Data(i).APFilter{1,1}(:,j)&...
                            Data(i).ncFilter(:,Data(i).ncFilterID==nc);

                        TotalProd(i,j,nc)=...
                            sum([Data(i).CompiledParticles{1,1}(ParticleFilter).TotalmRNA])/...
                            nanmean(Data(i).NEllipsesAP(FrameRange,j));                        
                        TotalProdError(i,j,nc)=...
                            sqrt(sum([Data(i).CompiledParticles{1,1}(ParticleFilter).TotalmRNAError].^2))/...
                            nanmean(Data(i).NEllipsesAP(FrameRange,j));
                        TotalProdN(i,j,nc)=...
                            length([Data(i).CompiledParticles{1,1}(ParticleFilter).TotalmRNA]);
                    end
                end
            end
        end
    end



    %Average over all embryos

    %Average the total amount produced over multiple embryos
    MeanTotalProd=nan(length(Data(1).APbinID),14);
    SDTotalProd=nan(length(Data(1).APbinID),14);
    SETotalProd=nan(length(Data(1).APbinID),14);

    for nc=1:14
        for j=1:length(Data(1).APbinID)
            TotalProdTemp=[];
            ErrorTotalProdTemp=[];
            for i=1:length(Data)
                %Check for the minimum number of particles and the size of the
                %AP window
                if (TotalProdN(i,j,nc)>=MinParticles)&(Data(i).APbinArea(j)>0)
                    TotalProdTemp=[TotalProdTemp,TotalProd(i,j,nc)];
                    ErrorTotalProdTemp=[ErrorTotalProdTemp,TotalProdError(i,j,nc)];
                end
            end


            %Check that we got the minimum number of embryos
            if length(TotalProdTemp)>=MinEmbryos
                %Calculate the mean
                MeanTotalProd(j,nc)=nanmean(TotalProdTemp);
                %Estimate the error from the SD
                SDTotalProd(j,nc)=nanstd(TotalProdTemp);
                %Standard error
                SETotalProd(j,nc)=nanstd(TotalProdTemp)/sqrt(length(TotalProdTemp));
            end
        end
    end
    
    


    %Average over all embryos

    %Average the total amount produced over multiple embryos
    MeanTotalProd=nan(length(Data(1).APbinID),14);
    SDTotalProd=nan(length(Data(1).APbinID),14);
    SETotalProd=nan(length(Data(1).APbinID),14);

    for nc=1:14
        for j=1:length(Data(1).APbinID)
            TotalProdTemp=[];
            ErrorTotalProdTemp=[];
            for i=1:length(Data)
                %Check for the minimum number of particles and the size of the
                %AP window
                if (TotalProdN(i,j,nc)>=MinParticles)&(Data(i).APbinArea(j)>0)
                    TotalProdTemp=[TotalProdTemp,TotalProd(i,j,nc)];
                    ErrorTotalProdTemp=[ErrorTotalProdTemp,TotalProdError(i,j,nc)];
                end
            end


            %Check that we got the minimum number of embryos
            if length(TotalProdTemp)>=MinEmbryos
                %Calculate the mean
                MeanTotalProd(j,nc)=mean(TotalProdTemp);
                %Estimate the error from the SD
                SDTotalProd(j,nc)=std(TotalProdTemp);
                %Standard error
                SETotalProd(j,nc)=std(TotalProdTemp)/sqrt(length(TotalProdTemp));
            end
        end
    end
    
    
else
    
    %Calculate the total production for each embryo
    TotalProd=nan(length(Data),14);
    TotalProdError=nan(length(Data),14);
    TotalProdN=nan(length(Data),14);
    for nc=StartNC:14
        for i=1:length(Data)
            
            %Filter for all particles in the nc and in the AP bins I'm
            %going to use
            jMin=min(find(~isnan(Data(i).APbinArea)));
            jMax=max(find(~isnan(Data(i).APbinArea)));
            
            %Find the frame range for this nc
            if nc==14
               FrameRange=median(Data(i).APDivision(nc,jMin:jMax)):length(Data(i).ElapsedTime);
            else
               FrameRange=Data(i).APDivision(nc,jMin:jMax):Data(i).APDivision(nc+1,jMin:jMax);
            end

            
            
            ParticleFilter=...
                sum(Data(i).APFilter(:,jMin:jMax),2)&Data(i).ncFilter(:,Data(i).ncFilterID==nc);
            
            ParticlesToSum=find(ParticleFilter);
            
            TotalProdTemp=[];
            TotalProdErrorTemp=[];
            for k=1:length(ParticlesToSum)
                TotalProdTemp=[TotalProdTemp,...
                    Data(i).CompiledParticles(ParticlesToSum(k)).TotalmRNA];
                TotalProdErrorTemp=[TotalProdErrorTemp,...
                    Data(i).CompiledParticles(ParticlesToSum(k)).TotalmRNAError];
            end
            
            
            
            
            
            TotalProd(i,nc)=sum(TotalProdTemp)/...
                    mean(sum(Data(i).NEllipsesAP(FrameRange,jMin:jMax),2));
            
            %I'm not sure this is useful
            TotalProdN(i,nc)=nan;
                %length(TotalProdTemp)/(jMax-jMin+1);
            
            TotalProdError(i,nc)=...
                sqrt(sum(TotalProdErrorTemp.^2));
   
        end
    end
    
    
    %Average over all embryos

    %Average the total amount produced over multiple embryos
    MeanTotalProd=nan(14,1);
    SDTotalProd=nan(14,1);
    SETotalProd=nan(14,1);
    
    for nc=StartNC:14
        nanFilter=~isnan(TotalProd(:,nc));
        
        if sum(nanFilter)>=MinEmbryos
    
            %Calculate the mean
            MeanTotalProd(nc)=mean(TotalProd(nanFilter,nc));
            %Estimate the error from the SD
            SDTotalProd(nc)=std(TotalProd(nanFilter,nc));
            %Standard error
            SETotalProd(nc)=std(TotalProd(nanFilter,nc))/sqrt(sum(nanFilter));
        end
    end
    
    

    
    
end
