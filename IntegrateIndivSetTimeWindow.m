function [Integral,SDIntegral,IntegralON,SDIntegralON]=IntegrateIndivSetTimeWindow(Data,MinParticles,...
    nc,TimeWindow,varargin)

%Integrates a particular data set over a specific time window in an nc.
%TimeWindow is a vector with the start and end time point measured with
%respect to the average start of the nc. The average start of the nc is
%obtained from APDivisions

%Integral: Integrated amount divided by the total number of nuclei in each
%AP window.
%IntegralON: Mean integrated amount of ON nuclei.


IntegrateAll=0;
if ~isempty(varargin)
    if strcmp(varargin{1},'IntegrateAll')
        IntegrateAll=1;
    end
    if length(varargin)==2
        Filter=varargin{2};
    end
end


if ~IntegrateAll
    Integral=nan(size(Data.APbinID));
    SDIntegral=nan(size(Data.APbinID));
    IntegralON=nan(size(Data.APbinID));
    SDIntegralON=nan(size(Data.APbinID));
    for j=1:length(Data.APbinID)

        if Data.APbinArea(j)>0

            %Find the frame we should start and end with for this AP window
            [Dummy,IndexStart]=...
                min((Data.ElapsedTime-(Data.ElapsedTime(Data.APDivision(nc,j))+TimeWindow(1))).^2);
            [Dummy,IndexEnd]=...
                min((Data.ElapsedTime-(Data.ElapsedTime(Data.APDivision(nc,j))+TimeWindow(2))).^2);


            ParticlesToSum=find((Data.ncFilter(:,Data.ncFilterID==nc))&Data.APFilter(:,j));

            if ~isempty(ParticlesToSum)
                clear IntegralTemp
                clear SDIntegralTemp
                for i=1:length(ParticlesToSum)
                    %Determine the AP bin of the particle so we can determine the time
                    %range with respect to this bin's division.
                    [IntegralTemp(i),SDIntegralTemp(i)]=...
                        IntegrateParticle(Data.CompiledParticles(ParticlesToSum(i)),...
                        [IndexStart:IndexEnd],Data.ElapsedTime);

                    %To check discrepancies related to particles not being assigned
                    %to the right nc
                    [IntegralTemp(i),Data.CompiledParticles(ParticlesToSum(i)).TotalmRNA];
                end

                %Check that we have enough particles
                if sum(IntegralTemp>0)>=MinParticles



                    %Add up the integrated amounts and divide by the number of nuclei
                    if nc==14
                        ncRange=Data.APDivision(nc,j):length(Data.ElapsedTime);
                    else
                        ncRange=Data.APDivision(nc,j):Data.APDivision(nc+1,j);
                    end

                    %To calculate the amount produced per ALL nuclei I average the
                    %production of all detected spots and then mulitply by the
                    %fraction of nuclei with detected expression
                    Integral(j)=mean(IntegralTemp(IntegralTemp~=0))*Data.EllipsesOnAP(j,nc-11)/Data.TotalEllipsesAP(j,nc-11);
                    SDIntegral(j)=sqrt(mean(SDIntegralTemp(IntegralTemp~=0).^2))*Data.EllipsesOnAP(j,nc-11)/Data.TotalEllipsesAP(j,nc-11);

                    %Note that I'm discarding those particles that give me zero
                    %integral
                    IntegralON(j)=mean(IntegralTemp(IntegralTemp~=0));
                    SDIntegralON(j)=sqrt(mean(SDIntegralTemp(IntegralTemp~=0).^2));
                end
            else
                Integral(j)=nan;
                SDIntegral(j)=nan;

                IntegralON(j)=nan;
                SDIntegralON(j)=nan;
            end

        end
    end
else

    %Find the frame we should start and end with for this AP window
    [Dummy,IndexStart]=...
        min((Data.ElapsedTime-(Data.ElapsedTime(Data.nc13)+TimeWindow(1))).^2);
    [Dummy,IndexEnd]=...
        min((Data.ElapsedTime-(Data.ElapsedTime(Data.nc13)+TimeWindow(2))).^2);


    ParticlesToSum=find((Data.ncFilter(:,Data.ncFilterID==nc))&Filter);

    
    if ~isempty(ParticlesToSum)
        clear IntegralTemp
        clear SDIntegralTemp
        for i=1:length(ParticlesToSum)
            %Determine the AP bin of the particle so we can determine the time
            %range with respect to this bin's division.
            [IntegralTemp(i),SDIntegralTemp(i)]=...
                IntegrateParticle(Data.CompiledParticles(ParticlesToSum(i)),...
                [IndexStart:IndexEnd],Data.ElapsedTime);

            %To check discrepancies related to particles not being assigned
            %to the right nc
            [IntegralTemp(i),Data.CompiledParticles(ParticlesToSum(i)).TotalmRNA];
        end

        %Check that we have enough particles
        if sum(IntegralTemp>0)>=MinParticles

            %Add up the integrated amounts and divide by the number of nuclei
            if nc==14
                ncRange=Data.nc14:length(Data.ElapsedTime);
            else
                ncRange=eval(['Data.nc',num2str(nc),':Data.nc',num2str(nc+1)]);
            end

            %To calculate the amount produced per ALL nuclei I average the
            %production of all detected spots and then mulitply by the
            %fraction of nuclei with detected expression
            Integral(j)=mean(IntegralTemp(IntegralTemp~=0))*Data.EllipsesOnAP(j,nc-11)/Data.TotalEllipsesAP(j,nc-11);
            SDIntegral(j)=sqrt(mean(SDIntegralTemp(IntegralTemp~=0).^2))*Data.EllipsesOnAP(j,nc-11)/Data.TotalEllipsesAP(j,nc-11);

            %Note that I'm discarding those particles that give me zero
            %integral
            IntegralON(j)=mean(IntegralTemp(IntegralTemp~=0));
            SDIntegralON(j)=sqrt(mean(SDIntegralTemp(IntegralTemp~=0).^2));
        end
    else
        Integral(j)=nan;
        SDIntegral(j)=nan;

        IntegralON(j)=nan;
        SDIntegralON(j)=nan;
    end

end

