clear all
T_vector = fliplr((15:.25:30) + 273.15); % units: K
Ref_index = find(single(T_vector) == single(25+273.15));
EmbryoLength = 554;
x_vector = (0:.005:1)*EmbryoLength; % units: um
Ea_vector = 20:2.5:100; % kJ/mol
Ea_Ref_index = find(single(Ea_vector) == single(60));

R =  8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
NumTemps = length(T_vector);
NumEnergies = length(Ea_vector);
NumXSteps = length(x_vector);

cmap = brewermap(length(T_vector),'Spectral');

D_Ref = 5;%7.4; % units: um^2/s (really 5-10 um^2/s)
Tau_Ref = 3000; % units: s
lambda_Ref = sqrt(D_Ref*Tau_Ref); % units: um
BcdMax_Ref = 60; % nM
Kd_Ref_nm = 0.08*BcdMax_Ref;
Kd_Ref = Kd_Ref_nm*10^(-9);
KdRatio_Ref= exp(5); % units: kJ/mol
r_Ref = 20/60; % units; Pol2/locus/second

%%
D_vector = D_Ref*T_vector./T_vector(Ref_index);% units: um^2/s


tau_mat = NaN(NumTemps, NumEnergies);
r_mat = NaN(NumTemps, NumEnergies);
lambda_mat = NaN(NumTemps, NumEnergies);
BcdMax_mat = NaN(NumTemps, NumEnergies, NumEnergies);
for T_index = 1:NumTemps
    tau_mat(T_index,:) = Tau_Ref*exp((Ea_vector./R)*...
        ((T_vector(Ref_index)-T_vector(T_index))/(T_vector(T_index)*T_vector(Ref_index))));
    r_mat(T_index,:) = r_Ref*exp(-1*(Ea_vector./R)*...
        ((T_vector(Ref_index)-T_vector(T_index))/(T_vector(T_index)*T_vector(Ref_index))));
    lambda_mat(T_index,:) = sqrt(D_vector(T_index)*tau_mat(T_index,:));
    for i = 1:length(Ea_vector)
        for j = 1:length(Ea_vector)
            DeltaEa = Ea_vector(i)-Ea_vector(j)/2;
            BcdMax_mat(T_index,i,j) = BcdMax_Ref*sqrt(D_Ref/D_vector(T_index)*exp(-(DeltaEa/R)*...
                ((T_vector(Ref_index)-T_vector(T_index))/(T_vector(T_index)*T_vector(Ref_index))));
        end
    end
end

%%
BcdProfiles = NaN(NumTemps, NumXSteps, NumEnergies*NumEnergies);
NormedBcdProfiles = NaN(NumTemps, NumXSteps, NumEnergies*NumEnergies);
BcdProfileEnergyIndices = zeros(2,  NumEnergies*NumEnergies, 'uint8');
BcdProfileEnergies = zeros(2,  NumEnergies*NumEnergies, 'single');

for T_index = 1:NumTemps
    counter = 1;
    for Ea_index = 1:NumEnergies
        for Ea_index2 = 1:NumEnergies
            BcdProfileEnergyIndices(1,counter) = Ea_index;
            BcdProfileEnergies(1,counter) = Ea_vector(Ea_index);
            BcdProfileEnergyIndices(2,counter) = Ea_index2;
            BcdProfileEnergies(2,counter) = Ea_vector(Ea_index2);
            BcdProfiles(T_index,:,counter) = BcdMax_mat(T_index,Ea_index,Ea_index2)*exp(-x_vector./lambda_mat(T_index,Ea_index2));
            NormedBcdProfiles(T_index,:,counter) = BcdProfiles(T_index,:,counter)./max(BcdProfiles(T_index,:,counter));
            counter = counter+1;
        end
    end
end

BcdRefProfile = BcdMax_Ref*exp(-x_vector/lambda_Ref);
NormedBcdRefProfile = BcdRefProfile/max(BcdRefProfile);

%%

Kd_vector_M = Kd_Ref.^(T_vector(Ref_index)./T_vector);
Kd_vector = Kd_vector_M*10^9;
KdRatio_vector = KdRatio_Ref.^(T_vector(Ref_index)./T_vector);

NormedHbProfiles = NaN(NumTemps, NumXSteps, NumEnergies*NumEnergies);

for T_index = 1:NumTemps
    for Ea_index = 1:NumEnergies
        for Ea_index2 = 1:NumEnergies
            counter = find(BcdProfileEnergyIndices(1,:) == Ea_index & BcdProfileEnergyIndices(2,:) == Ea_index2);
            Numerator = (1 + BcdProfiles(T_index,:,counter)/Kd_vector(T_index)).^6;
            Denominator = Numerator + KdRatio_vector(T_index);
            NormedHbProfiles(T_index,:,counter) = ...
                Numerator./Denominator;
            counter = counter+1;
        end
    end
end

NormedHbRefProfile = p_on(BcdRefProfile, Kd_Ref_nm, KdRatio_Ref);


%%
HbProfileEnergyIndices = zeros(3,  NumEnergies*NumEnergies*NumEnergies, 'uint8');
HbProfileEnergies = zeros(3,  NumEnergies*NumEnergies*NumEnergies, 'single');

HbProfiles = zeros(NumTemps, NumXSteps, NumEnergies*NumEnergies*NumEnergies,'single' );
for T_index = 1:NumTemps
    counter2 = 1;
    for Ea_index = 1:NumEnergies
        for Ea_index2 = 1:NumEnergies
            counter = find(BcdProfileEnergyIndices(1,:) == Ea_index & BcdProfileEnergyIndices(2,:) == Ea_index2);
            for Ea_index3 = 1:NumEnergies
                HbProfileEnergyIndices(1,counter2) = Ea_index;
                HbProfileEnergyIndices(2,counter2) = Ea_index2;
                HbProfileEnergyIndices(3,counter2) = Ea_index3;
                HbProfileEnergies(1,counter2) = Ea_vector(Ea_index);
                HbProfileEnergies(2,counter2) = Ea_vector(Ea_index2);
                HbProfileEnergies(3,counter2) = Ea_vector(Ea_index3);
                HbProfiles(T_index,:,counter2) = r_mat(T_index,Ea_index3)* NormedHbProfiles(T_index,:,counter);
                counter2 = counter2+1;
            end
        end
    end
end


%%
MinBcdLine = NaN(NumTemps, NumXSteps);
MaxBcdLine = NaN(NumTemps, NumXSteps);
MiddleBcdLine = NaN(NumTemps, NumXSteps);
MinHbLine = NaN(NumTemps, NumXSteps);
MaxHbLine = NaN(NumTemps, NumXSteps);
MiddleHbLine = NaN(NumTemps, NumXSteps);
MinProfIndex= NaN(NumTemps, NumXSteps);
MaxProfIndex= NaN(NumTemps, NumXSteps);
MiddleProfIndex= find((squeeze(BcdProfileEnergies(1,:)) == 60) & (squeeze(BcdProfileEnergies(2,:)) == 60));
for T_index = 1:NumTemps
    for i = 1:size(NormedHbProfiles, 2)
        MinProfIndex(T_index,i) = find(squeeze(NormedHbProfiles(T_index,i,:)) == min(squeeze(NormedHbProfiles(T_index,i,:))), 1);
        MinBcdLine(T_index,i,:) = BcdProfiles(T_index,i,MinProfIndex(T_index,i))/max(BcdProfiles(T_index,:,MinProfIndex(T_index,i)));
        MinHbLine(T_index,i,:) = NormedHbProfiles(T_index,i,MinProfIndex(T_index,i));
        MaxProfIndex(T_index,i) = find(squeeze(NormedHbProfiles(T_index,i,:)) == max(squeeze(NormedHbProfiles(T_index,i,:))), 1);
        MaxBcdLine(T_index,i,:) = BcdProfiles(T_index,i,MaxProfIndex(T_index,i))/max(BcdProfiles(T_index,:,MaxProfIndex(T_index,i)));
        MaxHbLine(T_index,i,:) = NormedHbProfiles(T_index,i,MaxProfIndex(T_index,i));
        MiddleBcdLine(T_index,i,:) = BcdProfiles(T_index,i,MiddleProfIndex)/max(BcdProfiles(T_index,:,MiddleProfIndex));
        MiddleHbLine(T_index,i,:) = NormedHbProfiles(T_index,i,MiddleProfIndex);
        %plot(NormedBcdProfiles(1, :, ),NormedHbProfiles(1,:,i))
    end
end


%%
HbHalfMaxPositions = NaN(NumTemps, NumEnergies*NumEnergies);
for T_index = 1:NumTemps
    for k = 1:NumEnergies*NumEnergies
        ProfVector = squeeze(NormedHbProfiles(T_index,:,k));
        DiffVector = abs(ProfVector-0.5);
        xIndex = find(DiffVector == min(DiffVector), 1);
        HbHalfMaxPositions(T_index, k) = x_vector(xIndex);
    end
end


%%
BcdLambdaPositions = NaN(NumTemps, NumEnergies*NumEnergies);
for T_index = 1:NumTemps
    for k = 1:NumEnergies*NumEnergies
        ProfVector = squeeze(NormedBcdProfiles(T_index,:,k));
        DiffVector = abs(ProfVector-exp(-1));
        xIndex = find(DiffVector == min(DiffVector), 1);
        BcdLambdaPositions(T_index, k) = x_vector(xIndex);
    end
end


%%
% warning('off', 'curvefit:fit:nonDoubleYData')
% Ea_fits = zeros(NumXSteps, NumEnergies^3, 'single');
% se_Ea_fits = zeros(NumXSteps, NumEnergies^3, 'single');
% LogA_fits = zeros(NumXSteps, NumEnergies^3, 'single');
% se_LogA_fits = zeros(NumXSteps, NumEnergies^3, 'single');
% alpha = 0.95;
% t = tinv((1+alpha)/2, 1);
% xfit = -1./(R*T_vector).';
% 
% startParallelPool(12, false, true);
% p = gcp();
% parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix');
% parfevalOnAll(gcp(), @warning, 0, 'off', 'curvefit:fit:nonDoubleYData');
% for TestIndex = 1:NumEnergies^3
%     if rem(TestIndex, 100)== 0
%         disp([num2str(TestIndex),'/',num2str(NumEnergies^3)])
%     end
%     parfor APIndex = 1:NumXSteps
%         
%         y = log(HbProfiles(:,APIndex,TestIndex));
%         
%         
%         [f, gof] = fit( xfit, y, 'poly1');
%         try
%             ci = confint(f);
%         catch
%             ci = NaN(2,2);
%         end
%         
%         
%         Ea_fits(APIndex,TestIndex) = f.p1;
%         ci_Ea = ci(:,1).';
%         se_Ea_fits(APIndex,TestIndex)  = (ci_Ea(2)-ci_Ea(1)) ./ (2*t); % Standard Error
%         
%         LogA_fits(APIndex,TestIndex) = f.p2;
%         ci_LogA = ci(:,2).';
%         se_LogA_fits(APIndex,TestIndex)  = (ci_LogA(2)-ci_LogA(1)) ./ (2*t); % Standard Error
%         
%     end
% end

save('S:/Gabriella/Dropbox\MatlabEnvironments\HbTheoryProfilesWithFits20210506.mat')
%%
% for i = 1:100
% TestIndex = randi([1 NumEnergies^3], 1);
% APIndex = 13;
% close all
% FigHandle = figure(1);
% FigAx = axes(FigHandle);
% set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.6]);
% set(gcf,'color','w');
% 
% plot(1./(R*T_vector), exp(LogA_fits(APIndex, TestIndex))*exp(-Ea_fits(APIndex,TestIndex)./(R*T_vector)), 'b-')
% hold on 
% scatter(1./(R*T_vector), HbProfiles(:,APIndex,TestIndex).', 80, 'r.')
% legend_label = strcat('$E_A = ', sprintf('%.1f ', -1*Ea_fits(APIndex,TestIndex)), ' \pm ', sprintf('%.1f ',  se_Ea_fits(APIndex,TestIndex)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
% legend_label2 = strcat('$\textrm{Pred. Prof.}: E_A^1=', sprintf('%.0f ', HbProfileEnergies(1,TestIndex)),'\,\textrm{kJ}\cdot\textrm{mol}^{-1}, E_A^2=', sprintf('%.0f ', HbProfileEnergies(2,TestIndex)),'\,\textrm{kJ}\cdot\textrm{mol}^{-1}, E_A^3=', sprintf('%.0f ', HbProfileEnergies(3,TestIndex)),'\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
% xlabel('1/RT (mol \cdot kJ^{-1})','Fontsize', 16)
% ylabel('[Hb] (nM)','Fontsize', 16)
% 
% set(FigAx,'Fontsize',16)
% legend(legend_label, legend_label2, 'Interpreter', 'latex', 'Location', 'Northwest', 'Fontsize', 14)
% hold off
% 
% saveas(FigHandle, ['S:/Gabriella/Dropbox\GMPlots\TemperatureAnalysis\Hunchback\Theory\20210502\',...
%     'APindex',num2str(APIndex), 'ValIndex',num2str(TestIndex),'.png'])
% 
% end
% 
% %%
% for i = 1:100
% TestIndex = randi([1 NumEnergies^3], 1);
% APIndex = 13;
% close all
% FigHandle = figure(1);
% FigAx = axes(FigHandle);
% set(FigHandle,'units', 'normalized', 'position',[0.05, 0.05, 0.5, 0.6]);
% set(gcf,'color','w');
% 
% 
% scatter(x_vector/EmbryoLength, abs(-Ea_fits(:,TestIndex)-HbProfileEnergies(3,TestIndex)).', 80, 'r.')
% %legend_label = strcat('$E_A = ', sprintf('%.1f ', -1*Ea_fits(APIndex,TestIndex)), ' \pm ', sprintf('%.1f ',  se_Ea_fits(APIndex,TestIndex)), '\,\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
% legend_label2 = strcat('$\textrm{Pred. Prof.}: E_A^1=', sprintf('%.0f ', HbProfileEnergies(1,TestIndex)),'\,\textrm{kJ}\cdot\textrm{mol}^{-1}, E_A^2=', sprintf('%.0f ', HbProfileEnergies(2,TestIndex)),'\,\textrm{kJ}\cdot\textrm{mol}^{-1}, E_A^3=', sprintf('%.0f ', HbProfileEnergies(3,TestIndex)),'\,\textrm{kJ}\cdot\textrm{mol}^{-1}$');
% ylabel('$\Delta E_{A}\,\, (\textrm{kJ} \cdot \textrm{mol}^{-1})$','Interpreter', 'latex','Fontsize', 16)
% xlabel('x/L','Fontsize', 16)
% ylim([0, max(abs(-Ea_fits(:,TestIndex)-HbProfileEnergies(3,TestIndex)))*1.2])
% set(FigAx,'Fontsize',16)
% legend(legend_label2, 'Interpreter', 'latex', 'Location', 'Northwest', 'Fontsize', 14)
% grid on 
% box on 
% hold off
% 
% saveas(FigHandle, ['S:/Gabriella/Dropbox\GMPlots\TemperatureAnalysis\Hunchback\Theory\20210502\',...
%      'CrossAPValIndex',num2str(TestIndex),'.png'])
% 
% end