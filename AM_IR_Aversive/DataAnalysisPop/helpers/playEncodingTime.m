
whichClass   = 'Full';

% Data settings
fn = set_paths_directories('','',1);

% AM
datadir_AM = fullfile(fn.figs,'ClassAM','AC',whichClass,'each');
q = load(fullfile(fn.processed,'Units'));
UnitInfo_AM = q.UnitInfo;
clear q

% Speech
datadir_Sp = fullfile(fn.figs,'ClassSpeech','Speech',whichClass,'each');
q = load(fullfile(fn.processed,'UnitsVS'));
UnitInfo_Sp = q.UnitInfo;
clear q

% Load SU data
q=load(fullfile(datadir_AM,'CR_each.mat'));
CReach_AM = q.CR;
clear q
q=load(fullfile(datadir_Sp,'CR_each.mat'));
CReach_Sp = q.CR;
clear q

% Also get CTTS
[CTTS_AM,theseCells_AM] = recallDataParams('AC','PopResp');
[CTTS_Sp,theseCells_Sp] = recallDataParams('Speech','PopResp');

PSTH_AM = permute(mean(CTTS_AM,3,'omitnan'),[1 2 4 3]);
PSTH_Sp = permute(mean(CTTS_Sp,3,'omitnan'),[1 2 4 3]);


nlags = 150;
lagvec = [-nlags:-1 1:nlags];

ET = nan(size(PSTH_AM,1),1);
figure;
hold on


for iUn = 90+(1:10)%size(PSTH_Sp,1)
    
    Data = reshape(PSTH_AM(iUn,:,:),[1 size(PSTH_AM,2)*size(PSTH_AM,3)]);
    
    halfACF = acf(Data',nlags);
    fullACF = [flipud(halfACF); halfACF];
    
    fullACF_Sh=[];
    for ii=1:10
        DataSh = Data(1,randperm(size(Data,2)));
        halfACF_Sh = acf(DataSh',nlags);
        fullACF_Sh = [fullACF_Sh [flipud(halfACF_Sh); halfACF_Sh]];
    end
    
    figure;
    hold on
    plot(lagvec,fullACF,'r','LineWidth',2)
    plot(lagvec,fullACF_Sh,'k')
    ylim([-0.2 1])
    pause(1)
% %     
% %     % Fit gaussian
% %     gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
% %     gaussEqn2 = 'g^2 + p * g * (1/sqrt(4*pi*s^2)) * exp(-(x^2)/(4*s^2))';
% %     
% %     f1 = fit(lagvec',fullACF,gaussEqn2,'Normalize','on','Start', [0.1 2 20]);
% %     
% %     yeval = f1.g^2 + f1.p * f1.g * (1/sqrt(4*pi*f1.s^2)) * exp(-(lagvec.^2)/(4*f1.s^2));
% %     
% %     
% % %     figure;
% % %     % plot(lagvec,fullACF)
% % %     % hold on
% %     plot(f1,lagvec, fullACF);
% % %     hold on
% % %     plot(lagvec,yeval,'g')
% %     
% % %     [r,p] = corr(yeval',fullACF)
% %     
% %     % Alt method
% %     % [m,s] = normfit(fullACF);
% %     
% %     ET(iUn) = 2*f1.s;
% %     
% %     if ET(iUn)>80
% %         figure;
% %         plot(f1,lagvec, fullACF);
% %         title(iUn)
% %     end
% %     
end

figure;
histogram(ET,'BinMethod','integers')



