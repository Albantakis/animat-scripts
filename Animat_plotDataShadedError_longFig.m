clear all
option = 2;
FA =  [0 128];
condNum = 3;

%% Zombie Data
%condition = 'c1a2_change_c23a14';
%condition = 'c3a4_36';
condition = 'c36a45_36'
%load(strcat(condition,'_ZombiedataAllC'));
load(strcat(condition,'_dataCB'));
fitnessBeforeChange = max(Fitness_level(:,59),[],2);

if option == 1
    %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
    %& (fitnessBeforeChange == 122)
    indF = find(max(BigPhiMip(:,55:59),[],2) > 0 & (fitnessBeforeChange >= 126));% & max(BigPhiMip(:,60:64),[],2) == 0) % find(mean(Fitness_level(:,end-10:end),2) > FA(i) & mean(Fitness_level(:,end-10:end),2) <= FA(i+1))
    %Those with 2 nodes.
%     indHU = [1     3     4     5     7     9    10    14    15    16    17    18    23    25    28    29    30    31    33    36    39    41    46    50    51    52    53    55    56    62];
%     indF = intersect(indF, indHU);
    gcolor = [0 0 0];
elseif option == 4
    %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
    indF = find(max(BigPhiMip(:,55:59),[],2) == 0 & (fitnessBeforeChange >= 126));% & max(BigPhiMip(:,60:64),[],2) == 0) % find(mean(Fitness_level(:,end-10:end),2) > FA(i) & mean(Fitness_level(:,end-10:end),2) <= FA(i+1))
    %Those with more than 2 nodes.
%     indHU = [2     6     8    11    12    13    19    20    21    22    24    26    27    32    34    35    37    38    40    42    43    44    45    47    48    49    54    57    58    59    60    61];
%     indF = intersect(indF, indHU);
    gcolor = [0 0 1];
elseif option == 2
    comp = 120.4760; %mean fitness of last 10 Fitness_level of c1a3
    BestFit = 128;
    Fit = mean(Fitness_level(:,end-9:end),2);
    [FitB, indFit] = sort(Fit, 'descend')
    for k = 1:50
        meanF = mean(FitB(1:k));
        %if abs(meanF-comp) < abs(BestFit-comp)
        if meanF < comp
            BestFit = meanF;
            indF = indFit(1:k);
            break;
        end
    end
    if numel(indF) == 1
        indF = indFit(1:10); %for Motor best
        gcolor = [0 0.5 0]; %oneMot best
    else
        %gcolor = [1 .5 0]; %oneSen best
        gcolor = [1 0 0];
    end
elseif option == 3
    load SimulatedFitnessNoise
    Fit = SimulatedFitness;
    [FitB, indFit] = sort(Fit, 'descend')
    indF = indFit(1:20);
    gcolor = [0.5,0,0.5]; %noise best
else
    indF = [1:size(Fitness_level, 1)];
    %gcolor = [0.75 0.75 0.5];
    %gcolor = [1/i,1-1/i,0];
    %gcolor = [1 .7 0]; %oneSen all
    %gcolor = [0 0.7 0]; %oneMot all
    %gcolor = [1,0,1]; %noise all
    gcolor = [ 0 0 0.5];
end

totsteps = 60000-1;


%MindHF = median(mean(Fitness_level(:,end-10:end),2)); 
%indHF = find(mean(Fitness_level(:,end-10:end),2) > MindHF);
%indLF = setdiff(1:50, indHF)';

length(indF)
inc = indF %[1:50];
MaxFitness = 128;
MeanFitness = mean(Fitness_level(inc,:),1);
MeanFitness = 100.*MeanFitness./MaxFitness;
Fitness_level = 100.*Fitness_level./MaxFitness;
%range = floor(range./1.5);

%Fitness_level(:,118)./Fitness_level(:,59)
%% Correlation Histograms
Zombie = load(strcat(condition,'_ZombiedataAllC'));
load(strcat(condition, '_dataCB'));

%range1S = 1:49;
%mean correlation coefficients
F = mean(Fitness_level(inc,:));
P = mean(BigPhiMip(inc,:));
EC = mean(MeanSizeComplex(inc,:));
C = mean(MeanNumConcepts(inc,:));
ZP = mean(Zombie.BigPhi(inc,:));
ZC = mean(Zombie.MeanNumConcepts(inc,:));
%CapME = mean(results.Capture_maxEnt);
%CapInd = mean(results.Capture_Gen1st);

%CorrMat = corrcoef([F; P; C; ZP; ZC; CapME; CapInd]');
[CorrMat, pCorr] = corr([F]',[ZC; ZP; EC; C; P]','type','Spearman');
CorrMat
pCorr

for i = 1:numel(inc)
    F = Fitness_level(inc(i),:);
    P = BigPhiMip(inc(i),:);
    EC = MeanSizeComplex(inc(i),:);
    C = MeanNumConcepts(inc(i),:);
    ZP = Zombie.BigPhi(inc(i),:);
    ZC = Zombie.MeanNumConcepts(inc(i),:);
    %CapME = results.Capture_Gen1st(i,:);
    %CapInd = results.Capture_Change(i,:);
    %CorrMat_trial = corrcoef([F; P; C; ZP; ZC; CapME; CapInd]');
    
    [CorrMat_trial, pCorrMat_trial] = corr([F]',[ZC; ZP; EC; C; P]','type','Spearman');
    CorrMat_trial(isnan(CorrMat_trial)) = 0;
    FitCorr(i,:) = CorrMat_trial;
    pFitCorr(i,:) = pCorrMat_trial;
    %BPhiCorr(i,:) = CorrMat_trial(4,:);
    
    [PhiCorrMat_trial, pPhiCorrMat_trial] = corr([P]',[F; ZC; ZP; EC; C]','type','Spearman');
    PhiCorrMat_trial(isnan(PhiCorrMat_trial)) = 0;
    PhiFitCorr(i,:) = PhiCorrMat_trial;
    pPhiFitCorr(i,:) = pPhiCorrMat_trial;

end

A = mean(FitCorr);
B = std(FitCorr)./sqrt(length(inc));
AB = [A; B]

C = mean(PhiFitCorr);
D = std(PhiFitCorr)./sqrt(length(inc));

CD = [C; D]

for k = 1:5
    indk = find(pFitCorr(:,k) <= 0.05 & FitCorr(:,k) >= 0);
    numCorr(k) = numel(indk);
    
    indPhik = find(pPhiFitCorr(:,k) <= 0.05 & PhiFitCorr(:,k) >= 0);
    numPhiCorr(k) = numel(indPhik);
end
numCorr
numPhiCorr
%%

figure(1000)
s(1) = subplot(5,4,1+condNum);
    hold on
    hist(FitCorr(:,1), -1:0.1:1);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','none','EdgeColor','k')
    hist(FitCorr(pFitCorr(:,1) <= 0.05,1), -1:0.1:1);
    h2 = findobj(gca,'Type','patch');
    set(h2,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(FitCorr(:,1)),mean(FitCorr(:,1))],[0,40]);
    xlim([-1.05,1.05])

s(2) = subplot(5,4,5+condNum)   ; 
    hold on
    hist(FitCorr(:,2), -1:0.1:1);
    h3 = findobj(gca,'Type','patch');
    set(h3,'FaceColor','none','EdgeColor','k')
    hist(FitCorr(pFitCorr(:,2) <= 0.05,2), -1:0.1:1);
    h4 = findobj(gca,'Type','patch');
    set(h4,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(FitCorr(:,2)),mean(FitCorr(:,2))],[0,40]);
    xlim([-1.05,1.05])
    
s(3) = subplot(5,4,9+condNum);
   hold on
    hist(FitCorr(:,3), -1:0.1:1);
    h5 = findobj(gca,'Type','patch');
    set(h5,'FaceColor','none','EdgeColor','k')
    hist(FitCorr(pFitCorr(:,3) <= 0.05,3), -1:0.1:1);
    h6 = findobj(gca,'Type','patch');
    set(h6,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(FitCorr(:,3)),mean(FitCorr(:,3))],[0,40]);
    xlim([-1.05,1.05])

s(4) = subplot(5,4,13+condNum);
     hold on
    hist(FitCorr(:,4), -1:0.1:1);
    h7 = findobj(gca,'Type','patch');
    set(h7,'FaceColor','none','EdgeColor','k')
    hist(FitCorr(pFitCorr(:,4) <= 0.05,4), -1:0.1:1);
    h8 = findobj(gca,'Type','patch');
    set(h8,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(FitCorr(:,4)),mean(FitCorr(:,4))],[0,40]);
    xlim([-1.05,1.05])

s(5) = subplot(5,4,17+condNum);
    hold on
    hist(FitCorr(:,5), -1:0.1:1);
    h9 = findobj(gca,'Type','patch');
    set(h9,'FaceColor','none','EdgeColor','k')
    hist(FitCorr(pFitCorr(:,5) <= 0.05,5), -1:0.1:1);
    h10 = findobj(gca,'Type','patch');
    set(h10,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(FitCorr(:,5)),mean(FitCorr(:,5))],[0,40]);
    xlim([-1.05,1.05])  
    
% s(6) = subplot(4,2,7);
%     hold on
%     hist(FitCorr(:,7), -1:0.1:1);
%     line([mean(FitCorr(:,7)),mean(FitCorr(:,7))],[0,10]);
%     xlim([-1.05,1.05])
% s(7) = subplot(4,2,6);
%     hold on
%     hist(BPhiCorr(:,6), -1:0.1:1);
%     line([mean(BPhiCorr(:,6)),mean(BPhiCorr(:,6))],[0,10]); 
%     xlim([-1.05,1.05])
% s(8) = subplot(4,2,8);
%     hold on
%     hist(BPhiCorr(:,7), -1:0.1:1);
%     line([mean(BPhiCorr(:,7)),mean(BPhiCorr(:,7))],[0,10]);
%     xlim([-1.05,1.05])
    
title(s(5),'$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12)
title(s(2),'$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12)
title(s(4),'#concepts')
title(s(1),'#concepts')
title(s(3),'#elements')

%% Correlations with Phi
figure(1001)
s(1) = subplot(5,4,1+condNum);
    hold on
    hist(PhiFitCorr(:,1), -1:0.1:1);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','none','EdgeColor','k')
    hist(PhiFitCorr(pPhiFitCorr(:,1) <= 0.05,1), -1:0.1:1);
    h2 = findobj(gca,'Type','patch');
    set(h2,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(PhiFitCorr(:,1)),mean(PhiFitCorr(:,1))],[0,40]);
    xlim([-1.05,1.05])

s(2) = subplot(5,4,5+condNum)   ; 
    hold on
    hist(PhiFitCorr(:,2), -1:0.1:1);
    h3 = findobj(gca,'Type','patch');
    set(h3,'FaceColor','none','EdgeColor','k')
    hist(PhiFitCorr(pPhiFitCorr(:,2) <= 0.05,2), -1:0.1:1);
    h4 = findobj(gca,'Type','patch');
    set(h4,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(PhiFitCorr(:,2)),mean(PhiFitCorr(:,2))],[0,40]);
    xlim([-1.05,1.05])
    
s(3) = subplot(5,4,9+condNum);
   hold on
    hist(PhiFitCorr(:,3), -1:0.1:1);
    h5 = findobj(gca,'Type','patch');
    set(h5,'FaceColor','none','EdgeColor','k')
    hist(PhiFitCorr(pPhiFitCorr(:,3) <= 0.05,3), -1:0.1:1);
    h6 = findobj(gca,'Type','patch');
    set(h6,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(PhiFitCorr(:,3)),mean(PhiFitCorr(:,3))],[0,40]);
    xlim([-1.05,1.05])

s(4) = subplot(5,4,13+condNum);
     hold on
    hist(PhiFitCorr(:,4), -1:0.1:1);
    h7 = findobj(gca,'Type','patch');
    set(h7,'FaceColor','none','EdgeColor','k')
    hist(PhiFitCorr(pPhiFitCorr(:,4) <= 0.05,4), -1:0.1:1);
    h8 = findobj(gca,'Type','patch');
    set(h8,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(PhiFitCorr(:,4)),mean(PhiFitCorr(:,4))],[0,40]);
    xlim([-1.05,1.05])

s(5) = subplot(5,4,17+condNum);
    hold on
    hist(PhiFitCorr(:,5), -1:0.1:1);
    h9 = findobj(gca,'Type','patch');
    set(h9,'FaceColor','none','EdgeColor','k')
    hist(PhiFitCorr(pPhiFitCorr(:,5) <= 0.05,5), -1:0.1:1);
    h10 = findobj(gca,'Type','patch');
    set(h10,'FaceColor',gcolor,'EdgeColor','k')
    line([mean(PhiFitCorr(:,5)),mean(PhiFitCorr(:,5))],[0,40]);
    xlim([-1.05,1.05])  
    
% s(6) = subplot(4,2,7);
%     hold on
%     hist(FitCorr(:,7), -1:0.1:1);
%     line([mean(FitCorr(:,7)),mean(FitCorr(:,7))],[0,10]);
%     xlim([-1.05,1.05])
% s(7) = subplot(4,2,6);
%     hold on
%     hist(BPhiCorr(:,6), -1:0.1:1);
%     line([mean(BPhiCorr(:,6)),mean(BPhiCorr(:,6))],[0,10]); 
%     xlim([-1.05,1.05])
% s(8) = subplot(4,2,8);
%     hold on
%     hist(BPhiCorr(:,7), -1:0.1:1);
%     line([mean(BPhiCorr(:,7)),mean(BPhiCorr(:,7))],[0,10]);
%     xlim([-1.05,1.05])
    
title(s(1),'Fitness')
title(s(3),'$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12)
title(s(5),'#concepts')
title(s(2),'#concepts')
title(s(4),'#elements')

%%
load(strcat(condition,'_ZombiedataAllC'));
indError = 1:length(range);
rangeA = range(indError)./10000;
figure(12)
subplot(6,4,1+condNum)
  hold on
  if numel(inc) == 1
    FitSEM = zeros(1, length(Fitness_level));
  else  
    FitSEM = std(Fitness_level(inc,:),1)/sqrt(length(inc));
  end  
    shadedErrorBar(rangeA, mean(Fitness_level(inc,indError),1)./(MaxFitness/100), FitSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0);
    xlim([0, max(rangeA)])
ylim([40,100])
    ylabel('Fitness (%)');

% subplot(4,2,3)
%     hold on
%     ConnecSEM = std(Num_Conn(inc,:))/sqrt(length(inc));
%     shadedErrorBar(range(indError), mean(Num_Conn(inc,(indError))),ConnecSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1.5},1)
%     xlim([1, max(range)])
% ylabel('#connections');
% xlabel('#Generations');

subplot(6,4,9+condNum)
    hold on
    if numel(inc) == 1
    BPSEM = zeros(1, length(Fitness_level));
    else  
    BPSEM = std(BigPhi(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(range)])
    shadedErrorBar(range(indError), mean(BigPhi(inc,indError),1),BPSEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
ylabel('$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12);
%xlabel('#Generations');
ylim([0,1])

subplot(6,4,5+condNum)
    hold on
    if numel(inc) == 1
        ConSEM = zeros(1, length(Fitness_level));
    else  
        ConSEM = std(MeanNumConcepts(inc,:))/sqrt(length(inc));
    end
    shadedErrorBar(range(indError)./10000, mean(MeanNumConcepts(inc,(indError)),1),ConSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
        
    HOConSEM = std(MeanHOConcepts(inc,:))/sqrt(length(inc));
    condition2 = strcat('HO',condition);
    %shadedErrorBar(range(indError), mean(MeanHOConcepts(inc,(indError))),HOConSEM(indError), {'color', gcolor,'DisplayName',condition2,'LineWidth', 1.5},1)  
    xlim([0, max(range)./10000])
ylim([0,4])
ylabel('#concepts');

    %%
load(strcat(condition, '_dataCB'));

indError = 1:length(range);


subplot(6,4,13+condNum)
  hold on
  if numel(inc) == 1
    ComplexSEM = zeros(1, length(Fitness_level));
  else  
    ComplexSEM = std(MeanSizeComplex(inc,:),1)/sqrt(length(inc));
  end  
    shadedErrorBar(range(indError)./10000, mean(MeanSizeComplex(inc, indError),1), ComplexSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    xlim([0, max(range)./10000])
ylabel('#elements');
ylim([0,1.5])

% subplot(4,2,4)
%     hold on
%    HOConSEM = std(MeanHOConcepts(inc,:))/sqrt(length(inc));
%     shadedErrorBar(range(indError), mean(MeanHOConcepts(inc,(indError))),HOConSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1.5},1)
%     xlim([1, max(range)])
% xlabel('#Generations');
% ylabel('#HO Concepts');
% ylim([0,0.3])

subplot(6,4,21+condNum)
  if numel(inc) == 1
    BPSEM = zeros(1, length(Fitness_level));
  else  
    BPSEM = std(BigPhiMip(inc,:))/sqrt(length(inc));
  end  
    xlim([0, max(range)])
    hold on
    shadedErrorBar(range(indError), mean(BigPhiMip(inc,(indError)),1),BPSEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    plot(range, mean(BigPhiMip_max), '-r')

ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);
xlabel('#Generations');
ylim([0,0.5])

subplot(6,4,17+condNum)
    hold on
  if numel(inc) == 1
    ConSEM = zeros(1, length(Fitness_level));
  else  
    ConSEM = std(MeanNumConcepts(inc,:))/sqrt(length(inc));
  end  
    shadedErrorBar(range(indError)./10000, mean(MeanNumConcepts(inc,(indError)),1),ConSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
ylabel('#concepts');
xlim([0, max(range)./10000])
ylim([0,1.5])


optimalConcepts = MeanNumConcepts./(2.^MeanSizeComplex-1);
%optimalConcepts(isnan(optimalConcepts)) = 0;
figure(103)
hold on
subplot(1,4,1+condNum)
 hold on
  if numel(inc) == 1
    optConSEM = zeros(1, length(Fitness_level));
  else  
    optConSEM = nanstd(optimalConcepts(inc,:))/sqrt(length(inc));
  end  
    shadedErrorBar(range(indError)./10000, nanmean(optimalConcepts(inc,(indError)),1),optConSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
ylabel('optimalConcepts');
xlim([0, max(range)./10000])
ylim([0,1])
%% I_SMMI and I_Pred

load(strcat(condition,'_ISenMot'));
I_SMMI = results.ISenMot;
HMot = results.HMot;
HSen = results.HSen;

% load(strcat(condition,'_ISenMot_2ts'));
% I_SMMI_2ts = results.ISenMot;
load(strcat(condition,'_IPred'));
I_Pred = results.IPred;
H_States = results.HStates;
% load(strcat(condition,'_IPred_2ts'));
% I_Pred_2ts = results.IPred;

figure(104)
subplot(6,4,1+condNum)
  hold on
  if numel(inc) == 1
    FitSEM = zeros(1, length(Fitness_level));
  else  
    FitSEM = std(Fitness_level(inc,:),1)/sqrt(length(inc));
  end  
    shadedErrorBar(rangeA, mean(Fitness_level(inc,indError),1)./(MaxFitness/100), FitSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0);
    xlim([0, max(rangeA)])
ylim([40,100])
    ylabel('Fitness (%)');

subplot(6,4,5+condNum)
    hold on
    if numel(inc) == 1
    ISMMI_SEM = zeros(1, length(Fitness_level));
    else  
    ISMMI_SEM = std(I_SMMI(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(rangeA)])
    shadedErrorBar(rangeA(indError), mean(I_SMMI(inc,indError),1),ISMMI_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    ylabel('I_{SMMI}')
    %xlabel('#Generations');
    %ylim([0,1])
    
subplot(6,4,9+condNum)
    hold on
    if numel(inc) == 1
    HSen_SEM = zeros(1, length(Fitness_level));
    else  
    HSen_SEM = std(HSen(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(rangeA)])
    shadedErrorBar(rangeA(indError), mean(HSen(inc,indError),1),HSen_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    ylabel('H_{Sen}')
    %xlabel('#Generations');
    %ylim([0,1])
    
subplot(6,4,13+condNum)
    hold on
    if numel(inc) == 1
    HMot_SEM = zeros(1, length(Fitness_level));
    else  
    HMot_SEM = std(HMot(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(rangeA)])
    shadedErrorBar(rangeA(indError), mean(HMot(inc,indError),1),HMot_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    
    ylabel('H_{Mot}')
    %xlabel('#Generations');
    %ylim([0,1])
    
% subplot(5,4,9+condNum)
%     hold on
%     if numel(inc) == 1
%     ISMMI_2ts_SEM = zeros(1, length(Fitness_level));
%     else  
%     ISMMI_2ts_SEM = std(I_SMMI_2ts(inc,:))/sqrt(length(inc));
%     end
%     xlim([0, max(rangeA)])
%     shadedErrorBar(rangeA(indError), mean(I_SMMI_2ts(inc,indError),1),ISMMI_2ts_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
%     ylabel('I_{SMMI}')
%     %xlabel('#Generations');
%     %ylim([0,1])

subplot(6,4,17+condNum)
    hold on
    if numel(inc) == 1
    IPred_SEM = zeros(1, length(Fitness_level));
    else  
    IPred_SEM = std(I_Pred(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(range)])
    shadedErrorBar(range(indError), mean(I_Pred(inc,indError),1),IPred_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    %plot(rangeA(indError), mean(H_States(inc,indError),1),'color',gcolor)
    ylabel('I_{Pred}')
    %xlabel('#Generations');
    %ylim([0,1])
    
% subplot(6,4,17+condNum)
%     hold on
%     if numel(inc) == 1
%     IPred_2ts_SEM = zeros(1, length(Fitness_level));
%     else  
%     IPred_2ts_SEM = std(I_Pred_2ts(inc,:))/sqrt(length(inc));
%     end
%     xlim([0, max(range)])
%     shadedErrorBar(range(indError), mean(I_Pred_2ts(inc,indError),1),IPred_2ts_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
%     ylabel('I_{Pred}')
%     xlabel('#Generations');

subplot(6,4,21+condNum)
    hold on
    if numel(inc) == 1
    H_States_SEM = zeros(1, length(Fitness_level));
    else  
    H_States_SEM = std(H_States(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(range)])
    shadedErrorBar(range(indError), mean(H_States(inc,indError),1),H_States_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    ylabel('H_{State}')
    xlabel('#Generations');
    %ylim([0,1])
    
%% IIT 2.0

load(strcat(condition,'_phi2008Noise'));
Phi2008 = results.Av_phi_val;
Phi2008MC = results.Av_numElMC;	

figure(106)
subplot(3,4,1+condNum)
  hold on
  if numel(inc) == 1
    FitSEM = zeros(1, length(Fitness_level));
  else  
    FitSEM = std(Fitness_level(inc,:),1)/sqrt(length(inc));
  end  
    shadedErrorBar(rangeA, mean(Fitness_level(inc,indError),1)./(MaxFitness/100), FitSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0);
    xlim([0, max(rangeA)])
    ylim([40,100])
    ylabel('Fitness (%)');
    
subplot(3,4,5+condNum)
    hold on
    if numel(inc) == 1
    Phi2008_SEM = zeros(1, length(Fitness_level));
    else  
    Phi2008_SEM = std(Phi2008(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(rangeA)])
    shadedErrorBar(rangeA(indError), mean(Phi2008(inc,indError),1),Phi2008_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    ylabel('\varphi IIT 2.0')
    %xlabel('#Generations');
    %ylim([0,1])
    
subplot(3,4,9+condNum)
    hold on
    if numel(inc) == 1
    Phi2008MC_SEM = zeros(1, length(Fitness_level));
    else  
    Phi2008MC_SEM = std(Phi2008MC(inc,:))/sqrt(length(inc));
    end
    xlim([0, max(rangeA)])
    shadedErrorBar(rangeA(indError), mean(Phi2008MC(inc,indError),1),Phi2008MC_SEM(indError),{'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
    ylabel('#MC elements IIT 2.0')
    %xlabel('#Generations');
    %ylim([0,1])
    
%% OLD    
%     
% %% Capture Mismatch
% load(strcat(condition, '_Capture_Elem'));
% Cap = results;
% load(strcat(condition, '_dataCB'));
% 
% figure(11)
%     hold on
%     CapSEM = nanstd(Cap.Capture_Change,1)/sqrt(size(Cap.Capture_Change,1));
%     shadedErrorBar(range(indError), nanmean(Cap.Capture_Change(:,indError),1),CapSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
%     %plot(range, Cap.Capture_Change(i,:), '-r')
%     xlim([0, max(range)])
%     ylabel('Matching');
%     xlabel('#Generations');
% 
% %  figure
% %     subplot(1,2,1)
% %     CapSEM = nanstd(Cap.Capture_Change,1)/sqrt(size(Cap.Capture_Change,1));
% %     plot(range, Cap.Capture_Change')
% %     %plot(range, Cap.Capture_Change(i,:), '-r')
% %     xlim([1, max(range)])
% %     ylabel('Matching');
% %     xlabel('#Generations');
% %     subplot(1,2,2)
% %     BPSEM = std(BigPhiMip(inc,:))/sqrt(length(inc));
% %     xlim([1, max(range)])
% %     hold on
% %     plot(range(indError), BigPhiMip(inc,(indError)))
% % ylabel('$_{\rm max}\Phi^{\rm MIP}$','Interpreter','latex','FontSize',12);
% % ylim([0,0.5])
% %%
% %load(strcat(condition,'_UsedNodes'));
% load(strcat(condition,'_ZombiedataAllC'));
% 
% Nodes = zeros(size(UsedNodes));
% for i = 1: size(UsedNodes,1)
%     for j = 1:size(UsedNodes,2)
%         Nodes(i,j) = numel(UsedNodes{i,j});
%     end
% end
% 
% figure(103)
% subplot(2,3,1+condNum)
%     hold on
%     if numel(inc) == 1
%         ConnSEM = zeros(1, length(Fitness_level));
%     else  
%         ConnSEM = std(Num_Conn(inc,:))/sqrt(length(inc));
%     end
%     shadedErrorBar(range(indError)./10000, mean(Num_Conn(inc,(indError)),1),ConnSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
%     xlim([0, max(range)./10000])
%     ylim([0,15])
% ylabel('#connections');
% 
% subplot(2,3,4+condNum)
%     hold on
%     if numel(inc) == 1
%         NodeSEM = zeros(1, length(Fitness_level));
%     else  
%         NodeSEM = std(Nodes(inc,:))/sqrt(length(inc));
%     end
%     shadedErrorBar(range(indError)./10000, mean(Nodes(inc,(indError)),1),NodeSEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0)
%     xlim([0, max(range)./10000])
%     ylim([0,8])
% ylabel('#elements');


