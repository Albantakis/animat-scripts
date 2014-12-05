clear all
option = 0;
FA =  [0 128]; %[0 96 104 112 120 128];
meanVal = 0;
XX = 1:118;
factor = 5;

% Zombie Data
condition = 'c36a45_36';

load(strcat(condition,'_ZombiedataAllC'));
load(strcat(condition,'_ISenMot'));
SMMI = results;
load(strcat(condition,'_IPred'));
IPred = results;
Match = load(strcat(condition,'_Matching'));
%%
if option == 2
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
else
    indF = [1:size(Fitness_level, 1)];
end
totsteps = 60000-1;
gcolor = [0,0.5,1];

%inc = [1:5 7:50]; 
inc = indF %[1:50];
   
%% Only last X Generations
load(strcat(condition,'_ZombiedataAllC'));
MaxFitness = 128;
MeanFitness = mean(Fitness_level(inc,:),1);
MeanFitness = 100.*MeanFitness./MaxFitness;
Fitness_level = 100.*Fitness_level./MaxFitness;

figure
subplot(4,2,1)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(MeanNumConcepts(inc,XX),1,[])]';
        [CorrCon, pCorrCon] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
        %'+','color', gcolor,'DisplayName',condition)
    else
        plot(MeanFitness(XX), mean(MeanNumConcepts(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('#concepts'); 

subplot(4,2,2)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(BigPhi(inc,XX),1,[])]';
        [CorrSumPhi pCorrSumPhi] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
                Counter = 0;

        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
        %'+','color', gcolor,'DisplayName',condition)
    else
        plot(MeanFitness(XX), mean(BigPhi(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12);
    
load(strcat(condition,'_dataCB'));
MaxFitness = 128;
MeanFitness = mean(Fitness_level(inc,:),1);
MeanFitness = 100.*MeanFitness./MaxFitness;
Fitness_level = 100.*Fitness_level./MaxFitness;
  
subplot(4,2,3)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(MeanNumConcepts(inc,XX),1,[])]';
        [CorrConPhi, pCorrConPhi] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
        plot(MeanFitness(XX), mean(MeanNumConcepts(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end    
    xlim([40, 100])
    ylabel('#concepts'); 
 
subplot(4,2,4)
hold on
     if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(BigPhiMip(inc,XX),1,[])]';
        [CorrPhi, pCorrPhi] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
    plot(MeanFitness(XX), mean(BigPhiMip(inc,:)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);

subplot(4,2,5)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(MeanSizeComplex(inc,XX),1,[])]';
        [CorrElem, pCorrElem] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
    plot(MeanFitness(XX), mean(MeanSizeComplex(inc,:)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('#elements');
    
subplot(4,2,6)
    hold on
    Match.Capture_Change(isnan(Match.PhiMatching)) = 0; % NaN are those with Elements\Sen = [] and thus should have Matching = 0;
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(Match.PhiMatching(inc,XX),1,[])]';
        [CorrMatch, pCorrMatch] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
    plot(MeanFitness(XX), mean(Match.PhiMatching(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('Matching');
    
subplot(4,2,7)
    hold on
         if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(SMMI.ISenMot(inc,XX),1,[])]';
        [CorrSMMI, pCorrSMMI] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
    plot(MeanFitness(XX), mean(SMMI.ISenMot(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('SMMI');
    
subplot(4,2,8)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]);  reshape(IPred.IPred(inc,XX),1,[])]';
        [CorrIPred, pCorrIPred] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
        plot(MeanFitness(XX), mean(IPred.IPred(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end    
    xlim([40, 100])
    ylabel('$I_{pred}$');
    
    
Corr = [CorrCon CorrSumPhi CorrConPhi CorrPhi CorrElem CorrMatch CorrSMMI CorrIPred;...
    pCorrCon pCorrSumPhi pCorrConPhi pCorrPhi pCorrElem pCorrMatch pCorrSMMI pCorrIPred]
  
%% Number of connections and Elements
load(strcat(condition,'_ZombiedataAllC'));
MaxFitness = 128;
MeanFitness = mean(Fitness_level(inc,:),1);
MeanFitness = 100.*MeanFitness./MaxFitness;
Fitness_level = 100.*Fitness_level./MaxFitness;

%load(strcat(condition,'_UsedNodes'));

Nodes = zeros(size(UsedNodes));
for i = 1: size(UsedNodes,1)
    for j = 1:size(UsedNodes,2)
        Nodes(i,j) = numel(UsedNodes{i,j});
    end
end

figure(8)
subplot(2,2,1)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(Num_Conn(inc,XX),1,[])]';
        [CorrConnect, pCorrConnect] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
        %'+','color', gcolor,'DisplayName',condition)
    else
        plot(MeanFitness(XX), mean(Num_Conn(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('#Connections'); 
    
subplot(2,2,2)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(Nodes(inc,XX),1,[])]';
        [CorrNodes, pCorrNodes] = corr(Data(:,1), Data(:,2),'type','Spearman');;
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
        %'+','color', gcolor,'DisplayName',condition)
    else
        plot(MeanFitness(XX), mean(Nodes(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('#Elements');     

    Corr2 = [CorrConnect CorrNodes; pCorrConnect pCorrNodes]; 
    

load(strcat(condition,'_UsedNodes'));   
Nodes = zeros(size(UsedNodes));
for i = 1: size(UsedNodes,1)
    for j = 1:size(UsedNodes,2)
        Nodes(i,j) = numel(UsedNodes{i,j});
    end
end

load(strcat(condition,'_data'));
MaxFitness = 128;
MeanFitness = mean(Fitness_level(inc,:),1);
MeanFitness = 100.*MeanFitness./MaxFitness;
Fitness_level = 100.*Fitness_level./MaxFitness;
    
subplot(2,2,3)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(Num_Conn(inc,XX),1,[])]';
        [CorrConnect, pCorrConnect] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
        %'+','color', gcolor,'DisplayName',condition)
    else
        plot(MeanFitness(XX), mean(Num_Conn(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('#Connections'); 

load(strcat(condition,'_UsedNodes'));

subplot(2,2,4)
    hold on
    if meanVal == 0
        Data = [reshape(Fitness_level(inc,XX),1,[]); reshape(Nodes(inc,XX),1,[])]';
        [CorrNodes, pCorrNodes] = corr(Data(:,1), Data(:,2),'type','Spearman');
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
        %'+','color', gcolor,'DisplayName',condition)
    else
        plot(MeanFitness(XX), mean(Nodes(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim([40, 100])
    ylabel('#Elements');   
    
Corr3 = [CorrConnect CorrNodes; pCorrConnect pCorrNodes];
CorrX = [Corr2 Corr3 Corr]

%% Matching vs Phi
load(strcat(condition,'_data'));
MaxFitness = 128;
MeanFitness = mean(Fitness_level(inc,:),1);
MeanFitness = 100.*MeanFitness./MaxFitness;
Fitness_level = 100.*Fitness_level./MaxFitness;

figure(4)
subplot(2,1,1)
hold on
    Match.Capture_Change(isnan(Match.Capture_Change)) = 0; % NaN are those with Elements\Sen = [] and thus should have Matching = 0;
    if meanVal == 0
        Data = [reshape(BigPhiMip(inc,XX),1,[]); reshape(Match.Capture_Change(inc,XX),1,[])]';
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
    plot(mean(Match.Capture_Change(inc,XX)),mean(BigPhiMip(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    %xlim([40, 100])
    ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);

subplot(2,1,2)
hold on
    if meanVal == 0
        Data = [reshape(BigPhiMip(inc,XX),1,[]); reshape(IPred.IPred(inc,XX),1,[])]';
        [A, B, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
    else
    plot(mean(IPred.IPred(inc,XX)),mean(BigPhiMip(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    %xlim([40, 100])
    ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);
