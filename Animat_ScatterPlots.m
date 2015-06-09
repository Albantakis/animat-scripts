function Animat_ScatterPlots
clear all
option = 0;
meanVal = 0;
MaxFitness = 128;

global XX
global factor
global condition

XX = 118;
factor = 5;
condition = 'c1a3_36';
condition = 'c1a2_36';
condition = 'c14a23_36';
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

if strcmp(condition, 'c1a3_36')
        gcolor = [0,0,0];
elseif strcmp(condition, 'c1a2_36')
        gcolor = [0,0,0.5];
elseif strcmp(condition, 'c14a23_36')
        gcolor = [0,0,1];
elseif strcmp(condition, 'c36a45_36')
        %gcolor = [0,0.5,1];
        gcolor = [0.5,0,0];
end

%inc = [1:5 7:50];
inc = indF %[1:50];
   
%% Only last X Generations
load(strcat(condition,'_ZombiedataAllC'));
Fitness_level = 100.*Fitness_level./MaxFitness;

figure(100)
subplot(4,2,1)
    hold on
    [Corr(1), pCorr(1)] = ScatterPlot(meanVal,Fitness_level,MeanNumConcepts,inc,[40, 100], gcolor);
    ylabel('#concepts'); 

subplot(4,2,2)
    hold on
    hold on
    [Corr(2), pCorr(2)] =ScatterPlot(meanVal,Fitness_level,BigPhi,inc,[40, 100],gcolor);    
    ylabel('$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12);
    
load(strcat(condition,'_dataCB'));
Fitness_level = 100.*Fitness_level./MaxFitness;
  
subplot(4,2,3)
    hold on
    [Corr(3), pCorr(3)] = ScatterPlot(meanVal,Fitness_level,MeanNumConcepts,inc,[40, 100], gcolor);
    ylabel('#concepts'); 
 
subplot(4,2,4)
    hold on
    [Corr(4), pCorr(4)] = ScatterPlot(meanVal,Fitness_level,BigPhiMip,inc,[40, 100], gcolor);
    ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);

subplot(4,2,5)
    hold on
    [Corr(5), pCorr(5)] = ScatterPlot(meanVal,Fitness_level,MeanSizeComplex,inc,[40, 100], gcolor);
    ylabel('#elements');
        
subplot(4,2,6)
    hold on
    [Corr(6), pCorr(6)] = ScatterPlot(meanVal,Fitness_level,SMMI.ISenMot,inc,[40, 100], gcolor);       
    ylabel('SMMI');
    
subplot(4,2,7)
    hold on
    [Corr(7), pCorr(7)] = ScatterPlot(meanVal,Fitness_level,IPred.IPred,inc,[40, 100], gcolor);       
    ylabel('$I_{pred}$');
    
% subplot(4,2,8)
%     hold on
%     [Corr(8), pCorr(8)] = ScatterPlot(meanVal,Fitness_level,Match.PhiMatching,inc,[40, 100], gcolor);
%     xlim([40, 100])
%     ylabel('Matching');
    
[Corr pCorr]
  
%% Number of connections and Elements
load(strcat(condition,'_ZombiedataAllC'));
Fitness_level = 100.*Fitness_level./MaxFitness;

%load(strcat(condition,'_UsedNodes'));

Nodes = zeros(size(UsedNodes));
for i = 1: size(UsedNodes,1)
    for j = 1:size(UsedNodes,2)
        Nodes(i,j) = numel(UsedNodes{i,j});
    end
end

figure(101)
subplot(2,2,1)
    hold on
    [Corr2(1), pCorr2(1)] = ScatterPlot(meanVal,Fitness_level,Num_Conn,inc,[40, 100], gcolor);       
    ylabel('#Connections'); 
    
subplot(2,2,2)
    hold on
    [Corr2(2), pCorr2(2)] = ScatterPlot(meanVal,Fitness_level,Nodes,inc,[40, 100], gcolor);       
    ylabel('#Elements');     

[Corr2 pCorr2]
    

load(strcat(condition,'_UsedNodes'));   
Nodes = zeros(size(UsedNodes));
for i = 1: size(UsedNodes,1)
    for j = 1:size(UsedNodes,2)
        Nodes(i,j) = numel(UsedNodes{i,j});
    end
end

load(strcat(condition,'_dataCB'));
Fitness_level = 100.*Fitness_level./MaxFitness;
    
subplot(2,2,3)
    hold on
    [Corr3(2), pCorr3(2)] = ScatterPlot(meanVal,Fitness_level,Num_Conn,inc,[40, 100], gcolor);       
    ylabel('#Connections'); 

load(strcat(condition,'_UsedNodes'));

subplot(2,2,4)
    hold on
    [Corr3(2), pCorr3(2)] = ScatterPlot(meanVal,Fitness_level,Nodes,inc,[40, 100], gcolor);       
    ylabel('#Elements');   
    
[Corr3 pCorr3]
%% ISMMI IPred vs Fitness and Phi
load(strcat(condition,'_ZombiedataAllC'));
Fitness_level = 100.*Fitness_level./MaxFitness;

figure(104)
subplot(1,6,3)
    hold on
    [Corr4(1), pCorr4(1)] = ScatterPlot(meanVal,Fitness_level,SMMI.HSen,inc,[40, 100], gcolor);       
    xlabel('Fitness')
    %ylabel('$I_{SMMI}$','Interpreter','latex','FontSize',12);
    ylabel('HSen')

subplot(1,6,4)
    hold on
    [Corr4(2), pCorr4(2)] = ScatterPlot(meanVal,Fitness_level,IPred.HStates,inc,[40, 100], gcolor);       
    xlabel('Fitness')
    %ylabel('$I_{pred}$','Interpreter','latex','FontSize',12);
    ylabel('H')


subplot(1,6,1)
    hold on
    [Corr4(3), pCorr4(3)] = ScatterPlot(meanVal,Fitness_level,MeanNumConcepts,inc,[40, 100], gcolor);
    xlabel('Fitness')
    ylabel('#concepts'); 

load(strcat(condition,'_dataCB'));
Fitness_level = 100.*Fitness_level./MaxFitness;

subplot(1,6,2)
    hold on
    [Corr4(4), pCorr4(4)] = ScatterPlot(meanVal,Fitness_level,BigPhiMip,inc,[40, 100], gcolor);
    xlabel('Fitness')
    ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);

    
subplot(1,6,5)
hold on
%     Match.Capture_Change(isnan(Match.Capture_Change)) = 0; % NaN are those with Elements\Sen = [] and thus should have Matching = 0;
%     [Corr4(1), pCorr4(1)] = ScatterPlot(meanVal,BigPhiMip,Match.Capture_Change,inc,[0, 1.5], gcolor)       
    [Corr4(5), pCorr4(5)] = ScatterPlot(meanVal,BigPhiMip,SMMI.HSen,inc,[0, 1.5], gcolor);      
    xlabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);
    %ylabel('$I_{SMMI}$','Interpreter','latex','FontSize',12);
    ylabel('HSen')
    
subplot(1,6,6)
hold on
    [Corr4(6), pCorr4(6)] = ScatterPlot(meanVal,BigPhiMip,IPred.HStates,inc,[0, 1.5], gcolor);       
    xlabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);
    %ylabel('$I_{pred}$','Interpreter','latex','FontSize',12);
    ylabel('H')
    [Corr4; pCorr4]
end
    
function [Corr, pCorr] = ScatterPlot(meanVal,plotdata1,plotdata2,inc,xRange,gcolor)
    global XX
    global factor
    global condition
    if meanVal == 0
        Data = [reshape(plotdata1(inc,XX),1,[]); reshape(plotdata2(inc,XX),1,[])]';
        [Corr, pCorr] = corr(Data(:,1), Data(:,2),'type','Spearman');
        pCorr = pCorr*40;
        [A, ~, C] = unique(Data,'rows');
        Counter = 0;
        for i = 1:size(A,1)
           Counter(i) = numel(find(C==i))*factor; 
        end
        scatter(A(:,1),A(:,2),Counter,gcolor)
        %'+','color', gcolor,'DisplayName',condition)
    else
        plot(mean(plotdata1(inc,XX)), mean(MeanNumConcepts(inc,XX)), '+', 'color', gcolor,'DisplayName',condition)
    end
    xlim(xRange)
end