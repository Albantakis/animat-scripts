function Animat2_plotDataShadedError_long
clear all
option = 4;
condNum = 0;
% condition = 'c14a23_36';
% condition = 'c23a14_36';
% condition = 'c36a45_36';
condition = 'c1a2_change_c14a23';
%% Data
cd('old_mixed_weighted_ave_mean')
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
    gcolor = [0.75 0.75 0.5];
    %gcolor = [0 0.5 1];
    %gcolor = [1/i,1-1/i,0];
    %gcolor = [1 .7 0]; %oneSen all
    %gcolor = [0 0.7 0]; %oneMot all
    %gcolor = [1,0,1]; %noise all
end

totsteps = 60000-1;

length(indF)
inc = indF %[1:50];
MaxFitness = 128;
MeanFitness = mean(Fitness_level(inc,:),1);
MeanFitness = 100.*MeanFitness./MaxFitness;
Fitness_level = 100.*Fitness_level./MaxFitness;

%% Zombie Data
load(strcat(condition,'_ZombiedataAllC'));
rangeB = eval('range');
indError = 1:length(rangeB);
rangeA = rangeB(indError)./10000;
figure(12)

subplot(8,4,1+condNum)
    plotVariable(inc, indError, rangeA, Fitness_level./(MaxFitness/100), gcolor, condition, [40,100])
    ylabel('Fitness (%)');
subplot(8,4,5+condNum)
    plotVariable(inc, indError, rangeA, MeanNumConcepts, gcolor, condition, [0,4])
    ylabel('#concepts');
subplot(8,4,9+condNum)
    plotVariable(inc, indError, rangeA, MeanHOConcepts, gcolor, condition, [0,1.5])
    ylabel('#HO concepts');
subplot(8,4,13+condNum)
    plotVariable(inc, indError, rangeB, BigPhi, gcolor, condition, [0,1])
    ylabel('$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12);

%% Complex Data
load(strcat(condition, '_dataCB'));
indError = 1:length(rangeB);

subplot(8,4,17+condNum)
    plotVariable(inc, indError, rangeA, MeanSizeComplex, gcolor, condition, [0,1.5])
    ylabel('#elements');
subplot(8,4,21+condNum)
    plotVariable(inc, indError, rangeA, MeanNumConcepts, gcolor, condition, [0,2])
    ylabel('#concepts');
subplot(8,4,25+condNum)
    plotVariable(inc, indError, rangeA, MeanHOConcepts, gcolor, condition, [0,0.3])
    ylabel('#HO concepts');
subplot(8,4,29+condNum)
    plotVariable(inc, indError, rangeB, BigPhiMip, gcolor, condition, [0,0.4])
    ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);
    xlabel('#Generations');

cd ..
%% I_SMMI and I_Pred
cd('I_sen_pred') 
load(strcat(condition,'_ISenMot'));
I_SMMI = results.ISenMot;
HMot = results.HMot;
HSen = results.HSen;

load(strcat(condition,'_IPred'));
I_Pred = results.IPred;
H_States = results.HStates;
cd ..

figure(102)
subplot(6,4,1+condNum)
    plotVariable(inc, indError, rangeA, Fitness_level./(MaxFitness/100), gcolor, condition, [40,100])
    ylabel('Fitness (%)');
subplot(6,4,5+condNum)
    plotVariable(inc, indError, rangeA, I_SMMI, gcolor, condition, [0,1])
    ylabel('I_{SMMI}')
subplot(6,4,9+condNum)
    plotVariable(inc, indError, rangeA, HSen, gcolor, condition, [0,2])
    ylabel('H_{Sen}')   
subplot(6,4,13+condNum)
    plotVariable(inc, indError, rangeA, HMot, gcolor, condition, [0,2])
    ylabel('H_{Mot}')
subplot(6,4,17+condNum)
    plotVariable(inc, indError, rangeA, I_Pred, gcolor, condition, [0,3])
    ylabel('I_{Pred}')
subplot(6,4,21+condNum)
    plotVariable(inc, indError, rangeB, H_States, gcolor, condition, [0,4])
    ylabel('H_{State}')
    xlabel('#Generations');
end
    
function plotVariable(inc, indError, rangeX, variable, gcolor, condition, yLimits)
  hold on
  if numel(inc) == 1
    SEM = zeros(1, length(variable));
  else  
    SEM = std(variable(inc,:),1)/sqrt(length(inc));
  end  
    shadedErrorBar(rangeX(indError), mean(variable(inc,indError),1), SEM(indError), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0);
    xlim([0, max(rangeX)])
    ylim(yLimits);
end   