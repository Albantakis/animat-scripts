function Animat_PredInfo
tic
clear all

Elem = 0:7;
numSen = 2;
Mot = [7 8];
%plotflag = 3;
trialnum = [0:49];
%trialnum = setdiff([0:9],[0 2 6 8 9]);
%trialnum = [0 5 8 10 12 13 14 15 19];
numtrials = numel(trialnum);
totsteps = 60000-1;
step = 512;
range = [0:step:totsteps];
%range = [0:step:totsteps/2 29984 30208:step:totsteps 59984];
%startpoints = [696, 189, 349, 1692, 2356, 4372];
%cond = 'c24a35_36';
%cond = 'c1a3_36';
%cond = 'c1a3_n001_90000_36';
%cond = 'c1b11a2b5_36';
%cond = 'c24a35_36';
%cond = 'c35a271_36';
cond = 'c1a3_n001_rep20_36';
%cond = 'c14a23_36';
%DPath2 = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
path = strcat(DPath, cond, '/trial');
path2 = path;
cond2 = cond;
path3 = path;
cond3 = cond;
% path2 = strcat(DPath, cond, '_2b', '/trial');
% path3 = strcat(DPath, cond, '_2', '/trial');
% cond2 = strcat(cond, '_2b');
% cond3 = strcat(cond, '_2');

%MaxFitness = 128;
Fitness_level = zeros(numtrials, length(range));
IPred = zeros(numtrials, length(range));
HStates = zeros(numtrials, length(range));

for t = 1:numtrials
    t
    if trialnum(t) > 9
        APath = strcat(path2, int2str(trialnum(t)), '_');
         if trialnum(t) > 19
            APath = strcat(path3, int2str(trialnum(t)), '_');
            cond = cond3;
         end
    else 
        APath = strcat(path, int2str(trialnum(t)), '_');
    end  
    for i = 1:length(range)
        %------------- get Fitness from Animat files---------------------------
         docname2 = strcat(APath, int2str(range(i)), '_KOdata.txt');
         Fitness = load(docname2);
         Fitness_level(t,i) = Fitness(1);  
        %Calculate probabilities for Sensory Motor predictive Information
        [p_States, p_Joint_t0t1, ~] = Animat_LifeTimeDist(range(i), APath, Mot);
        p_Prod_t0t1 = p_States*p_States';   %This works if the ps are column vectors and states are ordered 000 100 010 110 etc
        p_Prod_t0t1 = p_Prod_t0t1(:);
        IPred(t,i) = KLD(p_Joint_t0t1, p_Prod_t0t1); %Mutual information of Sensor and Motors
        
        p_States(p_States == 0) = 1;
        HStates(t,i) = DistEntropy(p_States);
    end
end

figure
subplot(2,1,1)
 hold on
    plot(mean(Fitness_level), mean(IPred), '+r')
    xlim([1, 128])
subplot(2,1,2)
    hold on
    plot(range, mean(IPred), '-r')
    xlim([1, max(range)])

%------------------------ Save results ------------------------------------ 
results.Fitness = Fitness_level;
results.range = range;
results.IPred = IPred;
results.HStates = HStates;
save(strcat(cond, '_IPred'), 'results');
   
toc
end

function [p_States, p_JointStates, SizeLife] = Animat_LifeTimeDist(gen, APath,Mot)  
    fullstatefile = strcat(APath, int2str(gen), '_', 'Lifetime' ,'LogicTable.txt');
    Full_tpm = importdata(fullstatefile,',', 1);
    % Input Probability
    %--------------------------------------------
    States = [Full_tpm.data(2:end,1:6) Full_tpm.data(1:end-1,Mot+9)];
    States = [Full_tpm.data(1,[1:6 Mot]); States];
    SizeLife = size(States, 1);
    [DistStates, ~, StMap] = unique(States, 'rows', 'First');
    for i = 1:size(DistStates, 1)
        %Order input States according to place in tpm
        StInd(i) = state2index(DistStates(i,:), 2.*ones(size(DistStates(1,:))));
    end
    StMap = StInd(StMap);
    States_distr = hist(StMap, 1:2^numel([1:6 Mot]))';
    p_States = States_distr./sum(States_distr,1);
       
    %Joint Probability
    %--------------------------------------------
    Joint_tpm = [States(1:end-1,:) States(2:end,:)];
    %Joint_tpm = [States(1:end-2,:) States(3:end,:)]; 2 time steps
    ind = 1:size(Joint_tpm,1);
    ind = setdiff(ind, 36:36:size(Joint_tpm,1));
    Joint_tpm = Joint_tpm(ind,:);
    [JointStates, ~, JointMap] = unique(Joint_tpm, 'rows', 'First');
    for i = 1:size(JointStates, 1)
        %Order output States according to place in tpm
        JointInd(i) = state2index(JointStates(i,:), 2.*ones(size(JointStates(1,:))));
    end
    JointMap = JointInd(JointMap);
      
    %Translate number of states into stateindex and make histogram
    JointStates_distr = hist(JointMap, 1:2^size(Joint_tpm,2))';
    p_JointStates = JointStates_distr./sum(JointStates_distr,1);
end    

function entropy = DistEntropy(Dist)
    entropy = -sum(Dist.*log2(Dist));
end

function index = state2index(state_vec, num_states_vec)
index = state_vec(1) + 1;
for i = 2:length(num_states_vec)
    index = index + state_vec(i)*prod(num_states_vec(1:i-1));  
end
end