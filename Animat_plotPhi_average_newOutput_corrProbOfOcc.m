clear all
Elem = 0:7;
plotflag = 3;
%trialnum = [3 4 15 24 27 35 37 45 49 50 51 53 56 59 67 72 76 81 82 83 84 89 90 95 100 103 104 118 130 131 132 136 145 147 157 163 170 175 182 186 187 191 193 195 197];
%trialnum = [0 3 4 5 6 10 18 27 28 29 31 32 35 38 47 48 49 52 57 58 61 62 64 65 66 68 70 72 74 76 79 83 85 91 92 96 100 110 111 118 125 129 132 133 134 135 141 143 150 152 153 157 158 162 164 167 171 174 180 185 186 192 199];
%trialnum = [0 2 4 7 8 10 12 19 21 23 33 40 43 44 45 48 49 50 52 54 61 70 83 86 87 88 93 109 111 112 115 117 123 126 132 138 147 148 152 159 162 167 170 172 182 191 193 197];
%trialnum = [24 36 47 78 82 92 102 124 137 141 156 163 168 206 215 232 260 264 292 331 346 349 363 390]; %128 and 126
%trialnum = [24 36 47 92 102 137 141 156 163 215 232 264 292 331 346 349 363 390];

trialnum = 0:300;
numtrials = numel(trialnum);
totsteps = 60000-1;
step = 512;
range = 0:step:totsteps;
cond = 'c3a1_change_c14a23';
weighted = 1;

%DPath = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
DPath = '~/dev/animats/results/work_';
path = strcat(DPath, cond, '/trial');
% path3 = strcat(DPath2, cond,'_2', '/trial');
%path2 = strcat(DPath2, cond, '/trial');
% cond2 = strcat(cond, '_2');
% cond3 = cond2;
path2 = path;
cond2 = cond;
%path2 = strcat(DPath2, cond, '_100','/trial');

MaxFitness = 128;

ElemUsed_temp = zeros(1, numel(Elem));
FitPhiCorr = [];

count = 0;
evaluatedTrials = [];
for t = 1:numtrials
    Foldername = strcat('Freeze_',cond, '_trial', int2str(trialnum(t)));
    if exist(Foldername,'dir') == 7 
        count = count+1;
        evaluatedTrials = [evaluatedTrials t-1];
        for i = 1:length(range)
            %------------- get Fitness from Animat files---------------------------
            if trialnum(t) > 99
                path = path2;
            end
            %docname = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_EdgeList.txt');
            docname2 = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_KOdata.txt');

            Fitness = load(docname2);
            Fitness_level(count,i) = Fitness(1);  
    
            %------------- get rest from Phi calculation --------------------------
            % only go to folder here, because I need to access functions from the
            % main folder then, so I need to go out of the folder after loading
            cd(Foldername)
            Animat_gen = range(i);
            FilenameA = strcat('Animat_gen', int2str(Animat_gen),'_ShortResults_CorrBound_noSL.mat');
            %FilenameB = strcat('Animat_gen', int2str(Animat_gen),'_Results.mat');
            if exist(FilenameA,'file') == 2 
                load(FilenameA)
                cd ..
                LSindex = results.LifeState_index;
                % initialize 
                subsys_numel = zeros(numel(LSindex),1);
                Phi = zeros(numel(LSindex),1);
                Phi_MIP = zeros(numel(LSindex),1);
                Num_Concepts = zeros(numel(LSindex),1);
                MConceptOrder = zeros(numel(LSindex),1);
                UsedNodes{count,i} = results.used_nodes;
                for j = 1:numel(LSindex)
                    %subsys = results.Complex{LSindex(j)};
                    subsys = results.Complex{j};
                    sys_results = results.state(j);
                    %Larissa 11/21/2014 Here there is still something wrong
                    %with some option in iit_run that leaves the variables []
                    %but only happens rarely...
                    Phi(j) = sys_results.Phi;
                    % Larissa: This is because it seems that Complex by default is 1
                    if Phi(j) > 0
                        subsys_numel(j) = numel(subsys);
                    else 
                        subsys_numel(j) = 0;
                    end    
                    subsys_ind = subsystem2index(subsys);
                    %sys_results = results.state(LSindex(j));

                    Phi_MIP(j) = sys_results.Phi_MIP;
                    IrrConcepts = [sys_results.concept(:).is_irreducible];
                    Num_Concepts(j) = nnz(IrrConcepts);
                    % higher level concepts
                    numeratorindex = find(IrrConcepts);
                    if Num_Concepts(j) > 0
                        ConceptOrder = zeros(size(numeratorindex));
                        for k = 1:length(numeratorindex)
                             ConceptOrder(k) = numel(index2concept(numeratorindex(k), subsys)) > 1;
                        end
                        MConceptOrder(j) = sum(ConceptOrder);
                    end
                end

                if weighted == 1
                    weights = results.p_LifeTransitions;
                else
                    weights = ones(numel(LSindex),1)./numel(LSindex);
                end
                
                BigPhi(count,i) = sum(weights .* Phi); %weighted ave/expected value
                BigPhiMip(count,i) = sum(weights .* Phi_MIP);  
                
                diffStates(count,i) = length(unique(Phi_MIP));
                NumLifeStates(count,i) = size(results.p_LifeTransitions, 1);
                Num_Conn(count,i) = results.numConn;
                MeanNumConcepts(count,i) = sum(weights .* Num_Concepts);
                MeanSizeComplex(count,i) = sum(weights .* subsys_numel);
                MeanHOConcepts(count,i) = sum(weights .* MConceptOrder);
                FitPhiCorr = [FitPhiCorr; [Fitness(1), BigPhiMip(count,i), range(i)]];
            else 
                cd ..
            end      
        end
    end
end    
%cd ..
%%
Fitness_Phi_Corr = corrcoef(mean(Fitness_level),mean(BigPhiMip))
MeanFitness = mean(Fitness_level,1);
MeanFitness = 100.*MeanFitness./MaxFitness;
if weighted == 1
    File = strcat(cond, '_dataCB_weighted');
else
    File = strcat(cond, '_dataCB_ave');
end
save(File, 'weights', 'evaluatedTrials','BigPhi','BigPhiMip','UsedNodes','Fitness_Phi_Corr', 'Fitness_level', 'MeanFitness', 'MeanNumConcepts', 'MeanHOConcepts', 'MeanSizeComplex', 'NumLifeStates', 'Num_Conn', 'cond', 'range');

if plotflag == 3
    figure
    subplot(3,1,1)
    hold on
    plot(range, MeanFitness)
    xlim([1, totsteps])
    
    subplot(3,1,2)
    hold on
    plot(range, mean(BigPhi,1), '-b')
    plot(range, mean(BigPhiMip,1), '-r')
    xlim([1, max(range)])

    subplot(3,1,3)
    hold on
    %plot(range, NumConcept, '.-b')
    plot(range, mean(NumLifeStates,1), '-b')
    plot(range, mean(Num_Conn,1), '-m')
    plot(range, mean(diffStates,1), '-r')
    plot(range, mean(MeanNumConcepts,1), '-k')
    plot(range, mean(MeanSizeComplex,1), '-g')
    plot(range, mean(MeanHOConcepts,1), '-r')
    xlim([1, max(range)])
    
else
    figure(1)
    subplot(2,1,1)
    hold on
    plot(range, mean(BigPhi,1), '-b')
    plot(range, mean(BigPhiMip,1), '-r')
    xlim([1, max(range)])

    subplot(2,1,2)
    hold on
    %plot(range, NumConcept, '.-b')
    plot(range, mean(diffStates,1), '-r')
    xlim([1, max(range)])
end
