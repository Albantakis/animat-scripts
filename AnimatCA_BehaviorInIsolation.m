function AnimatCA_BehaviorInIsolation
%Make struct with Life TPM decomposed into trials
clear all

EvoRun = [0:49];
trials = [1:length(EvoRun)];

condition = 'c14a23_36';

numNodes = 8;
generation = 59904;

% Set animat structure.
nodeTypes.sensors = [0, 1];
nodeTypes.motors = [6, 7];
numSen = numel(nodeTypes.sensors);

avTransientLength = zeros(length(trials),1);
maxTransientLength = zeros(length(trials),1);
avPeriodic = zeros(length(trials),1);
maxPeriodic = zeros(length(trials),1);
% -------------------------------------------------------------------------
for t = trials
    %Path to EdgeList Files
    %Path = strcat('/Users/larissa_alb/dev/animats/results/work_', cond, '/trial', int2str(EvoRun(t)), '_');
    Path = strcat('/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_', condition, '/trial', int2str(EvoRun(t)), '_');   
    load(strcat(condition, '_dataCB'))
    
    animatResults = Animat_Connectivity(length(nodeTypes.sensors), length(nodeTypes.motors), numNodes, generation, Path); 
    if ~isempty(animatResults)  
        % Convert TPM from binary to state index ------------------------------
        tpm = animatResults.tpm;
        % set Sensors to 0
        tpm(:,animatResults.usedNodes < numSen) = 0;
        % LOLI ordering, that is why we need to flip the TPM before bin2dec
        indexTPM = bin2dec(int2str(fliplr(tpm)));

        % Simulation to measure convergence to steady state
        % Set animat in all possible states. Keep Sensors the way they were and
        % evolve according to tpm.
        transientLength = zeros(length(indexTPM),1);
        periodicTemp = zeros(length(indexTPM),1);
        for s = 0:length(indexTPM)-1
            % Sensors are low index nodes, thus add 0, 1, 2, or 3, for sensor 
            % state 00, 10, 01, and 11 to indexTPM
            senState = mod(s, 4); 
            adjustedIndexTPM = indexTPM + senState;
            % LOLI ordering as in IIT 
            state = s;
            pastState = NaN;
            stateEvolution = state;
            while pastState ~= state && length(stateEvolution) <= length(indexTPM)
                state = adjustedIndexTPM(state+1);
                ind = find(stateEvolution == state, 1, 'first');
                if ~isempty(ind)
                    pastState = stateEvolution(ind);
                    if ind < length(stateEvolution)
                        periodicTemp(s+1) = length(stateEvolution)+1-ind;
                    end
                end
                stateEvolution = [stateEvolution; state];
            end
            if length(stateEvolution) > length(indexTPM)
                [EvoRun(t) s]
                stateEvolution'
            else
                transientLength(s+1) = length(stateEvolution)-2;
            end
        end
        avTransientLength(t) = mean(transientLength);
        maxTransientLength(t) = max(transientLength);
        avPeriodic(t) = mean(periodicTemp);
        maxPeriodic(t) = max(periodicTemp);
    end
end

figure
subplot(3,1,1)
plotCorrelation(Fitness_level(EvoRun+1,end), avTransientLength)
subplot(3,1,2)
plotCorrelation(BigPhiMip_B(EvoRun+1,end), avTransientLength)

%load(strcat(condition,'_ZombiedataAllC'))
subplot(3,1,3)
plotCorrelation(MeanNumConcepts(EvoRun+1,end), avTransientLength)

mean(avTransientLength)
A = [EvoRun' Fitness_level(EvoRun+1,end) BigPhiMip(EvoRun+1,end) maxTransientLength avTransientLength avPeriodic maxPeriodic];
end


% -------------------------------------------------------------------------
% Function to get connectivity matrix out of EdgeList
function results = Animat_Connectivity(numSen, numMot, numNodes, generation, AnimatPath) 
%TrialType = 'task1-3s-01102015'
%AnimatPath = strcat('/Users/larissa_alb/Dropbox/larissa-jaime/extrinsic-cause-info/', TrialType, '/trial', int2str(TrialNum), '_');
results = [];
J_tempfile = strcat(AnimatPath, int2str(generation), '_EdgeList.txt');   
%/Users/larissa_alb/dev/animats/results/work_c1a3_change_c36a45
if exist(J_tempfile,'file') == 2 
    J_temp = load(J_tempfile);
    J_temp = unique(J_temp, 'rows')+1;
    if ~isempty(J_temp)
        % MOTORS ARE SET TO 0 -> they don't actually have recurrent connections
        J_temp = J_temp(J_temp(:,1) <= numNodes-numMot,:); 
        %Sensors might have incoming connection, but because they actually don't do anything they shouldn't be taken into account
        J_temp = J_temp(J_temp(:,2) > numSen,:);

        J_sparse = sparse(J_temp(:,1), J_temp(:,2),1,numNodes,numNodes);

        J = full(J_sparse)';     

        [tpm, usedNodes] = Animat_ReducedTpmSmall(generation, AnimatPath, numSen, numMot, J, 0);
        tpm(:,usedNodes < numSen+1) = 0.5; %Sensors might be switched on through mechanism but that is overwritten by environment --> doesn't do anything

        J = J(usedNodes, usedNodes);
        if ~isempty(J)
            results.numConn = nnz(J);
            results.connect_mat = J;
            results.usedNodes = usedNodes-1;
            results.tpm = tpm; 
        end
    end
end
end

function plotCorrelation(data1, data2)
    plot(data1, data2, 'x')
    [CorrCoef p] = corr(data1,data2,'type','Spearman');
    [CorrCoef p*40]
end
