function Animat_JsonForJavaEvolution_Jaime
%Make struct with #Generation, Connectivity, Fitness, #Concepts and #Phi 
clear all

EvoRun = [3 4 15 24 27 35 37 45 49 50];
trials = [1:45];

cond = 'c1a3_change_c36a45';

numNodes = 8;
step = 512;
totsteps = 60000;
range = [0:step:totsteps];

% Set animat structure.
nodeTypes.sensors = [0, 1];
nodeTypes.hidden = [2, 3, 4, 5];
nodeTypes.motors = [6, 7];

% Label the variable you want to display.
dataLabels = {'Fitness'; 'Phi'; 'Number of Concepts'};
dataProperties = {'fitness'; 'phi'; 'numConcepts'};

dataAxes.(dataProperties{1}) = [0, 1];
dataAxes.(dataProperties{2}) = [0, 1.25];
dataAxes.(dataProperties{3}) = [0, 8];

% -------------------------------------------------------------------------
% Load files that have the variables you want to display in a trial (=EvoRun) x generation format
load(strcat(cond, '_dataCB'));
BigPhiMip_Complex = BigPhiMip;
load(strcat(cond, '_ZombiedataAllC'));
% -------------------------------------------------------------------------

for t = trials
    %Path to EdgeList Files
    Path = strcat('/Users/larissa_alb/dev/animats/results/work_', cond, '/trial', int2str(EvoRun(t)), '_');

    for g = 1:length(range)
        generations(g).generation = range(g); 

        generations(g).(dataProperties{1}) = Fitness_level(t,g)./128;
        generations(g).(dataProperties{2}) = BigPhiMip_Complex(t,g);
        generations(g).(dataProperties{3}) = MeanNumConcepts(t,g);

        %Connectivity matrix from function below.
        generations(g).connectivityMatrix = zeros(numNodes);
        generations(g).usedNodes = [];

        generationResults = Animat_Connectivity(length(nodeTypes.sensors), length(nodeTypes.motors), numNodes, range(g), Path); 
        if ~isempty(generationResults)
            %Make connectivity matrix 8x8
            generations(g).connectivityMatrix(generationResults.usedNodes+1, generationResults.usedNodes+1) = generationResults.connect_mat;
            generations(g).usedNodes = generationResults.usedNodes;
        end
    end

    data.nodeTypes = nodeTypes;
    data.dataLabels = dataLabels;
    data.dataAxes = dataAxes;
    data.dataProperties = dataProperties;
    data.generations = generations;

    save(strcat('Animat', int2str(EvoRun(t)), '_', cond), 'data') 
    savejson('',data, strcat('Animat', int2str(EvoRun(t)), '_', cond, '.json'))
    end
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
