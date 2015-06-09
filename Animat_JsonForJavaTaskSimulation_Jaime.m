function Animat_JsonForJavaTaskSimulation_Jaime
%Make struct with Life TPM decomposed into trials
clear all

blocksUsed = [3 4 6 5];
environmentHeight = 36;

EvoRun = [3 4 15 24 27 35 37 45 49 50];
trials = [1:length(EvoRun)];

cond = 'c1a3_change_c36a45';

numNodes = 8;
generation = 59904;

% Set animat structure.
nodeTypes.sensors = [0, 1];
nodeTypes.hidden = [2, 3, 4, 5];
nodeTypes.motors = [6, 7];

% -------------------------------------------------------------------------
for t = trials
    %Path to EdgeList Files
    Path = strcat('/Users/larissa_alb/dev/animats/results/work_', cond, '/trial', int2str(EvoRun(t)), '_');

    data.generation = generation; 

    %Connectivity matrix from function below.
    data.connectivityMatrix = zeros(numNodes);
    data.usedNodes = [];
    data.nodeTypes = nodeTypes;
    
    animatResults = Animat_Connectivity(length(nodeTypes.sensors), length(nodeTypes.motors), numNodes, generation, Path); 
    if ~isempty(animatResults)
        %Make connectivity matrix 8x8
        data.connectivityMatrix(animatResults.usedNodes+1, animatResults.usedNodes+1) = animatResults.connect_mat;
        data.usedNodes = animatResults.usedNodes;
    end
    
    %Life Time TPM decomposed into 128 individual trials
    docname = strcat(Path,  int2str(generation), '_LifetimeLogicTable.txt');
    LifeTimeAll = importdata(docname,',', 1);
    %Take sensors from before transition and hidden and motors from after
    %transitions
    Sensors = LifeTimeAll.data(:,nodeTypes.sensors+1);
    HiddenMotors = [LifeTimeAll.data(:, (10+length(nodeTypes.sensors)):end)];
    LifeTime = [Sensors, HiddenMotors];
    for g = 1:floor(length(LifeTime)/environmentHeight)
        data.Trial(g).trialNum = g;
        data.Trial(g).lifeTable = LifeTime((1+(g-1)*environmentHeight):(g*environmentHeight),:);
        %clear motors at first time step in new trial

        data.blockSize(g) = blocksUsed(ceil(g/(length(LifeTime)/(environmentHeight*length(blocksUsed)))));
    end

    save(strcat('AnimatLife', int2str(EvoRun(t)), '_', cond, '_', int2str(generation)), 'data') 
    savejson('',data, strcat('AnimatLife', int2str(EvoRun(t)), '_', cond, '_', int2str(generation), '.json'))
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
