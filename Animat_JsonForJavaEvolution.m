%Make struct with #Generation, Connectivity, Fitness, #Concepts and #Phi 
clear all

EvoRun = [32]
cond = 'c36a45_36';

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

% Load files that have the variables you want to display in a trial (=EvoRun) x generation format
load(strcat(cond, '_dataCB'));
BigPhiMip_Complex = BigPhiMip;
load(strcat(cond, '_ZombiedataAllC'));

%Path to files that have the connectivity matrix and used_nodes.
Path = strcat('~/AnimatIIT_FilesPaper1/Freeze_', cond, '_trial', int2str(EvoRun),'/Animat_gen');
 
for g = 1:length(range)
    generations(g).generation = range(g); 
    generations(g).(dataProperties{1}) = Fitness_level(EvoRun+1,g)./128;
    generations(g).(dataProperties{2}) = BigPhiMip_Complex(EvoRun+1,g);
    generations(g).(dataProperties{3}) = MeanNumConcepts(EvoRun+1,g);
    connectMat = zeros(8);
    FilenameA = strcat(Path, int2str(range(g)), '_ZombieAllN.mat');
    generations(g).usedNodes = [];
    if exist(FilenameA,'file') == 2 
        load(FilenameA);
        connectMat(results.used_nodes+1, results.used_nodes+1) = results.connect_mat;
        generations(g).usedNodes = results.used_nodes;
    end
    generations(g).connectivityMatrix = connectMat;
end

data.nodeTypes = nodeTypes;
data.dataLabels = dataLabels;
data.dataAxes = dataAxes;
data.dataProperties = dataProperties;
data.generations = generations;

save(strcat('Animat', int2str(EvoRun), '_Cond', cond), 'data') 
savejson('',data, strcat('Animat', int2str(EvoRun), '_Cond', cond, '.json'))
