%Make struct with #Generation, Connectivity, Fitness, #Concepts and #Phi 

EvoRun = 32;
cond = 'c36a45_36';

load(strcat(cond, '_dataCB'));
BigPhiMip_Complex = BigPhiMip;
load(strcat(cond, '_ZombiedataAllC'));
Path = strcat('~/AnimatIIT_FilesPaper1/Freeze_', cond, '_trial', int2str(EvoRun),'/Animat_gen');

for g = 1:length(range)
    Animat(g).generation = range(g);
    Animat(g).fitness = Fitness_level(EvoRun+1,g);
    Animat(g).numConcepts = MeanNumConcepts(EvoRun+1,g);
    Animat(g).bigPhiMip = BigPhiMip_Complex(EvoRun+1,g);
    connectMat = zeros(8);
    FilenameA = strcat(Path, int2str(range(g)), '_ZombieAllN.mat');
    Animat(g).usedNodes = [];
    if exist(FilenameA,'file') == 2 
        load(FilenameA);
        connectMat(results.used_nodes+1, results.used_nodes+1) = results.connect_mat;
        Animat(g).usedNodes = results.used_nodes;
    end
    Animat(g).connectivityMatrix = connectMat;
end

save(strcat('Animat', int2str(EvoRun)), 'Animat') 
savejson('./',Animat,strcat('Animat', int2str(EvoRun)))