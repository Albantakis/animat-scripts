%Make struct with #Generation, Connectivity, Fitness, #Concepts and #Phi 
clear all
%for EvoRun = [32]%[4, 15, 27, 42, 43, 47];
EvoRun = [32]
cond = 'c36a45_36';
numSen = 2;

load(strcat(cond, '_dataCB'));
BigPhiMip_Complex = BigPhiMip;
load(strcat(cond, '_ZombiedataAllC'));
Path = strcat('~/AnimatIIT_FilesPaper1/Freeze_', cond, '_trial', int2str(EvoRun),'/Animat_gen');
Path2 = strcat('/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_', cond, '/trial');
%%
% for g = 1:length(range)
%     Animat(g).generation = range(g);
%     Animat(g).fitness = Fitness_level(EvoRun+1,g);
%     Animat(g).numConcepts = MeanNumConcepts(EvoRun+1,g);
%     Animat(g).phi = BigPhiMip_Complex(EvoRun+1,g);
%     connectMat = zeros(8);
%     FilenameA = strcat(Path, int2str(range(g)), '_ZombieAllN.mat');
%     Animat(g).usedNodes = [];
%     if exist(FilenameA,'file') == 2 
%         load(FilenameA);
%         connectMat(results.used_nodes+1, results.used_nodes+1) = results.connect_mat;
%         Animat(g).usedNodes = results.used_nodes;
%     end
%     Animat(g).connectivityMatrix = connectMat;
% end
% 
% save(strcat('Animat', int2str(EvoRun)), 'Animat') 
% savejson('./',Animat,strcat('Animat', int2str(EvoRun), '.json'))
% clear all
% end 
%% Get LifeTimeTable
g = length(range);
blocksUsed = [3 4 6 5];
AnimatGen.generation = range(g);
connectMat = zeros(8);
FilenameA = strcat(Path, int2str(range(g)), '_ZombieAllN.mat');
AnimatGen.usedNodes = {'empty'};
if exist(FilenameA,'file') == 2 
    load(FilenameA);
    connectMat(results.used_nodes+1, results.used_nodes+1) = results.connect_mat;
    AnimatGen.usedNodes = results.used_nodes;
end
AnimatGen.connectivityMatrix = connectMat;
docname = strcat(Path2, int2str(EvoRun), '_', int2str(range(g)), '_LifetimeLogicTable.txt');
LifeTimeAll = importdata(docname,',', 1);
Sensors = LifeTimeAll.data(:,[1:numSen]);
HiddenMotors = [LifeTimeAll.data(:, (10+numSen):end)];
LifeTime = [Sensors, HiddenMotors];
for t = 1:floor(length(LifeTime)/36)
    AnimatGen.Trial(t).trialNum = t;
    AnimatGen.Trial(t).lifeTable = LifeTime((1+(t-1)*36):(t*36),:);
    %clear motors at first time step in new trial
    
    AnimatGen.blockSize(t) = blocksUsed(ceil(t/(length(LifeTime)/(36*4))));
end
save(strcat('AnimatGen', int2str(EvoRun)), 'AnimatGen') 
savejson('./',AnimatGen,strcat('AnimatBlockTrials', int2str(EvoRun),'_',int2str(range(g)), '.json'))