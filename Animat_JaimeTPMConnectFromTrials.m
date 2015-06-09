clear all
Elem = 0:7;
plotflag = 3;
trialnum = [0:9];
numtrials = numel(trialnum);
totsteps = 10500-1;
step = 512;
range = [0:step:totsteps];
cond = 'task1717-3s'

%Foldername = strcat('Freeze_', TrialType,'_trial', int2str(TrialNum));
%mkdir(Foldername)

for t = 1:numtrials
    for i = length(range)
%       ------------- get rest from Phi calculation --------------------------
%       only go to folder here, because I need to access functions from the
%       main folder then, so I need to go out of the folder after loading
        cd(strcat('Freeze_',cond, '_trial', int2str(trialnum(t))))
        Animat_gen = range(i);
        FilenameA = strcat('Animat_gen', int2str(Animat_gen),'_ZombieAllN.mat');
        if exist(FilenameA,'file') == 2 
            load(FilenameA)
            cd ..           
            tpm = results.tpm;
            connectivity_mat = results.connect_mat;
            used_nodes = results.used_nodes;
        else 
            cd ..
        end      
        File = strcat(cond, '_trial', int2str(t-1), '_', int2str(range(i)));
        save(File,'tpm', 'connectivity_mat', 'used_nodes');
    end
end 