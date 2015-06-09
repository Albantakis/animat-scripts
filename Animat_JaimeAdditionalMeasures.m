function Animat_JaimeAdditionalMeasures
clear all
Elem = 0:7;
loadflag = 0;
saveflag = 0;
plotflag = 3;
trialnum = [10:19];
numtrials = numel(trialnum);

cond = 'task7156331-2s'
measure = 'MI'

switch measure
    case 'MI'
        step = 512;
        totsteps = 10500-1;
        DPath = '/Users/larissa_alb/Dropbox/larissa-jaime/mutualinformation/';
    case 'Phi'
        step = 64;
        totsteps = 1050-1;
        DPath = '/Users/larissa_alb/Dropbox/larissa-jaime/bigPhi/'
    case 'extrCauseInfo'
        totsteps = 1050-1;
        step = 64;
        DPath = '/Users/larissa_alb/Dropbox/larissa-jaime/extrinsic-cause-info/';
end
range = [0:step:totsteps];

path = strcat(DPath, cond, '/trial');
path2 = path;
cond2 = cond;

MaxFitness = 128;

Fitness_level = zeros(numtrials,length(range));
HSen = zeros(numtrials,length(range));
MI = zeros(numtrials,length(range));
Phi = zeros(numtrials,length(range));
extrCauseInfo = zeros(numtrials,length(range));

UsedNodes = cell(numtrials,length(range));

for t = 1:numtrials
    for i = 1:length(range)
        %------------- get Fitness from Animat files---------------------------
        docname = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_KOdata.txt');        
        Fitness = load(docname);
        Fitness_level(t,i) = Fitness(1);
        switch measure
            case 'MI'
                HSen(t,i) = Fitness(3);
                MI(t,i) = Fitness(4);
            case 'Phi'
                HSen(t,i) = Fitness(2);
                MI(t,i) = Fitness(3);
                Phi(t,i) = Fitness(4);
            case 'extrCauseInfo'
                HSen(t,i) = Fitness(2);
                MI(t,i) = Fitness(3);
                extrCauseInfo(t,i) = Fitness(4);
        end
    end
end    

MeanFitness = mean(Fitness_level,1);
MeanFitness = 100.*MeanFitness./MaxFitness;

if saveflag == 1
if loadflag == 1
    HSenTemp = HSen; MITemp = MI;
    File = strcat(cond, '_ZombiedataAllC');
    load(File)
    clear HSen
    clear MI
    clear Phi
    clear extrCauseInfor
    HSen = HSenTemp; MI = MITemp;
        switch measure
            case 'MI'
                save(File,'Fitness_level','HSen', 'MI', '-append');

            case 'Phi'
                save(File,'Fitness_level','HSen', 'MI', 'Phi', '-append');
            case 'extrCauseInfo'
                save(File,'Fitness_level','HSen', 'MI', 'extrCauseInfo', '-append');
        end
else
    switch measure
            case 'MI'
                File = strcat(cond, '_MI');
                save(File,'Fitness_level','HSen', 'MI');
            case 'Phi'
                File = strcat(cond, '_Phi');
                save(File,'Fitness_level','HSen', 'MI', 'Phi');
            case 'extrCauseInfo'
                File = strcat(cond, '_extrCauseInfo');
                save(File,'Fitness_level','HSen', 'MI', 'extrCauseInfo');
    end
end
end
%%
if plotflag == 3
    figure
    subplot(3,1,1)
    hold on
    plot(range, MeanFitness)
    xlim([1, totsteps])
    
    subplot(3,1,2)
    plotPlusMean(range, HSen)
    
    subplot(3,1,3)
    switch measure
        case 'MI'
            plotPlusMean(range, MI)
        case 'Phi'
            plotPlusMean(range, Phi)
        case 'extrCauseInfo'
            plotPlusMean(range, extrCauseInfo)
    end    
end
end

function plotPlusMean(range, var)
    hold on
    plot(range, var)
    plot(range, mean(var,1), '-k')
    xlim([1, max(range)])
end