clear all
Elem = 0:7;
plotflag = 3;
trialnum = [0:49];
%trialnum = setdiff([0:9],[0 2 6 8 9]);
%trialnum = [0 5 8 10 12 13 14 15 19];
numtrials = numel(trialnum);
totsteps = 60000-1;
step = 512;
range = 0:step:totsteps;
%range = [0:step:totsteps/2 29984 30208:step:totsteps 59984];
%startpoints = [696, 189, 349, 1692, 2356, 4372];
%cond = 'c24a35_36';
%cond = 'c1a3_n001_rep20_36';
cond = 'c1a3_36'
%cond = 'c1a3_12Sen_36';
%cond = 'c1a3_n001_90000_36';
%cond = 'c1b11a2b5_36';
%cond = 'c23a45_36';
%cond = 'c35a271_36';
%cond = 'c56a47_36';
%cond = 'c1a3_12sen_36';
%DPath = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
path = strcat(DPath, cond, '/trial');
%  path3 = strcat(DPath2, cond,'_2', '/trial');
%  path2 = strcat(DPath, cond, '_2', '/trial');
%  cond2 = strcat(cond, '_2');
%  cond3 = cond2;
path2 = path;
path3 = path;
cond2 = cond;
cond3 = cond;

MaxFitness = 128;

for t =  [4 7 13 16 28 36];%1:numtrials [5 9 25 27 36 46 50];%
    figure(t)
        axis tight;
    count = 0;    
    for i = 1:length(range)
        %------------- get Fitness from Animat files---------------------------
        if trialnum(t) > 9
            cond = cond2;
            path = path2;
            if trialnum(t) > 19
            path = path3;
            cond = cond3;
            end
        end   
        docname2 = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_KOdata.txt');
        
        Fitness = load(docname2);
        Fitness_level(t,i) = Fitness(1);  
        %------------- get rest from Phi calculation --------------------------
        % only go to folder here, because I need to access functions from the
        % main folder then, so I need to go out of the folder after loading
        cd(strcat('Freeze_',cond, '_trial', int2str(trialnum(t))))
        Animat_gen = range(i);
        FilenameA = strcat('Animat_gen', int2str(Animat_gen),'_ShortResults_CorrBound_noSL.mat');
        if exist(FilenameA,'file') == 2
            if Fitness(1) == 128
            count = count+1;
            load(FilenameA)
            tpm = results.tpm;
            imagesc(tpm)           
            F(count) = getframe;
            end
            cd ..
        else 
            cd ..
        end
    end
    movie(F,1)
end
        