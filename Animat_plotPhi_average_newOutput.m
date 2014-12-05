clear all
Elem = 0:7;
plotflag = 3;
trialnum = [3 4 15 24 27 35 37 45 49 50 51 53 56 59 67 72 76 81 82 83 84 89 90 95 100 103 104 118 130 131 132 136 145 147 157 163 170 175 182 186 187 191 193 195 197];
%trialnum = [0 3 4 5 6 10 18 27 28 29 31 32 35 38 47 48 49 52 57 58 61 62 64 65 66 68 70 72 74 76 79 83 85 91 92 96 100 110 111 118 125 129 132 133 134 135 141 143 150 152 153 157 158 162 164 167 171 174 180 185 186 192 199];
%trialnum = [0 2 4 7 8 10 12 19 21 23 33 40 43 44 45 48 49 50 52 54 61 70 83 86 87 88 93 104 105 106 108 112 121 124 126 130 131 133 146 148 149 154 162 164 177 183 195 198];

numtrials = numel(trialnum);
totsteps = 60000-1;
step = 512;
range = 0:step:totsteps;
%range = [0:step:totsteps/2 29984 30208:step:totsteps 59984];
%startpoints = [696, 189, 349, 1692, 2356, 4372];
%cond = 'c24a35_36';
%cond = 'c1a3_n001_rep20_36';
%cond = 'c1a3_n001_rep20_36'
%cond = 'c1a3_12Sen_36';
%cond = 'c1a3_n001_90000_36';
%cond = 'c1b11a2b5_36';
%cond = 'c23a45_36';
%cond = 'c35a271_36';
cond = 'c1a3_change_c36a45';
%cond = 'c1a3_12sen_36';
DPath = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
DPath2 = '~/dev/animats/results/work_';
path = strcat(DPath, cond, '/trial');
% path3 = strcat(DPath2, cond,'_2', '/trial');
% path2 = strcat(DPath, cond, '_2', '/trial');
% cond2 = strcat(cond, '_2');
% cond3 = cond2;
%path2 = path;
cond2 = cond;
path2 = strcat(DPath2, cond,'_100', '/trial');

MaxFitness = 128;

BigPhi = zeros(numtrials, length(range));
BigPhiMip = zeros(numtrials,length(range));
BigPhiMip_B = zeros(numtrials,length(range));
BigPhiMip_max = zeros(numtrials,length(range));
diffStates = zeros(numtrials,length(range));
NumLifeStates = zeros(numtrials,length(range));
MeanNumConcepts = zeros(numtrials,length(range));
MeanSizeComplex = zeros(numtrials,length(range));
MeanHOConcepts = zeros(numtrials,length(range));
Fitness_level = zeros(numtrials,length(range));
Num_Conn = zeros(numtrials,length(range));

UsedNodes = cell(numtrials,length(range));

ElemUsed = zeros(numtrials,numel(Elem),length(range));
ElemUsed_temp = zeros(1, numel(Elem));
FitPhiCorr = [];

%cd Results
for t = 1:numtrials
    for i = 1:length(range)
        %------------- get Fitness from Animat files---------------------------
        if trialnum(t) > 99
            path = path2;
        end   
        %docname = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_EdgeList.txt');
        docname2 = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_KOdata.txt');
        
        Fitness = load(docname2);
        Fitness_level(t,i) = Fitness(1);  
%         EdgeList = load(docname);
%         if ~isempty(EdgeList)
%             Elements = unique(EdgeList); 
%             ElemUsed_temp(Elements+1) = 5;
%             Input = unique(EdgeList(:,1));
%             Output = unique(EdgeList(:,2));
%             ElemUsed_temp(setdiff(Elements, Output)+1) = 2;
%             ElemUsed_temp(setdiff(Elements, Input)+1) = 1;
%             ElemUsed(t,:,i) = ElemUsed_temp;
%         end
      
        %------------- get rest from Phi calculation --------------------------
        % only go to folder here, because I need to access functions from the
        % main folder then, so I need to go out of the folder after loading
        cd(strcat('Freeze_',cond, '_trial', int2str(trialnum(t))))
        Animat_gen = range(i);
        FilenameA = strcat('Animat_gen', int2str(Animat_gen),'_ShortResults_CorrBound_noSL.mat');
        FilenameB = strcat('Animat_gen', int2str(Animat_gen),'_Results.mat');
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
            UsedNodes{t,i} = results.used_nodes;
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


            BigPhi(t,i) = sum(results.p_LifeTransitions .* Phi); %weighted ave/expected value
            BigPhiMip(t,i) = sum(results.p_LifeTransitions .* Phi_MIP);  
            BigPhiMip_B(t,i) = mean(Phi_MIP);  
            BigPhiMip_max(t,i) = max(Phi_MIP);  

            diffStates(t,i) = length(unique(Phi_MIP));
            NumLifeStates(t,i) = size(results.p_LifeTransitions, 1);
            Num_Conn(t,i) = results.numConn;
            MeanNumConcepts(t,i) = mean(Num_Concepts);
            MeanSizeComplex(t,i) = mean(subsys_numel);
            MeanHOConcepts(t,i) = mean(MConceptOrder);
            FitPhiCorr = [FitPhiCorr; [Fitness(1), BigPhiMip(t,i), range(i)]];
        else 
            cd ..
        end      
    end
end    
%cd ..
%%
Fitness_Phi_Corr = corrcoef(mean(Fitness_level),mean(BigPhiMip))
MeanFitness = mean(Fitness_level,1);
MeanFitness = 100.*MeanFitness./MaxFitness;

File = strcat(cond, '_dataCB');
save(File, 'BigPhi','BigPhiMip','BigPhiMip_B','BigPhiMip_max','UsedNodes','Fitness_Phi_Corr', 'Fitness_level', 'MeanFitness', 'MeanNumConcepts', 'MeanHOConcepts', 'MeanSizeComplex', 'NumLifeStates', 'Num_Conn', 'cond', 'range');

if plotflag == 3
    figure
    subplot(3,1,1)
    hold on
    plot(range, MeanFitness)
    xlim([1, totsteps])
    
    subplot(3,1,2)
    hold on
    plot(range, mean(BigPhi,1), '-b')
    plot(range, mean(BigPhiMip_B,1), '-g')
    plot(range, mean(BigPhiMip_max,1), '-m')
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

% 
%  mBigPhi=[]; mBigPhiMip_B=[]; mBigPhiMip_max=[]; mBigPhiMip=[]; mNumLifeStates=[]; mdiffStates=[];
% for s = 0:50:totsteps-50 
%     mBigPhi = [mBigPhi mean(sum(BigPhi(:,s+1:s+50),2))];
%     mBigPhiMip_B = [mBigPhiMip_B mean(sum(BigPhiMip_B(:,s+1:s+50),2))];
%     mBigPhiMip_max = [mBigPhiMip_max mean(sum(BigPhiMip_max(:,s+1:s+50),2))];
%     mBigPhiMip = [mBigPhiMip mean(sum(BigPhiMip(:,s+1:s+50),2))];
%     mNumLifeStates = [mNumLifeStates mean(sum(NumLifeStates(:,s+1:s+50),2))];
%     mdiffStates = [mdiffStates mean(sum(diffStates(:,s+1:s+50),2))];
% end
% 
% if plotflag == 2
%     figure
%     subplot(3,1,1)
%     hold on
%     plot(mean(Fitness_level,1))
%     xlim([1, totsteps])
%     
%     subplot(3,1,2)
%     hold on
%     plot(25:50:totsteps, mBigPhi, '.-b')
%     plot(25:50:totsteps, mBigPhiMip_B, '.-g')
%     plot(25:50:totsteps, mBigPhiMip_max, '.-m')
%     plot(25:50:totsteps, mBigPhiMip, '.-r')
%     xlim([1, totsteps])
% 
%     subplot(3,1,3)
%     hold on
%     %plot(range, NumConcept, '.-b')
%     plot(25:50:totsteps, mNumLifeStates, '.-b')
%     plot(25:50:totsteps, mdiffStates, '.-r')
%     xlim([1, totsteps])
% end
% 
% 
if plotflag == 2
    figure
    subplot(2,1,1)
    hold on
    plot(FitPhiCorr(:,1), FitPhiCorr(:,2), '.k')
    
    [Fitind, Frow, Fcol] = unique(FitPhiCorr(:,1));
    for i= 1:length(Fitind)
        FPhi(i) = mean(FitPhiCorr(Fcol == i,2));
    end
    
    subplot(2,1,2)
    hold on
    plot(Fitind, FPhi)
end    
