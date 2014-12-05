clear all
Elem = 0:7;
numNodes = 8;
numMot = 2;
numSen = 2;

%Larissa: Simflag was added to simulate the animat in a different task than it 
%evolved to and find whether the matching of that task was less
simflag = 1; 
basis = [1 2 1 2];

plotflag = 3;
trialnum = [0:49];
SizeLife = 4608;
%trialnum = setdiff([0:9],[0 2 6 8 9]);
%trialnum = [0 5 8 10 12 13 14 15 19];
numtrials = numel(trialnum);
totsteps = 60000-1;
step = 512;
range = 0:step:totsteps;
%range = [0:step:totsteps/2 29984 30208:step:totsteps 59984];
%startpoints = [696, 189, 349, 1692, 2356, 4372];
cond = 'c14a23_36';
%cond = 'c1a3_n001_rep20_36';
%cond = 'c1a3_n001_rep20_36';
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
% path3 = strcat(DPath, cond, '_2', '/trial');
% path2 = strcat(DPath, cond, '_2b', '/trial');
% cond2 = strcat(cond, '_2');
% cond3 = cond;
path2 = path;
path3 = path;
cond2 = cond;
cond3 = cond;

MaxFitness = 128;

BigPhiMipWorld = zeros(numtrials, length(range));
BigPhiMipNoise = zeros(numtrials,length(range));
NconceptWorld = zeros(numtrials, length(range));
NconceptNoise = zeros(numtrials,length(range));
Matching = zeros(numtrials,length(range));
PhiMatching = zeros(numtrials,length(range));
diffStates = zeros(numtrials,length(range));
SimFitness = zeros(numtrials, length(range));
UsedNodes = cell(numtrials,length(range));

ElemUsed = zeros(numtrials,numel(Elem),length(range));
ElemUsed_temp = zeros(1, numel(Elem));
FitPhiCorr = [];

%cd Results
for t = 1:numtrials
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
        %------------- get rest from Phi calculation --------------------------
        % only go to folder here, because I need to access functions from the
        % main folder then, so I need to go out of the folder after loading
        cd(strcat('Freeze_',cond, '_trial', int2str(trialnum(t))))
        Animat_gen = range(i);
        FilenameA = strcat('Animat_gen', int2str(Animat_gen),'_ShortResults_CorrBound_noSL.mat');
        FilenameB = strcat('Animat_gen', int2str(Animat_gen),'_NoiseResults_CorrBound_noSL.mat');
        if exist(FilenameA,'file') == 2 
            t 
            i
            load(FilenameB);
            Noise = results;
            load(FilenameA)
            cd ..
            MaxComplex = [];
            for com = 1:numel(Noise.indONS)
                MaxComplex = union(MaxComplex, Noise.Complex{com});
            end
            for com = 1:numel(results.LifeState_index)
                MaxComplex = union(MaxComplex, results.Complex{com});
            end
                       
            UsedNodes{t,i} = results.used_nodes;
            used_nodes = UsedNodes{t,i}+1;
            MaxComplex = MaxComplex(used_nodes(MaxComplex) > 2);
            
            %get tpm ------------------------------------------------------
%             docname = strcat(path, int2str(trialnum(t)),'_', int2str(range(i)), '_FullLogicTable.txt');
%             FullLogic = importdata(docname,',', 1);
%             FullLogic = FullLogic.data(:,[1:numNodes,(numNodes+2):(2*numNodes+1)]);
% 
%             indMot0 = 1:2^(numNodes-numMot); %Motors are always last in TPM
%             tpm = repmat(FullLogic(indMot0,numNodes+(1:numNodes)),2^numMot,1);
%             FullLogic = [FullLogic(:,1:numNodes), tpm];

            %Exclude nodes with only in or outputs (only Outputs can be dangerous)
%             excluded_nodes = setdiff(1:numNodes,used_nodes);
%             ind0 = find(sum(FullLogic(:,excluded_nodes),2) == 0);
%             RedStates = FullLogic(ind0,used_nodes); 
%             tpm = tpm(ind0,used_nodes);  
%             tpm(:,used_nodes < numSen+1) = 0.5; %Sensors might be switched on through mechanism but that is overwritten by environment --> doesn't do anything
            if simflag ==1
                [SimTransitions, SimTrans_index, p_SimTransitions, SimFitness(t,i)] = AnimatSimulation_MatchingDiffTask020714(range(i),t,path,8,used_nodes ,basis);
            end

            if ~isempty(MaxComplex) 
            network = NetworkCalculation(results.tpm, used_nodes, results.connect_mat, results.options);           
            if simflag == 1              
                Sensors = SimTransitions(:,used_nodes < 3);
                p_SimTransitions_temp = repmat(p_SimTransitions, 1,nnz(used_nodes < 3));
                pSenONtemp = sum(p_SimTransitions_temp.*Sensors);
                
                pSenON = zeros(numSen,1);
                pSenON(used_nodes(used_nodes < 3)) = pSenONtemp;
                
                [NoiseSimTransitions, p_NoiseSimTransitions] = AnimatSimulation_MatchingConst1216(results.tpm, pSenON, used_nodes, 4608); 
                
                NoiseSimTrans_index = [];
                for k = 1:size(NoiseSimTransitions,1)
                    NoiseSimTrans_index(k) = state2index(NoiseSimTransitions(k,:),2.*ones(size(NoiseSimTransitions,2),1));
                end 
                
                [DiffStates2, iLife, iNoise] = union(SimTrans_index, NoiseSimTrans_index);
                [DiffStates, ~, ~] = union(Noise.LifeTrans_index, Noise.NoiseTrans_index);
                SetD = setdiff(DiffStates2,DiffStates);
                if ~isempty(SetD)
                    SetD
                    t
                    i
                    DiffStates = [DiffStates SetD];
                end
                %Noise.LifeTrans_index = SimTrans_index;
            else
                [DiffStates, iLife, iNoise] = union(Noise.LifeTrans_index, Noise.NoiseTrans_index);
            end
            
            subsys_numel = zeros(numel(DiffStates),1);
            Phi = zeros(numel(DiffStates),1);
            Phi_MIP = zeros(numel(DiffStates),1);
            Num_Concepts = zeros(numel(DiffStates),1);
                        
            counter = 0;
            AllindState = 0;
            indConcepts = [];
            phiAll = [];
            CDistP = [];
            CDistF = [];
            countOtherStates = 0;
            for j = 1:numel(DiffStates)
                PastCurrState = index2state(DiffStates(j), 2.*ones(2*numel(used_nodes),1))';
                network.past_state = PastCurrState(1:numel(used_nodes));

                indState = find(Noise.LifeTrans_index == DiffStates(j));
                
                if ~isempty(indState)
                    subsys = results.Complex{indState};
                    sys_results = results.state(indState);
                    AllindState(j) = indState;
                else
                    indState = find(Noise.OnlyNoiseStates == DiffStates(j));
                    if ~isempty(indState)   %Empty can only happen for simflag == 1
                        subsys = Noise.Complex{indState};
                        sys_results = Noise.state(indState);
                        Add2index = numel(Noise.LifeTrans_index);
                        AllindState(j) = indState+Add2index;
                    else
                        countOtherStates = countOtherStates + 1;
                        network.current_state = PastCurrState(numel(used_nodes)+1:end);
                        [Big_phi_M phi_M prob_M M_cell concept_MIP_M purviews_M network concept_MIP_M_subs] = big_phi_all(network, network.current_state);                                             
                        [Big_phi_MIP MIP Complex M_i_max BFCut Big_phi_MIP_M complex_MIP_M Big_phi_MIP_all_M complex_MIP_M_all BFCut_M] = ...
                            complex_search(Big_phi_M,M_cell, purviews_M, network.num_nodes,prob_M,phi_M,network.options,concept_MIP_M,network);

                        subsys = Complex;
                        sys_results = Animat_rewrap_data_short(Big_phi_M, phi_M, prob_M, M_cell, concept_MIP_M, purviews_M,...
                                    Big_phi_MIP, MIP, Complex, M_i_max, BFCut);
                                
                        Add2index = numel(Noise.LifeTrans_index)+numel(Noise.indONS);
                        AllindState(j) = countOtherStates+Add2index;
                    end
                end    
                    
                Phi(j) = sys_results.Phi;
                % Larissa: This is because it seems that Complex by default
                % is 1 if Phi = 0
                if Phi(j) > 0
                    subsys_numel(j) = numel(subsys);
                else 
                    subsys_numel(j) = 0;
                end  


                subsys_ind = sys_results.main_complex;
                %sys_results = results.state(LSindex(j));

                Phi_MIP(j) = sys_results.Phi_MIP;
                IrrConcepts = [sys_results.concept(:).is_irreducible];
                Num_Concepts(j) = nnz(IrrConcepts);
                % higher level concepts
                numeratorindex = find(IrrConcepts);
                %Larissa: This could be done more efficiently counting the same
                %concepts and not calculating them all over again but for
                %the moment it works
                
                whole_sys_state = PastCurrState(numel(used_nodes)+1:2*numel(used_nodes));
                for nc = 1:numel(numeratorindex)
                    counter = counter+1;
                    denomP = sys_results.concept(numeratorindex(nc)).denominator.backwards;
                    denomF = sys_results.concept(numeratorindex(nc)).denominator.forwards;
                    [phi, CDistTempP, CDistTempF] = CalculateConcepts(subsys, index2concept(numeratorindex(nc), subsys), denomP, denomF, whole_sys_state,network,MaxComplex);
                    phiAll(counter) = phi;
                    indConcepts(counter) = j;
                    CDistP(counter,:) = CDistTempP;
                    CDistF(counter,:) = CDistTempF;
                end
            end
            
            if ~isempty(indConcepts)
            % Probability times phi for World and Noise
            PWorld = zeros(numel(DiffStates),1); PNoise = zeros(numel(DiffStates),1);
            
            if simflag == 1
                for ser = 1:numel(AllindState)
                    if AllindState(ser) <= numel(Noise.p_LifeTransitions)
                        ind = find(SimTrans_index == Noise.LifeTrans_index(AllindState(ser))); 
                        indN = find(NoiseSimTrans_index == Noise.LifeTrans_index(AllindState(ser))); 
                        if ~isempty(ind)
                            PWorld(ser) = p_SimTransitions(ind);
                        end
                        if ~isempty(indN)
                            PNoise(ser) = p_NoiseSimTransitions(indN);
                        end
                    else
                        if AllindState(ser) > (numel(Noise.p_LifeTransitions)+numel(Noise.indONS))
                            ind = find(SimTrans_index == DiffStates(AllindState(ser)));
                            indN = find(NoiseSimTrans_index == DiffStates(AllindState(ser)));
                        else
                            ind = find(SimTrans_index == Noise.NoiseTrans_index(Noise.indONS(AllindState(ser)-numel(Noise.p_LifeTransitions))));
                            indN = find(NoiseSimTrans_index == Noise.NoiseTrans_index(Noise.indONS(AllindState(ser)-numel(Noise.p_LifeTransitions)))); 
                        end   
                        if ~isempty(ind)
                            PWorld(ser) = p_SimTransitions(ind);
                        end
                        if ~isempty(indN)
                            PNoise(ser) = p_NoiseSimTransitions(indN);
                        end
                    end
                end  
            else
                PWorld(1:numel(results.p_LifeTransitions)) = results.p_LifeTransitions; 
                PWorld = PWorld(AllindState);
                for ser = 1:numel(AllindState)
                    if AllindState(ser) <= numel(Noise.p_LifeTransitions)
                        ind = find(Noise.NoiseTrans_index == Noise.LifeTrans_index(AllindState(ser)));
                        if ~isempty(ind)
                            PNoise(ser) = Noise.p_NoiseTransitions(ind);
                        end 
                    else
                        PNoise(ser) = Noise.p_NoiseTransitions(Noise.indONS(AllindState(ser)-numel(Noise.p_LifeTransitions)));
                    end
                end
            end                                  
           
%            [PhiSorted indsort] = sort(Phi_MIP);
% 
%             figure(t)
%             hold on
%             subplot(3,1,1)
%             h = bar(1:length(DiffStates), PhiSorted)
%             set(gca,'XTickLabel',[]);
% 
%             subplot(3,1,2)
%             hold on
%             h = bar(1:length(DiffStates), Num_Concepts(indsort))
%             set(gca,'XTickLabel',[]);
% 
%             subplot(3,1,3)
%             hold on
%             h = bar(1:length(DiffStates),PWorld(indsort)-PNoise(indsort));
%             set(h,'facecolor','red');
%             YLIM = [-0.05 0.2];
%             set(gca,'ylim',YLIM);
            
            %set(gca,'XTick',1:length(DiffStates))
            %set(gca,'XTickLabel',num2str(DiffStates(indsort)','%d')); % uncomment this to have a binary valued x-axis
            %rotateXLabels( gca(), 90)
            
            PWorldxphi = reshape(PWorld(indConcepts),[],1).*phiAll'; PNoisexphi = reshape(PNoise(indConcepts),[],1).*phiAll';
            PhiPWorldxphi = reshape(Phi_MIP(indConcepts),[],1).*reshape(PWorld(indConcepts),[],1).*phiAll'; PhiPNoisexphi = reshape(Phi_MIP(indConcepts),[],1).*reshape(PNoise(indConcepts),[],1).*phiAll'; %Also multiplied by PhiMIP
            
            BigPhiMipWorld(t,i) = sum(Phi_MIP.*PWorld);  
            BigPhiMipNoise(t,i) = sum(Phi_MIP.*PNoise); 
            
            NconceptWorld(t,i) = sum(Num_Concepts.*PWorld);  
            NconceptNoise(t,i) = sum(Num_Concepts.*PNoise);
            
            DPhi(t) = BigPhiMipWorld(t,i)-BigPhiMipNoise(t,i);
            
            % Distance Matrix
            back_maxent = expand_prob([],MaxComplex,[]);
            forward_maxent = comp_pers_cpt(network.nodes,[],MaxComplex,whole_sys_state,'forward',[]); %System state doesn't matter here
            forward_maxent = forward_maxent(:);
        
            
            [UniqueConcepts, A, indUnique] = unique([CDistP CDistF], 'rows');

            VWorld = zeros(size(UniqueConcepts,1),1); VNoise = zeros(size(UniqueConcepts,1),1);
            PhiVWorld = zeros(size(UniqueConcepts,1),1); PhiVNoise = zeros(size(UniqueConcepts,1),1);
            
            for s = 1:size(UniqueConcepts,1)
                VWorld(s) = VWorld(s) + sum(PWorldxphi(indUnique == s)); 
                VNoise(s) = VNoise(s) + sum(PNoisexphi(indUnique == s));
                
                PhiVWorld(s) = PhiVWorld(s) + sum(PhiPWorldxphi(indUnique == s)); 
                PhiVNoise(s) = PhiVNoise(s) + sum(PhiPNoisexphi(indUnique == s));

            end
            
            VDiff = sum(VWorld)-sum(VNoise);
            PhiVDiff = sum(PhiVWorld)-sum(PhiVNoise);
            
            if VDiff > 0
                VWorld = [VWorld; 0];
                VNoise = [VNoise; VDiff];
            else
                VWorld = [VWorld; -VDiff];
                VNoise = [VNoise; 0];
            end
            
            if PhiVDiff > 0
                PhiVWorld = [PhiVWorld; 0];
                PhiVNoise = [PhiVNoise; PhiVDiff];
            else
                PhiVWorld = [PhiVWorld; -PhiVDiff];
                PhiVNoise = [PhiVNoise; 0];
            end
            
            ConceptDist = reshape(UniqueConcepts, size(UniqueConcepts,1), [], 2);
            DistMat = Matching_genEMDDistanceMatrix(ConceptDist, [back_maxent,forward_maxent],network.gen_dist_matrix); %past whole and cut distributions      
            
            Matching(t,i) = emd_hat_gd_metric_mex(VWorld,VNoise,DistMat);
            PhiMatching(t,i) = emd_hat_gd_metric_mex(PhiVWorld,PhiVNoise,DistMat);
            end
            end
        else 
            cd ..
        end      
    end
end    
%cd ..
%%
if simflag == 1
    File = strcat(cond, '_Matching', '_', int2str(basis));
else
    File = strcat(cond, '_Matching');
end
save(File, 'BigPhiMipWorld', 'BigPhiMipNoise', 'Matching', 'PhiMatching', 'SimFitness','NconceptWorld', 'NconceptNoise', 'cond', 'range');

if plotflag == 3
    figure
    subplot(2,1,1)
    hold on
    plot(range, mean(BigPhiMipWorld,1), '-b')
    plot(range, mean(BigPhiMipNoise,1), '-r')
    xlim([1, max(range)])

    subplot(2,1,2)
    hold on
    plot(range, mean(Matching,1), '-k')
    plot(range, mean(PhiMatching,1), '-b')
    xlim([1, max(range)])
end