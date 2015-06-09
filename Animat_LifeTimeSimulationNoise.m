function AnimatLifeTimeSimulationNoise

MaxFitness = 128;
totsteps = 60000-1;
step = 512;
range = 0:step:totsteps;
trialnum = [0:49];
numtrials = numel(trialnum);

%condi = 'c1a3_n001_rep20';
condi = 'c1a3';
cond = strcat(condi, '_36');
DPath2 = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
path = strcat(DPath, cond, '/trial');
% path2 = path;
% path3 = path;
% cond2 = cond;
% cond3 = cond;
path3 = strcat(DPath2, cond, '_2', '/trial');
path2 = strcat(DPath, cond, '_2', '/trial');
cond2 = strcat(cond, '_2');
cond3 = cond2;

AllCorr = zeros(numtrials,length(range));

basic = load(strcat(DPath, cond, '/basic.txt'));

for b = 1:length(basic)
    blockSize(b) = length(dec2bin(basic(b)));
end    

%cd Results
for t = 1:numtrials
    for i = length(range)
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
            load(FilenameA)
            cd ..         
            used_nodes = results.used_nodes+1;
            tpm = results.tpm;
            
            %figure
            for n = 1:50            
                count = 0;
                correct = 0;
            for x = 1:length(basic)
                for j = -1:2:1
                    for k = 1:16
                        agent_states = zeros(1, size(tpm, 2));
                        botPos = k-1;
                        blockPos = 0;
                        for l = 36:-1:1
                            count = count+1;
                            
                            sensor = zeros(1,2);
                            sensor(1)= ~isempty(intersect(botPos, mod(blockPos:blockPos+blockSize(x)-1,16))); 
                            sensor(2)= ~isempty(intersect(mod(botPos+2,16), mod(blockPos:blockPos+blockSize(x)-1,16)));
                            
                            Noise = rand(2,1) <= [0.01; 0.01];
                            sensor = abs(sensor - Noise');
                            
                            agent_states(used_nodes < 3) = sensor(used_nodes(used_nodes < 3));
                            
                            %applyNoise(agent, sensorNoise);
                    
                            agent_states(used_nodes > 6)=0; 
                            agent_states = updateStates(agent_states, tpm);
                    
                            action= sum([agent_states(used_nodes == 7), 2*agent_states(used_nodes == 8)]);
                            switch action
                            case 0
                            case 3
                            case 1
                                botPos=(botPos+1);
                            case 2
                                botPos=(botPos-1);
                            end    
                            botPos = mod(botPos, 16);

                            if j == -1
                                blockPos = mod(blockPos - 1, 16);
                            else    
                                blockPos = mod(blockPos + 1, 16);
                            end  
                            %LifeState(count,:) = agent_states;
                        end
                        
                        hit=0;
                        if ~isempty(intersect(mod(botPos:botPos+2,16), mod(blockPos:blockPos+blockSize(x)-1,16)))
                            hit = 1;
                        end    
                        
                        if mod(x,2) == 0
                            if hit == 0
                                correct = correct +1;
                            end
                        elseif mod(x,2) == 1
                            if hit == 1
                                correct = correct +1;
                            end
                        end
                        
                    end
                           
                end
            end
            Ncorr(n) = correct;
            end
        else
            cd ..
        end
    end
    AllCorr(t) = mean(Ncorr);
end
end

function new_state = updateStates(agent_state, tpm)
    indSt = state2index(agent_state, 2.*ones(size(agent_state,2),1));
    new_state = tpm(indSt,:);
end

function index = state2index(state_vec, num_states_vec)
index = state_vec(1) + 1;
for i = 2:length(num_states_vec)
    index = index + state_vec(i)*prod(num_states_vec(1:i-1));  
end
end
    