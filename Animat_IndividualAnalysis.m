%clear all
MaxFitness = 128;
totsteps = 60000-1;
op_Complex = true;
%condi = 'c1a3_n001_rep20';
condi = 'c1a3_change_c14a23';
%condi = 'c3a4';
%condi = 'c1a3_n001_rep20';
%condi = 'c36a45';
%cond = strcat(condi, '_36');
%cond = 'c35a271_36';
%cond = 'task7156331-3s'
cond = condi
%DPath = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
%DPath = '/Users/larissa_alb/Dropbox/larissa-jaime/mutualinformation/';
DPath = '/Users/larissa_alb/dev/animats/results/work_';
path = strcat(DPath, cond, '/trial');
path2 = path;
path3 = path;
cond2 = cond;
cond3 = cond;
path2 = strcat(DPath, cond, '/trial');

% path3 = strcat(DPath2, cond, '_2', '/trial');
% path2 = strcat(DPath, cond, '_2', '/trial');
% cond2 = strcat(cond, '_2');
% cond3 = cond2;

sim = 0;
numNodes = 8;
numMot = 2;
numSen= 2;
%trialnum = 0:49;
%trialnum = [193 90];
%trialnum = [0 2 4 8 12 21 33 40 43 44 49 52 61 70 83 86 87 88];
%trialnum = [3     4    15    24    27    35    37    45    49    50    51    53    56    59    67    72    76    81    82    83    84    89    90    95];%[0:49];
%trialnum = [3    53];
%trials = [ 1     4     5     6     7     9    11    15    17    19    28    29    30    32    33    36    39    41    42    48    49    50]; %[0:49]+1; %[4 5 9 13 16 25 27 28 36 38 46 50];%
%trialnum = [24 36 47 78 82 92 102 124 137 141 156 163 168]; trialnum = [206 215 232 260 264 292 331 346 349 363 390];
%%
if op_Complex
    load(strcat(cond, '_dataCB'));
end
Zombie = load(strcat(cond, '_ZombiedataAllC'));
range = Zombie.range;%./10000;
rangeS = 1:length(range); 
plotflag = 1;

trialnum = evaluatedTrials; 
trials = 1%1:length(trialnum);
%fitnessBeforeChange = max(Fitness_level(:,55:59),[],2);
%trials = trials(find(max(BigPhiMip(:,55:59),[],2) == 0 & fitnessBeforeChange == 120));
%trials = trials(find(max(BigPhiMip(:,59),[],2) == 0));
numel(trials)
%%
% Plot individual trials
if plotflag == 1;
for i = trials;%[0,2,3,8,10,11,12,13,16,18]+1
    figure(1000+i)
    subplot(3,2,1)
        hold on
        %plot(range, 100.*Fitness_level(i,rangeS)./128, 'k')
        plot(range, Fitness_level(i,rangeS), 'Color', 'k', 'DisplayName', int2str(trialnum(i)))
        %plot(range, 100.*mean(Fitness_level(:,rangeS))./128, 'color', [.6, .6, .6])
        plot(range, mean(Fitness_level(:,rangeS)), 'color', [.6, .6, .6])
        
        
        indFitPhi = find(BigPhiMip(i,rangeS) > 0);
        %plot(range(indFitPhi), Fitness_level(i,indFitPhi), '*r')
        
        xlim([0, max(range)])
        ylim([50,128])
    subplot(3,2,2)
        hold on
        plot(range, MeanSizeComplex(i,rangeS), 'k')
        plot(range, mean(MeanSizeComplex(:,rangeS)), 'color', [.6, .6, .6])
        xlim([0, max(range)])
    subplot(3,2,3)
        hold on 
        %plot(range, BigPhiMip(i,rangeS),  'k')
        %plot(range, Zombie2.BigPhi(i,:), 'b')
        plot(range, Zombie.BigPhi(i,:), 'r')
        xlim([0, max(range)])
        hold on
        %plot(range, mean(BigPhiMip(:,rangeS)), 'color', [.6, .6, .6])
        plot(range, mean(Zombie.BigPhi), 'color', [1, .5, .5])
    subplot(3,2,4)
        hold on 
        plot(range, BigPhiMip(i,rangeS),  'k')
        plot(range, BigPhiMip_max(i,rangeS), 'm')
        %plot(range, Zombie.BigPhi(i,:), 'r')
        xlim([0, max(range)])
        hold on
        plot(range, mean(BigPhiMip(:,rangeS)), 'color', [.6, .6, .6])
        %plot(range, mean(Zombie.BigPhi), 'color', [1, .5, .5])
    subplot(3,2,5)    
        hold on
        plot(range, Num_Conn(i,rangeS),'k')
        plot(range, mean(Num_Conn(:,rangeS)), 'color', [.6, .6, .6])
        %plot(range, MeanNumConcepts(i,rangeS),'k')
        %plot(range, mean(MeanNumConcepts(:,rangeS)), 'color', [.6, .6, .6])
        %plot(range, Zombie2.MeanNumConcepts(i,:),'b')
        plot(range, Zombie.MeanNumConcepts(i,:),'r')
        plot(range, mean(Zombie.MeanNumConcepts), 'color', [1, .5, .5])
        xlim([0, max(range)])
    subplot(3,2,6)    
        hold on
        plot(range, MeanNumConcepts(i,rangeS),'k')
        plot(range, mean(MeanNumConcepts(:,rangeS)), 'color', [.6, .6, .6])
        %plot(range, Zombie.MeanNumConcepts(i,:),'r')
        %plot(range, mean(Zombie.MeanNumConcepts), 'color', [1, .5, .5])
        xlim([0, max(range)])
end    
end
%%
range = (0:512:totsteps)/10000;
for t = trials;
gen = range(59)*10000;
%gen = range([59:60 end])*10000;
for g = gen;
    %------------- get Fitness from Animat files---------------------------
    if trialnum(t) > 99
        path = path2;
    end     
    
    J_tempfile = strcat(path, int2str(trialnum(t)),'_', int2str(g), '_EdgeList.txt')
    J_temp = load(J_tempfile);

    if ~isempty(J_temp)
        t
        J_temp = unique(J_temp, 'rows')+1;
        % MOTORS ARE SET TO 0 -> they don't actually have recurrent
        % connections
        J_temp = J_temp(J_temp(:,1) <= numNodes-numMot,:); 
        %Sensors might have incoming connection, but because they actually don't do anything they shouldn't be taken into account
        J_temp = J_temp(J_temp(:,2) > numSen,:);
        J_temp(J_temp(:,1) == J_temp(:,2),:)
        J_sparse = sparse(J_temp(:,1), J_temp(:,2),1,numNodes,numNodes);
        view(biograph(J_sparse))
    end
    if sim > 0
    basic = load(strcat(DPath, cond, '/basic.txt'));
    for b = 1:length(basic)
        blockSize(b) = length(dec2bin(basic(b)));
    end    
    docname = strcat(path, int2str(trialnum(t)), '_', int2str(gen), '_LifetimeLogicTable.txt');
    LogTable = importdata(docname,',', 1);
    Sensors = LogTable.data(:,1:2);
    Motors = LogTable.data(:, [end-1 end]);
    
    if sim == 1
        count = 0;
        count2 = 0;
        %figure
        for i = 1:length(basic)
            for j = -1:2:1
                for k = 1:16
                    botPos = k;
                    blockPos = 1;
                    for l = 36:-1:1
                        count = count+1;
                        if k == 8
                            %set(gca,'NextPlot','replaceChildren');       
                            figure(3000+trialnum(t))
                            m = 5+i*j;
                            if m >4 
                                m = m-1;
                            end    
                            subplot(2,4,m)
                            count2 = count2+1;
                            rectangle('Position',[blockPos,l,blockSize(i),1],'Curvature',[0,0],'FaceColor','b')
                            rectangle('Position',[botPos,-l,3,1],'Curvature',[0,0],'FaceColor','y')
                            if Sensors(count,1) == 1
                                rectangle('Position',[botPos,-l,1,1],'Curvature',[0,0],'FaceColor','r')
                            end    
                            if Sensors(count,2) == 1
                                rectangle('Position',[botPos+2,-l,1,1],'Curvature',[0,0],'FaceColor','r')
                            end    
                            xlim([0,18])
                            ylim([-36,37])
                            %F(count2) = getframe;
                            %clf;
                        end
                        if j == -1
                            blockPos = mod(blockPos - 1, 16);
                        else    
                            blockPos = mod(blockPos + 1, 16);
                        end    
                        botPos = botPos + Motors(count,1) - Motors(count,2);
                        botPos = mod(botPos, 16);
                    end
                end
                %close(1000)
            end
        end
    else 
        LifeStates = [LogTable.data(2:end,1:6) LogTable.data(1:end-1, [end-1,end])];
        figure
        x = 16;
        count = 0;
        for l = (x*36+1):(x*36+36) %length(LifeStates)
            count = count+1;
            currState = ([1 2 2 3; 1 2 2 3].*reshape(LifeStates(l,:),2,[]))';
            set(gca,'NextPlot','replaceChildren');
            imagesc(currState)
            caxis([0 4])
            F(count) = getframe;
            clf;
        end
    end
    end
end
end
%movie(F,5,2) % Play the movie twenty times

%%
indMaxF = find(Fitness_level(:,end)==128)-1
numel(indMaxF)
indMaxFPhhi = indMaxF(find(BigPhiMip(indMaxF+1,end)>0))
numel(indMaxFPhhi)
indFF = find(BigPhiMip(:,end)==0)-1
numel(indFF)
indFPhi = find(BigPhiMip(:,end)>0)-1
mean(Fitness_level(indFPhi+1,end))
mean(Fitness_level(indFF+1,end))

%% Correlation Histograms
%rangeS = 1:49;
%mean correlation coefficients
F = mean(Fitness_level(:,rangeS));
P = mean(BigPhiMip(:,rangeS));
EC = mean(MeanSizeComplex(:,rangeS));
C = mean(MeanNumConcepts(:,rangeS));
ZP = mean(Zombie.BigPhi(:,rangeS));
ZC = mean(Zombie.MeanNumConcepts(:,rangeS));
%CapME = mean(results.Capture_maxEnt);
%CapInd = mean(results.Capture_Gen1st);

%CorrMat = corrcoef([F; P; C; ZP; ZC; CapME; CapInd]');
[CorrMat, pCorr] = corrcoef([F; P; C; EC; ZP; ZC]');
CorrMat(1,:)
pCorr(1,:)

for i = 1:size(Fitness_level,1)
    F = Fitness_level(i,rangeS);
    P = BigPhiMip(i,rangeS);
    EC = MeanSizeComplex(i,rangeS);
    C = MeanNumConcepts(i,rangeS);
    ZP = Zombie.BigPhi(i,rangeS);
    ZC = Zombie.MeanNumConcepts(i,rangeS);
    %CapME = results.Capture_Gen1st(i,:);
    %CapInd = results.Capture_Change(i,:);
    %CorrMat_trial = corrcoef([F; P; C; ZP; ZC; CapME; CapInd]');
    
    [CorrMat_trial pCorrMat_trial] = corrcoef([F; P; C; EC; ZP; ZC]');
    CorrMat_trial(isnan(CorrMat_trial)) = 0;
    FitCorr(i,:) = CorrMat_trial(1,:);
    pFitCorr(i,:) = pCorrMat_trial(1,:);
    BPhiCorr(i,:) = CorrMat_trial(4,:);
end

figure(1000)
s(1) = subplot(5,1,1);
    hold on
    hist(FitCorr(:,2), -1:0.1:1);
    hist(FitCorr(pFitCorr(:,2) <= 0.05,2), -1:0.1:1);
    line([mean(FitCorr(:,2)),mean(FitCorr(:,2))],[0,20]);
    xlim([-1.05,1.05])  
s(2) = subplot(5,1,3)   ; 
    hold on
    hist(FitCorr(:,5), -1:0.1:1);
    hist(FitCorr(pFitCorr(:,5) <= 0.05,5), -1:0.1:1);
    line([mean(FitCorr(:,5)),mean(FitCorr(:,5))],[0,20]);
    xlim([-1.05,1.05])
s(3) = subplot(5,1,2);
    hold on
    hist(FitCorr(:,3), -1:0.1:1);
    hist(FitCorr(pFitCorr(:,3) <= 0.05,3), -1:0.1:1);
    line([mean(FitCorr(:,3)),mean(FitCorr(:,3))],[0,10]);
    xlim([-1.05,1.05])
s(4) = subplot(5,1,4);
    hold on
    hist(FitCorr(:,6), -1:0.1:1);
    hist(FitCorr(pFitCorr(:,6) <= 0.05,6), -1:0.1:1);
    line([mean(FitCorr(:,6)),mean(FitCorr(:,6))],[0,10]);
    xlim([-1.05,1.05])
s(5) = subplot(5,1,5);
    hold on
    hist(FitCorr(:,4), -1:0.1:1);
    hist(FitCorr(pFitCorr(:,4) <= 0.05,4), -1:0.1:1);
    line([mean(FitCorr(:,4)),mean(FitCorr(:,4))],[0,10]);
    xlim([-1.05,1.05])
% s(6) = subplot(4,2,7);
%     hold on
%     hist(FitCorr(:,7), -1:0.1:1);
%     line([mean(FitCorr(:,7)),mean(FitCorr(:,7))],[0,10]);
%     xlim([-1.05,1.05])
% s(7) = subplot(4,2,6);
%     hold on
%     hist(BPhiCorr(:,6), -1:0.1:1);
%     line([mean(BPhiCorr(:,6)),mean(BPhiCorr(:,6))],[0,10]); 
%     xlim([-1.05,1.05])
% s(8) = subplot(4,2,8);
%     hold on
%     hist(BPhiCorr(:,7), -1:0.1:1);
%     line([mean(BPhiCorr(:,7)),mean(BPhiCorr(:,7))],[0,10]);
%     xlim([-1.05,1.05])
    
title(s(1),'PhiMip')
title(s(2),'ZombiePhi')
title(s(3),'#Concepts')
title(s(4),'#ZombieConcepts')
title(s(5),'#Elements in Complex')

% title(s(5),'Capture ME')
% title(s(6),'Capture 1O')
% title(s(7),'Capture ME vs PhiMip')
% title(s(8),'Capture 1O vs PhiMip')

%%
Nodesexist = zeros(50,118);
for i = 1:size(UsedNodes,1)
    for j = 1:size(UsedNodes,2)
    Nodesexist(i,j) = nnz((UsedNodes{i,j} > 1).*(UsedNodes{i,j} < 6));
    end
end

%% Plot Phi vs Fitness

%set(0,'DefaultAxesColorOrder',[0:1/118:1; 1:-1/118:0; ones(1,119)]')
figure(8)
    hold on
    for j = trials'
        for i = 1:length(range)
            %line([Fitness_level(j,i) Fitness_level(j,i+1)], [BigPhi(j,i) BigPhi(j,i+1)], [range(i) range(i+1)], 'Color', [(i+1)/length(range) 0.5 1-(i+1)/length(range)])
            %line([Fitness_level(j,i) Fitness_level(j,i+1)], [BigPhiMip(j,i) BigPhiMip(j,i+1)], [range(i) range(i+1)], 'Color', [(i+1)/length(range) 0.5 1-(i+1)/length(range)])
            %plot(Fitness_level(j,i), Zombie.BigPhi(j,i), '.', 'Color', [i/length(range) 0.5*(1-i/length(range)) 1-i/length(range)],'MarkerSize', 15);
            plot(Fitness_level(j,i), Phi2008(j,i), '.', 'Color', [i/length(range) 0.5*(1-i/length(range)) 1-i/length(range)],'MarkerSize', 15);
        end
    end
xlim([56 128])

%% Plot 3D - Check whether increase in complexity is prior to increase in fitness.
FitnessDiff = [zeros(size(Fitness_level,1),1) diff(Fitness_level,1,2)];

%set(0,'DefaultAxesColorOrder',[0:1/118:1; 1:-1/118:0; ones(1,119)]')
%figure
    %hold on
    %plot3(Fitness_level, BigPhi, repmat(range,20,1))
    for j = trials
            figure
        for i = 1:length(range)-1
            %line([Fitness_level(j,i) Fitness_level(j,i+1)], [BigPhi(j,i) BigPhi(j,i+1)], [range(i) range(i+1)], 'Color', [(i+1)/length(range) 0.5 1-(i+1)/length(range)])
            %line([Fitness_level(j,i) Fitness_level(j,i+1)], [BigPhiMip(j,i) BigPhiMip(j,i+1)], [range(i) range(i+1)], 'Color', [(i+1)/length(range) 0.5 1-(i+1)/length(range)])
            line([Fitness_level(j,i) Fitness_level(j,i+1)], [FitnessDiff(j,i) FitnessDiff(j,i+1)],[BigPhiMip(j,i) BigPhiMip(j,i+1)], 'Color', [(i+1)/length(range) 0.5 1-(i+1)/length(range)])

            %line([Fitness_level(j,i) Fitness_level(j,i+1)], [BigPhiMip(j,i) BigPhiMip(j,i+1)], [range(i) range(i+1)], 'Color', [BigPhiMip(j,i)/max(max(BigPhiMip)) 1-BigPhiMip(j,i)/max(max(BigPhiMip)) 1-BigPhiMip(j,i)/max(max(BigPhiMip))])
            %line([Fitness_level(j,i) Fitness_level(j,i+1)], [MeanNumConcepts(j,i) MeanNumConcepts(j,i+1)], [range(i) range(i+1)], 'Color', [(i+1)/length(range) 0.5 1-(i+1)/length(range)])
        end
    end
%% Plot differences in Fitness    
FitnessDiff = diff(Fitness_level,1,2);
figure; plot(FitnessDiff', '-x')
FitnessDiff3D = zeros(size(FitnessDiff,1), size(FitnessDiff,2), size(FitnessDiff,2));
LogicBP = BigPhiMip > 0;
for i = 1:size(FitnessDiff,1)
    for j = 1:size(FitnessDiff,2)
        if LogicBP(i,j)
            FitnessDiff3D(i,j,1) = FitnessDiff(i,j);
        else
            FitnessDiff3D(i,1,j) = FitnessDiff(i,j);
        end
    end
end
figure
hold on
for i = 1:size(FitnessDiff,1)
    plot3(range(2:end), squeeze(FitnessDiff3D(i,:,1)),squeeze(FitnessDiff3D(i,1,:)));
end    
%%
figure
hold on
Fitness_levelP = Fitness_level;
Fitness_levelP(BigPhiMip == 0) = NaN;
for i = trials
    plot3(BigPhiMip(i,:), i.*ones(size(Fitness_level,2),1), Fitness_level(i,:)) 
    plot3(BigPhiMip(i,:), i.*ones(size(Fitness_level,2),1), Fitness_levelP(i,:), 'r') 
end

%%
figure
hold on
Fitness_levelP = Fitness_level;
Fitness_levelP(BigPhiMip == 0) = NaN;
for i = trials
    plot3(range, i.*ones(size(Fitness_level,2),1), Fitness_level(i,:)) 
    plot3(range, i.*ones(size(Fitness_level,2),1), Fitness_levelP(i,:), 'r') 
end
%% Make histogram of Fitness jumps
FitnessDiff(FitnessDiff == 0) = NaN;
FitnessChangePhi = FitnessDiff;
FitnessChangePhi0 = FitnessDiff;
FitnessChangePhi(BigPhiMip(:,1:end-1) == 0) = NaN;
FitnessChangePhi0(BigPhiMip(:,1:end-1) > 0) = NaN;

FCPhi = reshape(FitnessChangePhi,[],1);
FCPhi0 = reshape(FitnessChangePhi0,[],1);
nnz(~isnan(FCPhi))
nnz(~isnan(FCPhi0))
nnz(FCPhi > 0)
nnz(FCPhi0 > 0)

figure 
hold on 
A = histc(FCPhi, -40:1:60);
B = histc(FCPhi0, -40:1:60);
bar(-40:1:60,[A B])
%% find unique values of sum(smallphi) for a given Fitness and make histogram
F128 = reshape(Fitness_level,[],1) == 128;
Phi0 = reshape(BigPhiMip,[],1) == 0;
BigPhi128 = reshape(Zombie.BigPhi, [],1);
BigPhi128_phi0 = BigPhi128(F128.*Phi0==1);
BigPhi128_phi = BigPhi128(F128.*(1-Phi0)==1);
[un0, row0, col0] = unique(BigPhi128_phi0);
[unp, rowp, colp] = unique(BigPhi128_phi);
figure; 
subplot(2,1,1)
    hist(un0, 0:0.01:2)
subplot(2,1,2)
    hist(unp, 0:0.01:2)
    

