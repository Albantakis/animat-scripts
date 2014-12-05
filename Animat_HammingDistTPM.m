clear all
condi = 'c1a2_36';
cd Hamming
load(strcat('HammingDist_', condi))
cd ..
load(strcat(condi, '_dataCB'));

%% "nice figure" Plot all trials in which 128 Fitness is reached 
[rowP0, colP0] = find(BigPhiMip == 0); %what's found here is put to NaN later CHANGE == to >

Phi = sparse(rowP0,colP0,1,50,118);

trials = find(Fitness_level(:,end) == 128);

F128 = Fitness_level(trials,:);
HammingDist128 = HammingDist(trials,:);

F128Phi= Fitness_level;%.*Phi; 
F128Phi(logical(Phi)) = NaN;
F128Phi = F128Phi(trials,:);
HammingDist128Phi =HammingDist; 
HammingDist128Phi(logical(Phi)) = NaN;
HammingDist128Phi = HammingDist128Phi(trials,:);

% for i = [1 10]
% figure
% subplot(2,1,1)
% hold on 
% plot(range, F128(i,:)./128*100, '-k')
% plot(range, F128Phi(i,:)./128*100, '-b')
% 
% subplot(2,1,2)
% hold on
% plot(range, HammingDist128(i,:), '-k')
% plot(range, HammingDist128Phi(i,:), '-b')
% end


figure
for i = 1:length(trials)
subplot(ceil(length(trials)/2),2,i)
hold on 
[AX, H1, H2] = plotyy(range, [F128(i,:); F128Phi(i,:)], range, [HammingDist128(i,:); HammingDist128Phi(i,:)]);

set(H1,'Color','k');
set(H2,'Color','k');
set(AX(1), 'YLim', [50, 128]);
set(AX(2), 'YLim', [0, 600]);
%[AX1, H3, H4] = plotyy(range, F128Phi(i,:), range, HammingDist128Phi(i,:))
%plot(range, HammingDist128(i,:), '-k')
%plot(range, HammingDist128Phi(i,:), '-b')
end

%% plot fitness of trials where 128 is reached and hamming distance for trials of 128
% Note the first 128 one should actually be excluded cause that has the hamming distance to the one before that was not 128
[row, col] = find(Fitness_level(:,end) == 128 & BigPhiMip(:,end) == 0); %CHANGE > to ==
col = 118;
F128Phi = sparse(row,col,1,50,118);
HammingF128Phi = F128Phi.*HammingDist;
HammingF128Phi(HammingF128Phi == 0) = NaN;

nnz(F128Phi)

% figure;
% subplot(2,1,1)
%     
%     plot(range, Fitness_level(trials,:)')
% subplot(2,1,2)
%     plot(range, HammingF128Phi(trials,:)')
%% find unique ones and plot them
DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
path = strcat(DPath, condi, '/trial');
F128Phivect = reshape(F128Phi',[],1); %those that are unique and have F = 128 and Phi =0/>0
Allntpm2D128Phi = Allntpm2D;
Allntpm2D128Phi(logical(~F128Phivect),:) = 0;
[UniqueF128Phi, ub, uc] = unique(Allntpm2D128Phi, 'rows');
ubsort = sort(ub);
ubsort = ubsort(ubsort >1);
indUniqueF128Phi = [floor((ubsort-1)/118)+1, mod(ubsort-1,118)+1];

cd Hamming
cd(condi) 

indrow = unique(indUniqueF128Phi(:,1));
count = 0;
Motors = zeros(length(indUniqueF128Phi),9216);
LODRange = zeros(length(indUniqueF128Phi),2);
for i = indrow'
    %load(strcat('AllNormalizedTPM_trial', int2str(i-1)))
    ind = find(indUniqueF128Phi(:,1) == i);
    for j = ind'
        count = count+1;
        %figure
        %imagesc(Allntpm(:,:,col2(j)))
        docname = strcat(path, int2str(i-1), '_', int2str(range(indUniqueF128Phi(j,2))), '_LifetimeLogicTable.txt');
        LogTable = importdata(docname,',', 1);
        MotorsTemp = LogTable.data(:, [end-1 end]);
        MotorsTemp(logical(sum(MotorsTemp,2) ==2),:) = 0;
        MotorsTemp = reshape(MotorsTemp, 1, []);
        Motors(count,:) = MotorsTemp; 
        LODRange(count,:) = [i, indUniqueF128Phi(j,2)];
    end
end
%%
[uM, aM, bM] = unique(Motors, 'rows');
size(uM, 1)
for i = 1:size(Motors,1)
    for j = 1:size(Motors,1)
        HammingMotor(i,j) = sum(Motors(i,:) ~= Motors(j,:));
    end
end
UPTHammingMotor = triu(HammingMotor);
UPTHammingMotor(UPTHammingMotor == 0) = NaN;
averageHammingMotor = nanmean(reshape(UPTHammingMotor, [],1))/9216
stdHammingMotor = nanstd(reshape(UPTHammingMotor, [],1))/9216
figure; hist(reshape(UPTHammingMotor, [],1)/9216, 50)
%% find unique ones and their distance from each other
F128Phivect = reshape(F128Phi',[],1); %those that are unique and have F = 128 and Phi =0/>0
Allntpm2D = Allntpm2D(logical(F128Phivect),:);
[Allntpm2D, ubtpm2D, uctpm2D] = unique(Allntpm2D, 'rows');
size(Allntpm2D,1)
HammingMat = zeros(size(Allntpm2D,1));
for i = 1:size(Allntpm2D,1)
    for j = 1:size(Allntpm2D,1)
        HammingMat(i,j) = sum(Allntpm2D(i,:) ~= Allntpm2D(j,:));
    end
end
UPTHamming = triu(HammingMat);
UPTHamming(UPTHamming == 0) = NaN;
averageHamming = nanmean(reshape(UPTHamming, [],1))/size(Allntpm2D,2)
stdHamming = nanstd(reshape(UPTHamming, [],1))/size(Allntpm2D,2)
figure; hist(reshape(UPTHamming, [],1)/size(Allntpm2D,2), 50)
cd ..
%% find unique connectivity matrices and their distances from each other
load(strcat('HammingDistConn_', condi))
F128Phivect = reshape(F128Phi',[],1); %those that are unique and have F = 128 and Phi =0/>0
Allncm2D = Allncm2D(logical(F128Phivect),:);
[Allncm2D, ubcm2D, uccm2D] = unique(Allncm2D, 'rows');
size(Allncm2D,1)
HammingMatConn = zeros(size(Allncm2D,1));
for i = 1:size(Allncm2D,1)
    for j = 1:size(Allncm2D,1)
        HammingMatConn(i,j) = sum(Allncm2D(i,:) ~= Allncm2D(j,:));
    end
end
UPTHammingConn = triu(HammingMatConn);
UPTHammingConn(UPTHammingConn == 0) = NaN;
averageHammingConn = nanmean(reshape(UPTHammingConn, [],1))/size(Allncm2D,2)
stdHammingConn = nanstd(reshape(UPTHammingConn, [],1))/size(Allncm2D,2)
figure; hist(reshape(UPTHammingConn, [],1)/size(Allncm2D,2), 50)

for i = 1:size(Allncm2D, 1)
    Conn = reshape(Allncm2D(i,:)',8,8)';
    J_Sparse = sparse(Conn);
    J = full(J_Sparse);
    indDiag = find(diag(J) == 1)
    view(biograph(J_Sparse));
end
%%
numNodes = 8;
numMot = 2;
numSen= 2;
sim = 1;
UniqueMotLODRange = LODRange(aM,:);

for c = 1:50%size(UniqueMotLODRange,1)
    t = UniqueMotLODRange(c,1);
    g = range(UniqueMotLODRange(c,2));
   
    J_tempfile = strcat(path, int2str(t-1),'_', int2str(g), '_EdgeList.txt');
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
        %view(biograph(J_sparse))
    end
    if sim == 1
    basic = load(strcat(DPath, condi, '/basic.txt'));
    for b = 1:length(basic)
        blockSize(b) = length(dec2bin(basic(b)));
    end    
    docname = strcat(path, int2str(t-1), '_', int2str(g), '_LifetimeLogicTable.txt');
    LogTable = importdata(docname,',', 1);
    Sensors = LogTable.data(:,1:2);
    Motors = LogTable.data(:, [end-1 end]);
    count = 0;
    count2 = 0;
    figure
    title(strcat(int2str(t), '__', int2str(g)))
    hold on
    for i = 1:length(basic)
        for j = -1:2:1
            for k = 1:16
                botPos = k;
                blockPos = 1;
                for l = 36:-1:1
                    count = count+1;
                    if k == 8
                        %set(gca,'NextPlot','replaceChildren');       
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
    end
end

%movie(F,20) % Play the movie twenty times