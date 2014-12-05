clear all
option = 0;
FA =  [0 128]; %[0 96 104 112 120 128];
meanVal = 0;
XX = 1:118;
factor = 5;

% Zombie Data
condition = 'c36a45_36';

Zombie = load(strcat(condition,'_ZombiedataAllC'));
load(strcat(condition,'_dataCB'));

%%
[PhiHistF, ind] = histc(reshape(Fitness_level, [],1), 64:4:128);
ZConcepts = reshape(Zombie.MeanNumConcepts, [], 1);
PhiMip = reshape(BigPhiMip, [], 1);
NConn = reshape(Num_Conn, [], 1);
figure 
boxplot(PhiMip, ind)
figure
boxplot(ZConcepts, ind)
% figure
% plot(NConn, PhiMip,'*k')
%%


for i = 1:17
figure(5)
hold on
    plot(NConn(ind == i), PhiMip(ind == i), '*k')
    xlim([1, 25])
    ylim([0,1.5])
end
