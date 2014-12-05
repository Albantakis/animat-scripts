condition = 'c14a23_36';
load(strcat(condition, '_dataCB'))
Sim = load(strcat(condition,'_Matching_1212'));
Match = load(strcat(condition,'_Matching'));
%%
for i = 20%1:50
    figure(100+i)
    subplot(3,1,1)
    hold on
    plot(Fitness_level(i,:), '-b')
    plot(Sim.SimFitness(i,:), '-r')
    subplot(3,1,2)
    hold on
    plot(Match.BigPhiMipWorld(i,:)-Match.BigPhiMipNoise(i,:), '-b')
    plot(Sim.BigPhiMipWorld(i,:)-Sim.BigPhiMipNoise(i,:), '-r')
    subplot(3,1,3)
    hold on
    plot(Match.PhiMatching(i,:),'-b')
    plot(Sim.PhiMatching(i,:), '-r')
end
%%
figure
hold on
%plot((Match.PhiMatching-Sim.PhiMatching)')
hist((Match.PhiMatching(:,end)-Sim.PhiMatching(:,end))', -1.5:0.01:1)