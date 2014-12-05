lim = 99:118;
inc = 27
Delta = Ent_NoiseComplex - Ent_LifeComplex;
DeltaIn = Ent_NoiseInput - Ent_LifeInput;

figure
subplot(1,3,1)
hold on
plot(mean(BigPhiMip(:,lim),2), mean(Delta(:,lim),2), '*k');
plot(mean(BigPhiMip(inc,lim),2), mean(Delta(inc,lim),2), '*r');


subplot(1,3,2)
hold on
plot(mean(MeanNumConcepts(:,lim),2), mean(Delta(:,lim),2), '*k');
plot(mean(MeanNumConcepts(inc,lim),2), mean(Delta(inc,lim),2), '*r');



subplot(1,3,3)
hold on
plot(mean(Fitness_level(:,lim),2), mean(Delta(:,lim),2), '*k');
plot(mean(Fitness_level(inc,lim),2), mean(Delta(inc,lim),2), '*r');