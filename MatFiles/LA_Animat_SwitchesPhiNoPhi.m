Switches = BigPhiMip; 
Switches(Fitness_level ~= 128) = 0;
Switches(Switches > 0) = 1;
DiffSwitches = diff(Switches,1,2);
DiffAllSwitchesPhi = diff(BigPhiMip > 0, 1, 2);
PhiSwitches = abs(DiffAllSwitchesPhi).*DiffSwitches;

figure; imagesc(PhiSwitches)

nnz(PhiSwitches)
numel(find(PhiSwitches == 1)) %from 0 to > 0
numel(find(PhiSwitches == -1)) %from > 0 to 0
