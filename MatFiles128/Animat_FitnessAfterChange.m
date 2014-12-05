% plot Fitness against each other 
plotflag =3;

condi = 'c1a2'
condA = strcat(condi, '_change_c14a23');
condB = strcat(condi, '_change_c36a45');

load(strcat(condA, '_dataCB'));
range1 = range./10000;
range1S = 1:length(range);
%Zombie = load(strcat(condi, '_36_ZombiedataAllC'));
PhiMip1 = BigPhiMip;
Fitness1 = Fitness_level;
clear Fitness_level

load(strcat(condB, '_dataCB'));
PhiMip2 = BigPhiMip;
Fitness2 = Fitness_level;

if plotflag == 1
    % plot Fitness against each other
    figure
    hold on
    indP = [];
    indP0 = [];
    for i = 1:size(Fitness1,1)
        if max(BigPhiMip(i,55:59)) > 0
            indP = [indP i];
            plot(Fitness1(i,60:end)', Fitness2(i,60:end)', 'k')
            plot(Fitness1(i,60)', Fitness2(i,60)', 'ok')
            plot(Fitness1(i,end)', Fitness2(i,end)', '*k')
        else
            indP0 = [indP0 i];
            plot(Fitness1(i,60:end)', Fitness2(i,60:end)', 'b')
            plot(Fitness1(i,60)', Fitness2(i,60)', 'ob')
            plot(Fitness1(i,end)', Fitness2(i,end)', '*b')
        end
    end
    plot(mean(Fitness1(indP,60:end)), mean(Fitness2(indP,60:end)),'r')
    plot(mean(Fitness1(indP0,60:end)), mean(Fitness2(indP0,60:end)),'g')

elseif plotflag == 2
%% plot PhiMip against each other
    figure
    hold on
    indP = [];
    indP0 = [];
    for i = 1:size(PhiMip1,1)
        if max(BigPhiMip(i,55:59)) > 0
            indP = [indP i];
            plot(PhiMip1(i,60:end)', PhiMip2(i,60:end)', 'k')
            plot(PhiMip1(i,60)', PhiMip2(i,60)', 'ok')
            plot(PhiMip1(i,end)', PhiMip2(i,end)', '*k')
        else
            indP0 = [indP0 i];
            plot(PhiMip1(i,60:end)', PhiMip2(i,60:end)', 'b')
            plot(PhiMip1(i,60)', PhiMip2(i,60)', 'ob')
            plot(PhiMip1(i,end)', PhiMip2(i,end)', '*b')
        end
    end
    plot(mean(PhiMip1(indP,60:end)), mean(PhiMip2(indP,60:end)),'r')
    plot(mean(PhiMip1(indP0,60:end)), mean(PhiMip2(indP0,60:end)),'g')

elseif plotflag == 3
    color = [max(PhiMip1(:,55:59),[],2), 0.5.*ones(size(Fitness1, 1),1), max(max(PhiMip1(:,55:59),[],2)) - max(PhiMip1(:,55:59),[],2)]./max(max(PhiMip1(:,55:59),[],2));

    figure
    subplot(1,2,1)
    hold on
    for i = 1:size(Fitness1,1)
        if max(PhiMip1(i,55:59)) > 0
            plot(range, Fitness1(i,:), 'color', [color(i,:)])
        else
            plot(range, Fitness1(i,:), '-k')
        end
    end
    
    subplot(1,2,2)
    hold on
    for i = 1:size(Fitness2,1)
        if max(PhiMip1(i,55:59)) > 0
            plot(range, Fitness2(i,:), 'color', [color(i,:)])
        else
            plot(range, Fitness1(i,:), '-k')
        end
    end
end