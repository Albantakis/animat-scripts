function Animat_RandomIIT_Analysis

Ntrials = 256;
NumProb = 0;
NumEl = 2;

BigPhi = NaN(Ntrials,NumProb);
BigPhiMip = NaN(Ntrials,NumProb);
NumConcepts = NaN(Ntrials,NumProb);
SizeComplex = NaN(Ntrials,NumProb);
HOConcepts = NaN(Ntrials,NumProb);
MaxStates = 1;
count = 0;
for prob0 =0% 0.1:0.1:0.9 %[0.05 0.1:0.1:0.9 0.95]
    count = count +1;
%for r = 1:Ntrials
for r = 1:Ntrials

    cd(strcat('Random_trials', int2str(prob0*100)))
    if prob0 == 0.05 || prob0 == 0.95
        FilenameA = strcat('Random', int2str(r), '_MC', int2str(prob0*1000), '.mat');
    else
        if NumEl == 4
            FilenameA = strcat('Random', int2str(r), '_MC', int2str(prob0*10), '.mat');        
        elseif NumEl == 2
            FilenameA = strcat('Random', int2str(NumEl), '_', int2str(r-1), '_MC', int2str(prob0*1000), '.mat');
        else
            FilenameA = strcat('Random', int2str(NumEl), '_', int2str(r), '_MC', int2str(prob0*1000), '.mat');
        end
    end

    if exist(FilenameA,'file') == 2 
        load(FilenameA)
        cd ..

        Nstates = numel(results.state);
        % initialize 
        subsys_numel = NaN(MaxStates,1);
        Phi = NaN(MaxStates,1);
        Phi_MIP = NaN(MaxStates,1);
        Num_Concepts = NaN(MaxStates,1);
        MConceptOrder = NaN(MaxStates,1);
        for j = 1:Nstates
            subsys = results.Complex{j};
            sys_results = results.state(j);
            if isempty(sys_results.Phi)
                Phi(j) = NaN;
                Phi_MIP(j) = NaN;
                Num_Concepts(j) = NaN;
                subsys_numel(j) = NaN;
                MConceptOrder(j) = NaN;
            else
                Phi(j) = sys_results.Phi;
                % Larissa: This is because it seems that Complex by default is 1
                if Phi(j) > 0
                    subsys_numel(j) = numel(subsys);
                else 
                    subsys_numel(j) = 0;
                end    

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
        end

        BigPhi(r,count) = Phi; %weighted ave/expected value
        BigPhiMip(r,count) = Phi_MIP;  
        %Num_Conn(r,:) = results.numConn;
        NumConcepts(r,count) = Num_Concepts;
        SizeComplex(r,count) = subsys_numel;
        HOConcepts(r,count) = MConceptOrder;
    else 
        cd ..
    end      

end 
end

File = strcat('Random_data', int2str(NumEl));
save(File, 'BigPhi','BigPhiMip','NumConcepts','HOConcepts', 'SizeComplex');



figure(500)
subplot(3,1,1)
    hold on
    plot(1:length(BigPhiMip), BigPhiMip)
    xlim([1 100])
subplot(3,1,2)    
    hold on
    plot(1:length(BigPhiMip), NumConcepts)
    xlim([1 100])
subplot(3,1,3)
    hold on
    plot(1:length(BigPhiMip), SizeComplex)
    xlim([1 100])


