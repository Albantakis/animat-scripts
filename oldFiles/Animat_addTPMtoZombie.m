% add tpm to zombie files
clear all
%% Animat specificities
numSen = 2;
numMot = 2;
numNodes = 8;
tic
% -------------------------------------------------------------------------
for TrialNum = 4
TrialType = 'c1a3_36'
TrialNum
%AnimatPath = strcat('~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_', TrialType , '/trial', int2str(TrialNum),'_');
AnimatPath = strcat('/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_', TrialType , '/trial', int2str(TrialNum),'_');

% -------------------------------------------------------------------------
%% options for one state and KLD and sum of small phis
%in_options = [3     1     2     1     1     0     0     1     1     0     0     0     1     1     1     0     0     1     0];
%all states KLD and sum of small phis
%in_options = [0     1     1     0     0     1     1     0     0     3     1     2     1     0     0     1     1     0];
%all states L1 and L1 norm for complexes
%in_options = [0     1     1     1     1     1     1     0     2     0     0     1     1     0     0     1     1     0];
%all states EMD and L1 norm for constellations only (9: parfor,10:strongconn, 11:freeze)
%Without normalization ([6,7] = 0)
in_options =  [0     1    1     2     2     0     0     0     1     0     0     1     1     0     0     1     1     0];

  
%Larissa: TODO check this and options in general!
op_average = in_options(2); % 0: use a specific current state 1: average over all possible current states
op_write = 2; % 1 --> write output to file

if op_write > 0
    Foldername = strcat('Freeze_', TrialType,'_trial', int2str(TrialNum));
    mkdir(Foldername)
end 

%% parallel computing
% if a pool is open, close it
if matlabpool('size')
    matlabpool close force;
end
% if parallel option is on, open a new pool
op_parallel = in_options(1);
if op_parallel
    matlabpool;
end

%% begin timer and disp notification
tic
fprintf('\nRunning...\n\n')

%% for loop across generations
for g = 58368%[0:512:60000-1];
    results = [];
    
    Animat_gen = g;
        
    J_tempfile = strcat(AnimatPath, int2str(Animat_gen), '_EdgeList.txt');
    J_temp = load(J_tempfile);

    if ~isempty(J_temp)
        J_temp = unique(J_temp, 'rows')+1;
        % MOTORS ARE SET TO 0 -> they don't actually have recurrent
        % connections
        J_temp = J_temp(J_temp(:,1) <= numNodes-numMot,:); 
        %Sensors might have incoming connection, but because they actually don't do anything they shouldn't be taken into account
        J_temp = J_temp(J_temp(:,2) > numSen,:);
        
        %%Larissa: If one sensor only, sensor 2 is always set to 0
%         if g < 30000
%        J_temp = J_temp(J_temp(:,1) ~= 2,:);
%         J_temp = J_temp(J_temp(:,2) ~= 8,:);
%         end
                
        J_sparse = sparse(J_temp(:,1), J_temp(:,2),1,numNodes,numNodes);
            
        J = full(J_sparse)';     
        
        %TODO CHECK IF TPM Is correct with large J
        [tpm, used_nodes] = Animat_ReducedTpmSmall(Animat_gen, AnimatPath, numSen, numMot, J, in_options(3));
        tpm(:,used_nodes < numSen+1) = 0.5; %Sensors might be switched on through mechanism but that is overwritten by environment --> doesn't do anything
