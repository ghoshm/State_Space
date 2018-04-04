% State_Space_6 

load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'score');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'fish_tags');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'sleep_cells');

% Remember to deal with NaN Values! 

% Settings
reps = 1; % set the number of repetitions
k_max = 1:20; % set values of k (clusters) to try
options = statset('MaxIter',1000); % Hard coded number of iterations
knee_dim = 3; % pca dimensions to keep  

X{1,1} = score(:,1:knee_dim); % active bouts 
X{2,1} = sleep_cells(:,3); % inactive bouts 

% Calculate Regularization
score_values = unique(X{1,1}(:)')'; % Find unique scores
score_zero = knnsearch(score_values,0); % Find the closest to zero
rv = abs(score_values(score_zero)); % Regularization value

% submit fish seperately 
X{1,1} = X{1,1}(fish_tags{1,1} == 578,:); 
X{2,1} = X{2,1}(fish_tags{2,1} == 578,:); 

% Function 

% allocate 
GMModels = cell(2,k_max(end)); % models 
idx = cell(2,k_max(end)); % cluster assignments 
P = cell(2,k_max(end)); % posterior probabilities 
BIC = zeros(2,k_max(end),'single'); % info criteria 

% Loop
for s = 1:2 % for active/inactive 
    
    tic 
    for k = k_max

        GMModels{s,k} = fitgmdist(X{s,1},k,...
            'Options',options,'RegularizationValue',...
            rv,'Replicates',5); % Fit K gaussians
        
        [idx{s,k},~,P{s,k}] = cluster(GMModels{s,k},X{s,1});
        P{s,k} = max(P{s,k},[],2); % Keep only assigned probabilities (helps memory)
        
        % Information Criteria
        BIC(s,k)= GMModels{s,k}.BIC; % Extract BIC
        
        disp(num2str(k)); 
    end
    toc
    
end

% save best structures 

