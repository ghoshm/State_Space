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

%% RCE
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'score');
%load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'sleep_cells');

X = datasample(score(:,1:3),10000)'; % 40000 limit
%X = datasample(sleep_cells(:,3),10000)'; X = [X ; X];
N = size(X,2); % number of samples

% set the fuzzifier constant to 1.4
m = 1.4;

% Optimize the swarm using 80% resampling rate and mahalanobis distance
swarm = RCE(X, 15, 'distance','euclidean','fuzzifier',m, 'display','text', ...
    'swarm',200, 'subsprob',0.03, 'maxiter',100,'resampling_rate',0.8,'calculate_labels', false);

% input vectors using the Swarm
[softlabels, crisplabels, numlabels] = swarm_cluster(X,swarm);

% plot the fuzzy voronoi cells on the 1st and 2nd dimension
visualize_swarm(X,swarm,1,2,m,200)

% Perform fuzzy evidence accumulation on the swarm
ensemble = EnsembleAggregate(softlabels,'complete',true);

% plot the scatter matrix
figure('name','scatterplot');
gplotmatrix(X(:,:)',[],ensemble.ensemble_labels)

%% Work on Results 
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\Results.mat')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'i_experiment_tags')
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'i_group_tags')
tags =i_group_tags(i_experiment_tags == 9); 

cmap = lbmap(max(tags),'RedBlue'); 

clf; hold on; 
for t = 1:max(tags) 
    plot(diff(BIC{1,1}(tags == t,:))','color',cmap(t,:)); 
    pause(3)
end 

% Raw 
for f = 1:size(tags,1) % for each fish 
    scrap(f,1) = knee_pt(AIC{2,1}(f,:)); 
    scrap(f,2) = knee_pt(BIC{2,1}(f,:));
end 

% Smoothed 
for f = 1:size(tags,1) % for each fish 
    scrap(f,1) = knee_pt(smooth(AIC{2,1}(f,:),3)); 
    scrap(f,2) = knee_pt(smooth(BIC{2,1}(f,:),3));
end 

%% Simple Example
%https://pdfs.semanticscholar.org/04d5/4946c3438fca1ed41341f8f4d2bac1fbb00a.pdf

% settings 
sample = 1000; % number of points to sample 
its = 100; % number of iterations 
t = 0.2; % adjacency threshold 

% data 
idx = nan(sample,its,'single'); % allocate 
tic
for i = 1:sample % for each sample 
    if i <= sample/2 % for the first half of the points 
        idx(i,:) = datasample([1 2 3],its);
    else
        idx(i,:) = datasample([3 4 5],its);
    end
end
toc 

% Clustering 
swarm = RCE(idx,2) 

% Evidence Acumulation
tic
idx_ea = sparse(size(idx,1),size(idx,1)); % allocate (n x n)
for i = 1:size(idx,1) % for each point
    v = zeros(1,size(idx,1),'single');
    for ii = (i + 1):size(idx,1) % for each comparison
        v(1,ii) = sum(idx(i,:) == idx(ii,:));
    end
    v = v/its; % divide by its 
    v(v < t) = 0; % threshold
    idx_ea(i,i:end) = v(1,i:end); % fill 
    idx_ea(i:end,i) = v(1,i:end); % fill symetrically 
    
    if mod(i,100) == 0
       disp(num2str(i)); 
    end 
end
toc 

% Clustering the EA Matrix 
[C] = SpectralClustering(idx_ea, 2,1); 
clf; plot(C); 

% SVD 
% [U,S,V] = svds(idx_ea,10); 
% Y = U*(S^-1); 

% https://www.quora.com/Should-I-use-the-U-or-V-matrix-returned-by-U-S-V-svd-data-in-MATLAB-to-project-my-data-to-a-lower-dimensional-space
% https://uk.mathworks.com/matlabcentral/answers/379546-pca-and-data-projection-issue
% X = X-mean(X);          % 'center' the data around zero
% A = (X'*X) / length(X); % compute the covariance matrix (normalised by the num of elements)
% [V,~] = eig(A);         % compute the eigenvectors -- this results in a increasing ordering
% V = fliplr(V);          % flip the eigenvectors so that the most significant come first
% V = V(:,1:D);           % take only those eigenvectors required
% Y = X * V;              % project the original data to the new coordinate system

