% State_Space_6 

load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'score');
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'sleep_cells');
knee_dim = 2; % pca dimensions to keep  
% Handling NaN Values 
    % Note that there are so few NaN values that giving them "fake" values 
    % for the clustering won't make a difference 
sleep_cells_nan_track = isnan(sleep_cells(:,3)); % store nan locations  
sleep_cells(sleep_cells_nan_track,3) = 1; % set NaN's to 1 (the mode ~= 18% of data) 
X{1,1} = score(:,1:knee_dim); % active bouts 
X{2,1} = sleep_cells(:,3); % inactive bouts 

% Example Data 
% X{1,1} = [normrnd(1,1,1000,2) ; normrnd(10,1,1000,2)]; 
% X{2,1} = [normrnd(1,1,1000,1) ; normrnd(10,1,1000,1)]; 

% Settings
reps = 100; % set the number of repetitions
k_vals = 2:40; % set values of k (clusters) to try
a_size = 10000; % number of points to check  
s_size = 100000; % number of points to sample 
GMM_reps = 5; % number of GMM Models to fit per iteration 
max_its = 1000; % Hard coded number of iterations
method = 'average'; % linkage measure 
nn = 9; % number of nearest neighbours

% Calculate Regularization
score_values = unique(X{1,1}(:)')'; % Find unique values
score_zero = knnsearch(score_values,0); % Find the closest to zero
rv = abs(score_values(score_zero)); % Regularization value for GMM 

% Loop
for s = 1:2 % for active & inactive
    tic 
    [ea{s,1}, idx{s,1}, idx_cts{s,1}, ~, ...
        ea_links{s,1}, ea_idx{s,1}, ~, th(s,1), sample_a{s,1}] = ...
        gmm_sample_ea(X{s,1},reps,k_vals,a_size,s_size,rv,GMM_reps,max_its,method,nn);
    toc 
end

% Cluster Colormap
numComp = [max(ea_idx{1,1}) max(ea_idx{2,1})]; % Choose active & inactive numbers of clusters
scrap = lbmap(sum(numComp),'RedBlue'); % Color Scheme
cmap_cluster{1,1} = scrap(1:numComp(1),:); % generate a colormap for the clusters 
cmap_cluster{2,1} = scrap((numComp(1)+1):end,:); % generate a colormap for the clusters 

clear s score_values score_zero scrap 

%% Clustering Figures 

% Evidence Accumulation
for s = 1:2 % for active & inactive
    figure; box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
    set(gca,'FontName','Calibri');
    
    % Dendrogram
    subplot('position',[0.0500    0.8178    0.9000    0.1322]);
    [H,~,perm] = dendrogram(ea_links{s,1},size(ea{s,1},1),'colorthreshold',th(s,1)); % dendrogram
    lineColours = cell2mat(get(H,'Color'));
    colourList = unique(lineColours,'rows');
    
    for c = 2:size(colourList,1) % for each colour
        i = ismember(lineColours, colourList(c,:),'rows');
        lineColours(i, :) = repmat(cmap_cluster{s,1}(c-1,:),sum(i),1);
    end
    for l = 1:size(H,1)
        set(H(l), 'Color', lineColours(l,:));
    end
    axis off;
    
    line(get(gca,'xlim'), [th(s,1) th(s,1)],'LineStyle',':','LineWidth',1.5);
    text(1,double(th(s,1)),'Maximum Lifetime Cut','verticalalignment','bottom',...
        'FontName','Calibri','FontSize',16);
    
    % EA Matrix
    subplot('position', [0.0500    0.0500    0.9000    (0.7294+0.0184)]);
    imagesc(ea{s,1}(perm,perm) );
    axis off
    c = colorbar;
    c.Label.String = 'E.A. Index';
    c.Location = 'southoutside';
    c.FontSize = 16;
end

clear s H perm lineColours colourList c i l 

%% Centroids Scatter 
figure; hold on; 
scatter(X{1,1}(:,1),X{1,1}(:,2),'.k'); % scatter data 
scatter(idx_cts{1,1}(:,1),idx_cts{1,1}(:,2),'+','r'); 

%% Using Posterior Probabilities 

tic
idx_cts = [];
for r = 1:reps % for each iteration
    
    clear sample k GMModels idx P;
    
    sample = datasample(X{s,1},s_size); % sample points
    k = datasample(k_vals,1); % choose a value for k
    
    GMModels = fitgmdist(sample,k,...
        'Options',options,'RegularizationValue',...
        rv,'Replicates',5); % Fit k gaussians
    
    [idx,~,P] = cluster(GMModels,X{s,1}(sample_a(:,s),:)); % cluster the answer points
    idx_cts = [idx_cts ; grpstats(X{s,1}(sample_a(:,s),:),idx,'mean')]; % answer centroids
    
    ea = ea + (P*P'); % add to the EA Matrix
    
    % Information Criteria
    disp(num2str(r));
end
toc

D = diag(diag(ea)); % diagonals
ea = D^-0.5 * ea * D^-0.5; % rescale

% construct the dendrogram
ea_dist = 1-ea; % invert
ea_dist(eye(size(ea))==1) = 0;
Z = linkage(squareform(ea_dist),'average');

%% RCE
load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'score');
%load('D:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\New\180111.mat', 'sleep_cells');

X = datasample(score(:,1:3),10000)'; % 40000 limit
%X = datasample(sleep_cells(:,3),10000)'; X = [X ; X];
N = size(X,2); % number of samples

% set the fuzzifier constant to 1.4
m = 1.4;

% Note: 180409 Start by evaluating cluster solutions 
   % https://uk.mathworks.com/help/stats/evalclusters.html 

% Optimize the swarm using 80% resampling rate and mahalanobis distance
swarm = RCE(X, 6, 'distance','euclidean','fuzzifier',m, 'display','text', ...
    'swarm',100, 'subsprob',0.03, 'maxiter',100,'resampling_rate',0.8,'calculate_labels', false);

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

