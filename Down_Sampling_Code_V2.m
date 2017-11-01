% Down-Sampling Script V2

%% Required Scripts 

% dir2 - marcus.ghosh.11@ucl.ac.uk 

% Nat Sort Files -
    %http://uk.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort 
    
%% Load in Data

load('F:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\Test.mat')

%% Down-Sampling Curve 

% Settings 
clear X; 
X = zscore(wake_cells(:,3:end));  % z-score the data 
[coeff,score,~,~,explained,~] = pca(X); % pca 
[knee_dim] = knee_pt(explained); % Choose this many dimensions 
disp(horzcat('Reduced data to ',num2str(knee_dim),' dimensions')); 
X = score(:,1:knee_dim);  

% Hard
ds = [0.01 0.1 1 10 25 50 75]; % Set your down-sample %'s 
options = statset('MaxIter',1000); % Hard coded number of iterations 

% Soft 
ds = ds/100; 
X(:,end+1) = 1:size(X,1); % Track each bout number

% Pre-allocation
AIC_ds = zeros(1,size(ds,2)); 
BIC_ds = zeros(1,size(ds,2)); 
GMModels_ds = cell(size(ds,2),1);
idx_ds = cell(size(ds,2),1);
nlogl_ds = cell(size(ds,2),1);
P_ds = cell(size(ds,2),1);
idx_numComp_sorted_ds = cell(size(ds,2),1); 

for d = 1:size(ds,2) % For each downsample 
    sample = []; clear score_values score_zero rv mean_cluster_length c O; 
  
    % Calculate the number of samples to take per condition: 
        % experiments, days/nights, groups   
    sample_per_c = round((size(X,1)*(ds(d)))/...
        (max(experiment_tags{1,1})*max(parameter_indicies{1,1})*...
            max(group_tags{1,1}))); 
    
    % Sample X equally across conditions 
    for e = 1:max(experiment_tags{1,1}) % For each experiment 
        for t = 1:max(parameter_indicies{1,1}) % For each time point 
            for g = 1:max(group_tags{1,1}) % For each group 
                sample = [sample ; datasample(X(experiment_tags{1,1} == e &...
                    parameter_indicies{1,1} == t & group_tags{1,1} == g,:),...
                        sample_per_c,'Replace',false)]; 
            end 
        end 
    end 
    
    % Keep track of which bout's you've sampled 
    sample_tags = sample(:,end); 
    sample = sample(:,1:end-1); 
    
    % GM Fitting & Clustering 
    % Set-up 
    score_values = unique(sample); % Find unique scores
    score_zero = knnsearch(score_values,0); % Find the closest to zero
    rv = score_values(score_zero); % Regularization value
    
    GMModels_ds{d} = fitgmdist(sample,numComp(1),...
        'Options',options,'RegularizationValue',...
        abs(rv),'Replicates',5); % Fit K gaussians
    
    AIC_ds(d)= GMModels_ds{d}.AIC; % Extract AIC
    BIC_ds(d)= GMModels_ds{d}.BIC; % Extract BIC
    
    % Cluster using this mixing
    [idx_ds{d},nlogl_ds{d},P_ds{d}] = cluster(GMModels_ds{d},sample);

    % Sorting Clusters
    mean_cluster_length = nan(1,numComp(1),'single'); % pre-allocate
    for c = 1:numComp(1) % For each cluster
        mean_cluster_length(c) = nanmean(idx_numComp_sorted{1,1}...
            (sample_tags(idx_ds{d,1}==c),1)); 
        % Calculate mean X1 value
    end
    
    [~,O] = sort(mean_cluster_length); % Sort by increasing bout length
    
    idx_numComp_sorted_ds{d,1} = nan(size(idx_ds{d,1}),'single'); % Pre-allocate
    
    for c = 1:numComp(1)  % For each cluster
        idx_numComp_sorted_ds{d,1}(idx_ds{d,1} == O(c),:) = c;
        % Re-assign cluster numbers
    end
    
    % Assign the remaining data to clusters 
    
    disp(horzcat('Finished down-sampling & analysis using ',...
        num2str(ds(d)*100),'% of the data'));
end 

%% Figure Workings 

figure; hold on; 
for c = 1:numComp(1)  % For each cluster 
    scatter(sample(idx_numComp_sorted_ds{1,1}==c,1),sample(idx_numComp_sorted_ds{1,1}==c,2),...
        'markerfacecolor',cmap_cluster{1,1}(c,:),...
        'markeredgecolor',cmap_cluster{1,1}(c,:));
    xlabel('PC1','Fontsize',12); 
    ylabel('PC2','Fontsize',12); 
    pause(3)
end 

% Overlay 
figure; hold on
for c = 1:numComp(1) % For each cluster
    subplot(3,4,c); hold on;
      scatter(sample(idx_numComp_sorted_ds{1,1}==c,1),sample(idx_numComp_sorted_ds{1,1}==c,2),...
          'markerfacecolor',cmap_cluster{1,1}(c,:),...
         'markeredgecolor',cmap_cluster{1,1}(c,:));
         scatter(X(idx_numComp_sorted{1,1}==c,1),X(idx_numComp_sorted{1,1}==c,2),...
        'markerfacecolor','k',...
        'markeredgecolor','k');
    xlabel('PC1','Fontsize',12);
    ylabel('PC2','Fontsize',12);
end

%% Sorting Clusters
[a,b] = min(min(pdist2(scrap_centroids,scrap_centroids_two))); 


    
%% Finding a feature to sort the clusters by

for c = 1:numComp(1) % for each cluster 
    scrap_centroids(c,:) = nanmean(X(idx_numComp_sorted{1,1} == c,1:2)); 
    scrap_centroids_two(c,:) = nanmean(X(sample_tags(idx_numComp_sorted_ds{d,1}==c),1:2)); 
end 

figure; hold on; axis([-5 5 -5 5])
for c = 1:numComp(1) 
    text(scrap_centoids(c,1),scrap_centroids(c,2),{num2str(c)})
    text(scrap_centroids_two(c,1),scrap_centroids_two(c,2),{num2str(c)},...
        'Color','r')
end 

figure; hold on; 
for c = 1:numComp(1)
    scatter(X(idx_numComp_sorted{1,1}==c,1),X(idx_numComp_sorted{1,1}==c,2),...
        'markerfacecolor',cmap_cluster{1,1}(c,:),...
        'markeredgecolor',cmap_cluster{1,1}(c,:));
end 
for c = 1:numComp(1)
    text(scrap(c,1),scrap(c,2),{num2str(c)},'fontsize',20,'color','k'); 
end 

for c = 1:numComp(1) 
    cluster_boundaries{c} = boundary(X(idx_numComp_sorted{1,1}==c,1),...
        X(idx_numComp_sorted{1,1}==c,2),1); 
end 

figure; hold on; 
for c = 1:numComp(1)
    clear b_x b_y; 
%     scatter(X(idx_numComp_sorted{1,1}==c,1),X(idx_numComp_sorted{1,1}==c,2),...
%         'markerfacecolor',cmap_cluster{1,1}(c,:),...
%         'markeredgecolor',cmap_cluster{1,1}(c,:));
    b_x = X(idx_numComp_sorted{1,1}==c,1); 
    b_y = X(idx_numComp_sorted{1,1}==c,2); 
    
    plot(b_x(cluster_boundaries{c}),b_y(cluster_boundaries{c}),...
        'linewidth',3,'color','k'); 
end 

