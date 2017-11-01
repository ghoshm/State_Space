% Down-Sampling Script 

%% Required Scripts 

% dir2 - marcus.ghosh.11@ucl.ac.uk 

% Nat Sort Files -
    %http://uk.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort 
    
%% Load data 

folder_path = uigetdir; % Choose a folder 
folder_open = dir2(folder_path); % Open this folder

num_clusters = size(folder_open,1); % Calculate the number of clusters 
disp(horzcat('Number of clusters = ',num2str(num_clusters))); % Feedback  

% Pre-allocate 
GMModels = cell(num_clusters,1);
AIC = zeros(1,num_clusters); 
BIC = zeros(1,num_clusters); 
idx = cell(num_clusters,1);
nlogl = cell(num_clusters,1);
P = cell(num_clusters,1);
sheet_names = cell(size(folder_open,1),1); % Excel sheet names  

% Ordering by File Name 
for f = 1:size(folder_open,1) % For each excel file 
    sheet_names{f} = folder_open(f).name; % Take it's name 
end % Note that these will be in "computer" order 
    % Ie. 1-10, 100, 1000 etc 

[~,O] = natsortfiles(sheet_names); % Use natsortfiles to sort by file name
    clear sheet_names; % Clear Sheet names 
    
% Merge data 
counter = 1; % start a counter 
for k = O' % For each cluster 
    
    load(strcat(folder_path,'\',folder_open(k).name)); % load data
    %GMModels{counter,1} = results{1}(counter,1); 
    AIC(counter) = results{2}(1,counter); 
    BIC(counter) = results{3}(1,counter); 
    
    if counter == 1
        idx{counter,1} = results{4};
        nlogl{counter,1} = results{5};
        P{counter,1} = results{6};
    else
        idx{counter,1} = results{4}{counter,1};
        nlogl{counter,1} = results{5}{counter,1};
        P{counter,1} = results{6}{counter,1};
    end
    
    counter = counter + 1; % add to counter 
    clear results 
end 

clear counter f folder_open folder_path GMModels O

%% Load in Data + Post Score Data (From PCA) 
load('F:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Legion_Version\data.mat'); 
load('F:\Behaviour\SleepWake\Re_Runs\Clustered_Data\Legion_Version\Wake_cells_2d.mat'); 

[~,numComp] = min(AIC); 
cmap_cluster = hsv(numComp); % generate a colormap for the clusters 

%% Sorting Clusters 

% PC1 - seems to use length 
mean_cluster_length = nan(1,numComp,'single'); % pre-allocate
for c = 1:numComp % For each cluster 
    mean_cluster_length(c) = max(wake_cells(idx{numComp,1}==c,3));
        % Calculate max bout length 
end 

[~,O] = sort(mean_cluster_length); % Sort by increasing bout length 

idx_numComp_sorted = nan(size(idx{numComp,1}),'single'); % Pre-allocate 

for c = 1:numComp  % For each cluster
    idx_numComp_sorted(idx{numComp,1} == O(c),:) = c; 
        % Re-assign cluster numbers 
end 

clear mean_cluster_length c O; 

%% Down-Sampling Curve 

% Settings 
ds = [0.01 0.1 1 10 25 50 75]; % Set your down-sample %'s 
ds = ds/100; 
options = statset('MaxIter',1000); % Hard coded number of iterations 
X(:,end+1) = 1:size(X,1); % Track each bout number (% Start here!) 

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
    sample_per_c = round((size(idx{numComp,1},1)*(ds(d)))/...
        (max(experiment_tags{1,1})*max(parameter_indicies{1,1})*max(group_tags{1,1}))); 
    
    % Sample X equally across conditions 
    for e = 1:max(experiment_tags{1,1}) % For each experiment 
        for t = 1:max(parameter_indicies{1,1}) % For each time point 
            for g = 1:max(group_tags{1,1}) % For each group 
                sample = [sample ; datasample(X(experiment_tags{1,1} == e &...
                    parameter_indicies{1,1} == t & group_tags{1,1} == g,:),...
                    sample_per_c,'Replace','false')]; 
            end 
        end 
    end 
    
    % GM Fitting & Clustering 
    % Set-up 
    score_values = unique(sample); % Find unique scores
    score_zero = knnsearch(score_values,0); % Find the closest to zero
    rv = score_values(score_zero); % Regularization value
    
    GMModels_ds{d} = fitgmdist(sample,numComp,...
        'Options',options,'RegularizationValue',...
        abs(rv),'Replicates',5); % Fit K gaussians
    
    AIC_ds(d)= GMModels_ds{d}.AIC; % Extract AIC
    BIC_ds(d)= GMModels_ds{d}.BIC; % Extract BIC
    
    % Cluster using this mixing
    [idx_ds{d},nlogl_ds{d},P_ds{d}] = cluster(GMModels_ds{d},sample);

    % Sorting Clusters
    mean_cluster_length = nan(1,numComp,'single'); % pre-allocate
    for c = 1:numComp % For each cluster
        mean_cluster_length(c) = max(wake_cells(idx_ds{d,1}==c,3));
        % Calculate mean X1 value
    end
    
    [~,O] = sort(mean_cluster_length); % Sort by increasing bout length
    
    idx_numComp_sorted_ds{d,1} = nan(size(idx_ds{d,1}),'single'); % Pre-allocate
    
    for c = 1:numComp  % For each cluster
        idx_numComp_sorted_ds{d,1}(idx_ds{d,1} == O(c),:) = c;
        % Re-assign cluster numbers
    end
    
    % Assign the remaining data to clusters 
    
    disp(horzcat('Finished down-sampling using ',num2str(ds(d)*100),...
        '% of the data'));
end 

%% Figure Workings 

figure; hold on; 
for c = 1:numComp  % For each cluster 
    scatter(sample(idx_numComp_sorted_ds{1,1}==c,1),sample(idx_numComp_sorted_ds{1,1}==c,2),...
        'markerfacecolor',cmap_cluster(c,:),...
        'markeredgecolor',cmap_cluster(c,:));
    xlabel('PC1','Fontsize',12); 
    ylabel('PC2','Fontsize',12); 
end 

% Scattering bouts in PCA Space  
figure; hold on; 
for c = 1:numComp  % For each cluster 
    scatter(X(idx_numComp_sorted==c,1),X(idx_numComp_sorted==c,2),...
        'markerfacecolor',cmap_cluster(c,:),...
        'markeredgecolor',cmap_cluster(c,:));
    xlabel('PC1','Fontsize',12); 
    ylabel('PC2','Fontsize',12); 
    pause(2); 
end 

% Overlay 
figure;
for c = 1:numComp % For each cluster
    subplot(3,4,c); hold on;
    scatter(X(idx_numComp_sorted==c,1),X(idx_numComp_sorted==c,2),...
        'markerfacecolor','k',...
        'markeredgecolor','k');
%     scatter(sample(idx_numComp_sorted_ds{1,1}==c,1),sample(idx_numComp_sorted_ds{1,1}==c,2),...
%         'markerfacecolor',cmap_cluster(c,:),...
%         'markeredgecolor',cmap_cluster(c,:));
    xlabel('PC1','Fontsize',12);
    ylabel('PC2','Fontsize',12);
end