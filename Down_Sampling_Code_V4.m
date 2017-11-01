% Down-Sampling Script V4

    % V3 - Using Issac's idea of centroid consistency 
    % V4 - Samples evenly from the bout data space 
    
%% Required Scripts 

% dir2 - marcus.ghosh.11@ucl.ac.uk 

% Nat Sort Files -
    %http://uk.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort 
    
%% Load in Data

load('F:\Behaviour\SleepWake\Re_Runs\Post_State_Space_Data\WT\Test.mat')

%% Settings 

clear X; 
X = zscore(wake_cells(:,3:end));  % z-score the data 
[coeff,score,~,~,explained,~] = pca(X); % pca 
[knee_dim] = knee_pt(explained); % Choose this many dimensions 
disp(horzcat('Reduced data to ',num2str(knee_dim),' dimensions')); 
X = score(:,1:knee_dim);  

% Hard
ds = [0.01 0.1 1 5 10 15 25 50]; % Set your down-sample %'s 
%ds = [0.01 0.1 1 5 10 15]; % first pass  
%ds = [25 50]; % second pass 
reps = 3; % set your repetitions 
options = statset('MaxIter',1000); % Hard coded number of iterations 

% Soft 
ds = ds/100; 

% Pre-allocation
GMModels_ds = cell(1);
idx_ds = cell(1);
cluster_centroids = cell(size(idx_ds)); 
time_sheet = nan(size(ds,2),reps,'single');

%% Ensure Even Sampling 

% Caluclate probability density in 2d 
[N,~,~,binX,binY] = histcounts2(X(:,1),X(:,2),...
    [100 100],'Normalization','probability');

% Assign a probability to each bout 
bout_probs = nan(size(binX),'single'); % pre-allocate 
for ex = 1:max(binX) % for each x-bin
    for ey = 1:max(binY) % for each y-bin
        if N(ex,ey) > 0 % for probs > 0 
            bout_probs(binX == ex & binY == ey,1) = N(ex,ey);  
                % assign probs to bouts 
        end
    end
    disp(horzcat('Assigned Bout probabilities ',num2str(ex),...
        ' of ',num2str(max(binX)))); % report progress 
end

% Invert bout probabilities 
bout_probs_lin = 1 - bout_probs; 
    % Note that this works best of various version's I've tried 
    
%% Calculation
for d = 1:size(ds,2) % For each downsample
        
    for r = 1:reps % for each repeat
        tic
        sample = []; %sample_tags = []; 
        GMModels_ds = cell(1);
        idx_ds = cell(1);
        clear score_values score_zero rv O; 
        
        % Weighted Sample
        [sample,~] = datasample(X,round(size(X,1)*ds(d)),...
            1,'replace',false,'weights',bout_probs_lin);
            % My intuition is that sampling without replacement is the
            % correct approach, though this may be worth revisiting in
            % future (170926 - MG). 
            
        % GM Fitting & Clustering
        % Set-up
        score_values = unique(sample); % Find unique scores
        score_zero = knnsearch(score_values,0); % Find the closest to zero
        rv = score_values(score_zero); % Regularization value
        
        GMModels_ds{1} = fitgmdist(sample,numComp(1),...
            'Options',options,'RegularizationValue',...
            abs(rv),'Replicates',5); % Fit K gaussians
        
        % Cluster using this mixing
        idx_ds{1} = cluster(GMModels_ds{1},sample);
        
        % Store cluster centroids
        for c = 1:numComp(1) % for each cluster
            cluster_centroids{d,r}(c,:) = ...
                nanmean(sample(idx_ds{1} == c,:));
        end
        
        % Sort cluster centroids by x-axis 
        [~,O] = sort(cluster_centroids{d,r}(:,1)); 
        cluster_centroids{d,r} = cluster_centroids{d,r}(O,:);
        
        % Report progress 
        disp(horzcat('Finished down-sampling & clustering using ',...
            num2str(ds(d)*100),'% of the data, repetition ',...
            num2str(r)));
        time_sheet(d,r) = toc; % record time for each iteration 
    end

end

%% Re-ordering "Ground Truth" clusters  
mean_cluster_length = nan(1,numComp(1),'single'); % pre-allocate
for c = 1:numComp(1) % For each cluster 
    mean_cluster_length(c) = nanmean(X(idx_numComp_sorted{1,1} == c,1));
        % Calculate max bout length 
end 

[~,O] = sort(mean_cluster_length); % Sort by increasing bout length 

idx_numComp_re_sorted{1,1} = nan(size(Wake.idx,1),1,'single'); % Pre-allocate 

for c = 1:numComp(1)  % For each cluster
    idx_numComp_re_sorted{1,1}(idx_numComp_sorted{1,1} == O(c),1) = c; 
        % Re-assign cluster numbers 
end 

clear mean_cluster_length c O; 

%% Figures 

% DS Data (no color)
figure; 
for d = 1:size(ds,2) % for each downsampling
    subplot(2,4,d); hold on; set(gca,'Fontsize',12); 
    title(horzcat(num2str(ds(d)*100),'% SubSamples'),'Fontsize',18); 
    for r = 1:reps % for each repetition
        for c = 1:12 % for each cluster
            plot(cluster_centroids{d,r}(c,1),cluster_centroids{d,r}(c,2),...
                'o','linewidth',1.5,'MarkerSize',6,'color','k');
            if r == reps % for the last repetition
                plot(nanmean(X(idx_numComp_re_sorted{1,1}==c,1)),...
                    nanmean(X(idx_numComp_re_sorted{1,1}==c,2)),'x','linewidth',3,...
                    'MarkerSize',12,'color',cmap_cluster{1,1}(c,:));
            end
        end
    end
    xlabel('PC 1','Fontsize',16); 
    ylabel('PC 2','Fontsize',16); 
    axis([-4 8 -3 4]);
    
end

%% Time Comparison 
figure; hold on; title('Time Taken for each Subsample','Fontsize',18)
spread_cols = plotSpread(time_sheet'/60,'showMM',4); 
spread_cols{2}(1).LineWidth = 3; % Change marker width 
spread_cols{2}(2).LineWidth = 3; % Change marker width 
set(gca,'Fontsize',12);
xticklabels({num2str(ds'*100)});
xlabel('Percentage of Data Sampled','Fontsize',16); 
ylabel('Time Taken (minutes)','Fontsize',16); 

%% Comparison 
for c = 1:numComp(1) % for each cluster 
    cluster_centroids_gt(c,:) = nanmean(X(idx_numComp_re_sorted{1,1}==c,:)); 
end 

% Compare "paired" distance to a ground truth cluster
s_dists = cell(size(ds,2),1); 
for d = 1:size(ds,2) % for each downsampling 
    for r = 1:reps % for each repetition
        s_dists{d,1}(:,r) = diag(pdist2(cluster_centroids_gt,cluster_centroids{d,r}));  
    end 
end 

for d = 1:size(ds,2) % for each downsampling 
    plot(nanmean(s_dists{d,1},2),'linewidth',3); 
    pause(3); 
end 
