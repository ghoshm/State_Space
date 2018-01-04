% State_Space_5 

% V4 - Aims to use previous state space versions as a rough frame work to now 
    % work with data from Bed_frames (frame by frame data)
% V5 - Fit's a GMM to a sub-sample of the data, then assigns all data 
    % to cluster's based on posterior probabilities 
    
%% Required scripts 

% Lbmap 
%https://uk.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps

% Knee_pt 
%https://uk.mathworks.com/matlabcentral/fileexchange/35094-knee-point
%% Load bout structure data (Parameter Extracted Data)
tic
% Load Data - using multiselect
[filename, pathname] = uigetfile('*.mat', 'Select files','MultiSelect','on'); %Select files
if isequal(filename,0) %If no file is selected
    error('No Files Selected') %Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) %Show selected filenames
end

% Data Structures (Concatenate Variables)
% Note for (2,1) cells - 1 = wake, 2 = sleep
wake_cells = [];
sleep_cells = [];
i_group_tags = [];
i_experiment_tags = [];
experiment_tags = cell(2,1);
group_tags = cell(2,1);
fish_tags = cell(2,1);
parameter_indicies = cell(2,1);

counter = 1; % Start a counter
for f = 1:size(filename,2) %For each file
    clear experiment;
    experiment = load(strcat(pathname,filename{f})); %Load the mat file
    
    % Nab variables
    days_crop{f} = experiment.days_crop; % days crop
    nights_crop{f} = experiment.nights_crop; % nights crop
    parameters{f} = experiment.parameters; % parameters
    cmap{f} = experiment.cmap; % color map
    cmap_2{f} = experiment.cmap_2; % expanded color map
    night_color{f} = experiment.night_color; % night color
    geno_list{f} = experiment.geno_list; % group names
    units{f} = experiment.units; % units
    unit_conversion{f} = experiment.unit_conversion; % unit conversion
    days{f} = experiment.days; % days
    nights{f} = experiment.nights; % nights
    first_night{f} = experiment.first_night; % first night
    time_window{f} = experiment.time_window; % time window 
    fps{f} = experiment.fps; % fps 
    lb{f} = experiment.lb; % note assumes the same number of light boundaries for each experiment
    lb_sec{f} = experiment.lb_sec; % note assumes the same number of light boundaries for each experiment
    
    % Concatenate variables
    for i = 1:size(experiment.wake_cells,2) % For each fish
        wake_cells = [wake_cells ; experiment.wake_cells{1,i}]; % wake cells
        sleep_cells = [sleep_cells ; experiment.sleep_cells{1,i}]; % sleep cells
        parameter_indicies{1,1} = [parameter_indicies{1,1} ; experiment.parameter_indicies{1,i}]; % wake bout windows
        parameter_indicies{2,1} = [parameter_indicies{2,1} ; experiment.parameter_indicies{2,i}]; % sleep bout windows
        group_tags{1,1} = [group_tags{1,1} ; ones(size(experiment.wake_cells{1,i},1),1)*...
            experiment.group_tags(i,1)]; % wake group tags
        group_tags{2,1} = [group_tags{2,1} ; ones(size(experiment.sleep_cells{1,i},1),1)*...
            experiment.group_tags(i,1)]; % sleep group tags
        experiment_tags{1,1} = [experiment_tags{1,1} ; ones(size(experiment.wake_cells{1,i},1),1)*f]; % wake experiment tags
        experiment_tags{2,1} = [experiment_tags{2,1} ; ones(size(experiment.sleep_cells{1,i},1),1)*f]; % sleep experiment tags
        fish_tags{1,1} = [fish_tags{1,1} ; ones(size(experiment.wake_cells{1,i},1),1)*counter]; % wake fish tags
        fish_tags{2,1} = [fish_tags{2,1} ; ones(size(experiment.sleep_cells{1,i},1),1)*counter]; % sleep fish tags
        counter = counter + 1; % add to fish counter
    end
    
    %delta_px_sq{1,f} = experiment.delta_px_sq;
    i_group_tags = [i_group_tags ; experiment.group_tags]; % individual fish group tags
    i_experiment_tags = [i_experiment_tags ; ones(size(experiment.wake_cells,2),1)*f]; % individual fish experiment tags
    
end

clear experiment f i counter;
toc
%% GMM Settings

% Wake
tic
X{1,1} = []; % empty X
% Z-score each fishes data
for f = 1:max(fish_tags{1,1}) % for each fish
    X{1,1} = [X{1,1} ; zscore(wake_cells(fish_tags{1,1} == f,3:end))];
    if mod(f,100) == 0 % report every 100 fish 
        disp(horzcat('Completed ',num2str(f),' fish of ',...
            num2str(max(fish_tags{1,1}))));
    end
end
toc

[coeff,score,~,~,explained,~] = pca(X{1,1}); % pca
% Note: By default Matlab Centers data for PCA by subtracting
% the mean from each column (as the means are not quite zero this is
% appropriate)
[knee_dim] = knee_pt(explained); % Choose this many dimensions
disp(horzcat('Reduced Wake data to ',num2str(knee_dim),' dimensions',...
    ' Explains ',num2str(sum(explained(1:knee_dim))),' % '));
wake_cells_norm = X{1,1}; % store a z-scored version of wake_cells 
X{1,1} = score(:,1:knee_dim);

% Sleep
% Handling NaN Values 
    % Note that zscore cannot handle NaN values 
    % Note that there are so few NaN values that giving them "fake" values 
    % for the clustering won't make a difference 
sleep_cells_nan_track = isnan(sleep_cells(:,3)); % store nan locations  
sleep_cells(sleep_cells_nan_track,3) = 1; % set NaN's to 1 (the mode ~= 18% of data) 

tic
X{2,1} = []; % empty X
% Z-score each fishes data
for f = 1:max(fish_tags{2,1}) % for each fish
    X{2,1} = [X{2,1} ; zscore(sleep_cells(fish_tags{2,1} == f,3))];
    if mod(f,100) == 0 % report every 100 fish 
        disp(horzcat('Completed ',num2str(f),' fish of ',...
            num2str(max(fish_tags{2,1}))));
    end
end
toc

% Hard Coded Parameters
%ds = [0.01 0.1 1 5 10 15 25 50]; % set the number of data points to use for each clustering iteration
ds = 2000000; % two million
reps = 1; % set the number of repetitions
k_max = 1:20; % set values of k to try
options = statset('MaxIter',1000); % Hard coded number of iterations

% Save path for data
save_pathname = 'D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\New';

% Pre-allocation - {state} clusters x rep x downsample {s}(k,r,d)
for s = 1:2 % for sleep/wake
    AIC{s,1} = zeros(size(k_max,2),reps,size(ds,2));
    BIC{s,1} = zeros(size(k_max,2),reps,size(ds,2));
    cluster_centroids{s,1} = cell(size(k_max,2),reps,size(ds,2));
    time_sheet{s,1} = nan(size(k_max,2),reps,size(ds,2),'single');
end

clear f s
%% Ensure Even Sampling - Wake 

% % Caluclate probability density in 2d 
% 
% [N,~,~,binX,binY] = histcounts2(X{1,1}(:,1),X{1,1}(:,2),...
%     'Normalization','pdf');
% tic
% % Assign a probability to each bout 
% bout_probs = nan(size(binX),'single'); % pre-allocate 
% for ex = 1:max(binX) % for each x-bin
%     for ey = 1:max(binY) % for each y-bin
%         if N(ex,ey) > 0 % for probs > 0 
%             bout_probs(binX == ex & binY == ey,1) = N(ex,ey);  
%                 % assign probs to bouts 
%         end
%     end
%     if mod(ex,100) == 0
%     disp(horzcat('Assigned Bout probabilities ',num2str(ex),...
%         ' of ',num2str(max(binX)))); % report progress 
%     end 
% end
% toc
% % Invert bout probabilities 
% bout_probs_lin{1,1} = (max(bout_probs) + 1) - bout_probs; 
% % Note that even the minimum weight (1) will be chosen 
% %bout_probs_lin{1,1} = 1 - bout_probs; 
%     
% clear N binX binY ex ey bout_probs 

 %% Ensure Even Sampling - Sleep 
% 
% [N,~,binX] = histcounts(X{2,1},'Normalization','pdf');
% 
% % Assign a probability to each bout
% tic
% bout_probs = nan(size(binX),'single'); % pre-allocate 
% for ex = 1:max(binX) % for each x-bin
%     if N(ex) > 0 % for probs > 0
%         bout_probs(binX == ex,1) = N(ex);
%         % assign probs to bouts
%     end
%     if mod(ex,1000) == 0 
%     disp(horzcat('Assigned Bout probabilities ',num2str(ex),...
%         ' of ',num2str(max(binX)))); % report progress
%     end 
% end
% toc 
% 
% % Invert bout probabilities 
% bout_probs_lin{2,1} = (max(bout_probs) + 1) - bout_probs; 
% %bout_probs_lin{2,1} = 1 - bout_probs; 
%      
% clear N binX bout_probs ex 

%% Calculation

for s = 1:2 % for wake/sleep
    for d = 1:size(ds,2) % for each down-sample
        for r = 1:reps % for each repeat
            
            sample = []; sample_tags = []; 
            clear score_values score_zero rv;
            
            % Un Weighted Sample 
            [sample,sample_tags] = datasample(X{s,1},ds(d),...
                1,'replace',false);
            
            % Weighted Sample
%             [sample,~] = datasample(X{s,1},ds(d),...
%                 1,'replace',false,'weights',bout_probs_lin{s,1});
            
            % Calculate Regularization
            score_values = unique(sample); % Find unique scores
            score_zero = knnsearch(score_values,0); % Find the closest to zero
            rv = score_values(score_zero); % Regularization value
            
            % GMM
            for k = k_max % for each value of k 
                % Re-set variables 
                GMModels = cell(1);
                idx = cell(1);
                P = cell(1);
                clear O;
                
                tic % start timer
                GMModels{1} = fitgmdist(sample,k,...
                    'Options',options,'RegularizationValue',...
                    abs(rv),'Replicates',5); % Fit K gaussians
                % To the down-sampled data
                
                % Cluster the full data-set
                % Using this mixture
                [idx{1},~,P{1}] = cluster(GMModels{1},X{s,1});
                P{1} = max(P{1},[],2); % Keep only assigned probabilities (helps memory)
                
                time_sheet{s,1}(k,r,d) = toc; % record time for each iteration
                
                % Information Criteria
                AIC{s,1}(k,r,d)= GMModels{1}.AIC; % Extract AIC
                BIC{s,1}(k,r,d)= GMModels{1}.BIC; % Extract BIC
                
                % Store cluster centroids
                cluster_centroids{s,1}{k,r,d} = GMModels{1}.mu;                
                
                % Sort cluster centroids by ascending PC1 
                [~,O] = sort(cluster_centroids{s,1}{k,r,d}(:,1)); 
                cluster_centroids{s,1}{k,r,d} = cluster_centroids{s,1}{k,r,d}(O,:);
                
                % Save Data
                save(strcat(save_pathname,'\',num2str(s),'s_',num2str(ds(d)),'d_',num2str(r),'r_',num2str(k),'k','.mat'),...
                    'idx','P','GMModels','time_sheet','AIC','BIC','cluster_centroids','sample','sample_tags');
                
                % Report progress
                disp(horzcat('Finished clustering (State = ',num2str(s),') using ',...
                    num2str(ds(d)),' Points. ',' Repetition = ',num2str(r),...
                    ', k = ',num2str(k),', Time Taken = ',...
                    num2str(((time_sheet{s,1}(k,r,d)/60)/60)),' hours'));
            end
        end
    end
end

clear sample sample_tags score_values score_zero rv O ...
    GMModels idx P s d r k c;

% Save Data
save(strcat(save_pathname,'\Full_Set.mat'),'-v7.3');

%% Remember to Handle NaN Values in Sleep Cells! 

%% Choice of Sample Size 
    % Re-worked from Down_Sampling_Code_V4
    
cmap_cluster{1,1} = lbmap(max(k_max),'RedBlue'); % generate a colormap for the clusters 

% DS Data (k,r,d)
figure; 
for d = 1:size(ds,2) % for each downsampling
    subplot(2,4,d); hold on; set(gca,'Fontsize',12); 
    title(horzcat(num2str(ds(d)*100),'% SubSamples'),'Fontsize',18); 
    
    for r = 1:reps % for each repetition
        scatter(cluster_centroids{1,r,d}(:,1),cluster_centroids{1,r,d}(:,2),72,...
            cmap_cluster{1,1},'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
    end
    
    xlabel('PC 1','Fontsize',16); 
    ylabel('PC 2','Fontsize',16); 
    axis([-3 9 -3 9]);
    
end

% Comparison 
cluster_centroids_comp = nan(size(ds,2),max(k_max),'single'); % ds x k 
for d = 1:size(ds,2) % for each downsampling
    for c = 1:max(k_max) % for each cluster
        scrap = []; 
        for r = 1:reps % for each repetition
            scrap = [scrap ; cluster_centroids{1,r,d}(c,:)]; 
        end
        cluster_centroids_comp(d,c) = sum(pdist(scrap))/2;
    end
end

figure; subplot(2,1,1); set(gca,'Fontsize',12); hold on; 
title('Variation in Centroid Position','Fontsize',18); 
plot(1:size(ds,2),sum(cluster_centroids_comp,2),'color',cmap(1,:),'linewidth',3)
ylabel('Total Distance (Euclidean)','Fontsize',16); 
xticklabels({num2str(ds'*100)});
xlabel('Percentage of Data Sampled','Fontsize',16); 

% Time Figure 
ax = subplot(2,1,2); hold on; ax.FontSize = 12; 
title('Time Taken for each Subsample','Fontsize',18);
spread_cols = plotSpread(squeeze(time_sheet)/60,'showMM',4,'distributionColors',cmap_2(2,:)); 
spread_cols{2}(1).LineWidth = 3; spread_cols{2}(2).LineWidth = 3; % Change marker width 
spread_cols{2}(1).Color = cmap_2(1,:); spread_cols{2}(2).Color = cmap_2(1,:); 
xticklabels({num2str(ds'*100)});
xlabel('Percentage of Data Sampled','Fontsize',16); 
ylabel('Time Taken (minutes)','Fontsize',16); 

%% Choice of Number of Clusters (Can Alternatively Load Data from Legion Here)  

% Load data
Wake = load('F:\Behaviour\SleepWake\Re_Runs\Clustered_Data\WT\Wake_Data\Auto\Wake_Data.mat');
pause(30); % Helps with memory 
Sleep = load('F:\Behaviour\SleepWake\Re_Runs\Clustered_Data\WT\Sleep_Data\Auto\Sleep_Data.mat');
pause(30); % Helps with memory 
load('F:\Behaviour\SleepWake\Re_Runs\Clustered_Data\WT\Wake_Data\Auto\Wake_Data.mat');

% Info Criteria Figure 
figure; 
ax(1) = subplot(1,2,1); hold on; title('Active');
plot(Wake.AIC,'linewidth',3); 
plot(Wake.BIC,'linewidth',3); 

box off; set(gca, 'Layer','top'); set(gca,'Fontsize',24);
[~,icons,plots,~] = legend('AIC','BIC','Location','northeast');
legend('boxoff'); 
set(icons(1:2),'Fontsize',20) ; set(plots,'LineWidth',3);
xlabel('Number of Clusters','Fontsize',20);
ylabel('Criteria Value','Fontsize',20);

ax(2) = subplot(1,2,2); hold on; title('Inactive'); 
plot(Sleep.AIC,'linewidth',3); 
plot(Sleep.BIC,'linewidth',3);

box off; set(gca, 'Layer','top'); set(gca,'Fontsize',24);
[~,icons,plots,~] = legend('AIC','BIC','Location','northeast');
legend('boxoff'); 
set(icons(1:2),'Fontsize',20) ; set(plots,'LineWidth',3);
xlabel('Number of Clusters','Fontsize',20);
ylabel('Criteria Value','Fontsize',20);

numComp = [14 6]; % Choose wake and sleep numbers of clusters
cmap_cluster{1,1} = haxby(numComp(1)); % generate a colormap for the clusters 
cmap_cluster{2,1} = cool(numComp(2)); % generate a colormap for the clusters 

% Add to Figure 
subplot(ax(1))
scatter(numComp(1),Wake.AIC(numComp(1)),72,'k','filled'); 
subplot(ax(2)); 
scatter(numComp(2),Sleep.AIC(numComp(2)),72,'k','filled'); 

%% Bout Data Figures 
    % Should sort by something (e.g. length)
fps = 25; 

figure; 
imagesc(zscore(wake_cells(:,3:end)),[-3 3])
colormap autumn; 
c = colorbar; c.Label.String = 'Z-Score'; 
set(gca,'XTickLabels',parameters(1:6))
set(gca,'Fontsize',18)
xlabel('Parameters','Fontsize',32) 
ylabel('Bouts','Fontsize',32) 
title('Active Bout Data','Fontsize',32)

figure; 
imagesc(sleep_cells(:,3)/fps,[0 3]); 
colormap cool; 
c = colorbar; c.Label.String = 'Time (seconds)';
set(gca,'XTick',1); 
set(gca,'XTickLabels',parameters{10})
set(gca,'Fontsize',18)
xlabel('Parameters','Fontsize',32) 
ylabel('Bouts','Fontsize',32) 
title('Inactive Bout Data','Fontsize',32)

clear c fps 

%% Posterior Probabilities 
p_assign_dist{1,1} = nan(size(Wake.idx{numComp(1),1},1),1,'single'); % pre-allocate 
p_assign_dist{2,1} = nan(size(Sleep.idx{numComp(1),1},1),1,'single'); % pre-allocate 

% Wake 
for o = 1:size(p_assign_dist{1,1},1) % For each observation
    p_assign_dist{1,1}(o,1) = single(Wake.P{numComp(1),1}...
        (o,Wake.idx{numComp(1),1}(o,1))); % Pull posterior probabilities
    % For each cluster assignment
end

% Sleep 
for o = 1:size(p_assign_dist{2,1},1) % For each observation
    p_assign_dist{2,1}(o,1) = single(Sleep.P{numComp(2),1}...
        (o,Sleep.idx{numComp(2),1}(o,1))); % Pull posterior probabilities
    % For each cluster assignment
end

% Posterior Probabilities Figure 
pd = fitdist(p_assign_dist{1,1},'Kernel','Width',0.1); % Fit distribution 
p_assign_dist{1,2} = pdf(pd,0:0.1:1); clear pd; % Save as pdf  
pd = fitdist(p_assign_dist{2,1},'Kernel','Width',0.1); % Fit distribution
p_assign_dist{2,2} = pdf(pd,0:0.1:1); % Save as pdf 

figure; hold on; title('Posterior Probabilities'); 
plot(0:0.1:1,p_assign_dist{1,2},'linewidth',3)
plot(0:0.1:1,p_assign_dist{2,2},'linewidth',3)

box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32);
[~,icons,plots,~] = legend('Active Bouts','Inactive Bouts','Location','northwest');
legend('boxoff'); 
set(icons(1:2),'Fontsize',28) ; set(plots,'LineWidth',3);
xlabel('Posterior Probability','Fontsize',28);
ylabel('Probability Density','Fontsize',28);
clear ax s o pd icons plots; 

% Clear uncessary data 
Wake.idx = Wake.idx{numComp(1),1}; 
Wake.nlogl = Wake.nlogl{numComp(1),1}; 
Wake.P = Wake.P{numComp(1),1}; 

Sleep.idx = Sleep.idx{numComp(2),1}; 
Sleep.nlogl = Sleep.nlogl{numComp(2),1}; 
Sleep.P = Sleep.P{numComp(2),1}; 


%% Sorting Clusters 
    % Note that if clustering using State_Space_5 - the clusters are
    % already ordered by PC1 - so no need to re-order 
    
% Wake 
% PC1 - seems to use length 
mean_cluster_length = nan(1,numComp(1),'single'); % pre-allocate
for c = 1:numComp(1) % For each cluster 
    mean_cluster_length(c) = max(wake_cells(Wake.idx==c,3));
        % Calculate max bout length 
end 

[~,O] = sort(mean_cluster_length); % Sort by increasing bout length 

idx_numComp_sorted{1,1} = nan(size(Wake.idx,1),1,'single'); % Pre-allocate 

for c = 1:numComp(1)  % For each cluster
    idx_numComp_sorted{1,1}(Wake.idx == O(c),:) = c; 
        % Re-assign cluster numbers 
end 

clear mean_cluster_length c O; 

% Sleep 
mean_cluster_length = nan(1,numComp(2),'single'); % pre-allocate
for c = 1:numComp(2) % For each cluster 
    mean_cluster_length(c) = max(sleep_cells(Sleep.idx==c,3));
        % Calculate max bout length 
end 

[~,O] = sort(mean_cluster_length); % Sort by increasing bout length 

idx_numComp_sorted{2,1} = nan(size(Sleep.idx,1),1,'single'); % Pre-allocate 

for c = 1:numComp(2)  % For each cluster
    idx_numComp_sorted{2,1}(Sleep.idx == O(c),:) = c; 
        % Re-assign cluster numbers 
end 

clear mean_cluster_length c O; 

%% Mean "Goodness" of fit for each cluster
p_assign_dist_mean{1,1} = nan(max(fish_tags{1,1}),numComp(1),'single'); % pre-allocate  
p_assign_dist_mean{2,1} = nan(max(fish_tags{2,1}),numComp(2),'single'); % pre-allocate  
    % fish x clusters 

for f = 1:max(fish_tags{1,1}) % For each fish  
    % Wake 
    for c = 1:numComp(1) % For each active cluster 
        p_assign_dist_mean{1,1}(f,c) = nanmean(p_assign_dist{1,1}...
            (idx_numComp_sorted{1,1} == c & fish_tags{1,1} == f));  
    end 
    % Sleep 
    for c = 1:numComp(2) % For each inactive cluster 
        p_assign_dist_mean{2,1}(f,c) = nanmean(p_assign_dist{2,1}...
            (idx_numComp_sorted{2,1} == c & fish_tags{2,1} == f));  
    end 
end 

% Active Figure 
figure; hold on; % Figure
for c = 1:numComp(1) % For each active cluster
    subplot(3,5,c); title(horzcat('Active Cluster ',num2str(c)));
    for e = 1:max(experiment_tags{1,1}) % For each experiment
        spread_cols = plotSpread(p_assign_dist_mean{1,1}(i_experiment_tags == e,c),...
            'distributionIdx',i_group_tags(i_experiment_tags == e),...
            'distributionColors',cmap+(1-cmap)*(1-(1/e^.5)),'showMM',2);
        spread_cols{2}.LineWidth = 3; % Change marker width
    end
    ylabel('Posterior Probability','Fontsize',12);
    set(gca,'xticklabel',geno_list.colheaders,'Fontsize',12); % Name each group
end

% Inactive Figure 
figure; hold on; % Figure
for c = 1:numComp(2) % For each inactive cluster
    subplot(2,3,c); title(horzcat('Inactive Cluster ',num2str(c)));
    for e = 1:max(experiment_tags{2,1}) % For each experiment
        spread_cols = plotSpread(p_assign_dist_mean{2,1}(i_experiment_tags == e,c),...
            'distributionIdx',i_group_tags(i_experiment_tags == e),...
            'distributionColors',cmap+(1-cmap)*(1-(1/e^.5)),'showMM',2);
        spread_cols{2}.LineWidth = 3; % Change marker width
    end
    ylabel('Posterior Probability','Fontsize',12);
    set(gca,'xticklabel',geno_list.colheaders,'Fontsize',12); % Name each group
end

% Stats 
% Active 
for c = 1:numComp(1) % For each active cluster
    if max(experiment_tags{1,1}) > 1 % With experiment tags 
        [twa.gf.active.p(1:3,c),~,twa.gf.active.stats{c}] = anovan...
            (p_assign_dist_mean{1,1}(:,c),{i_group_tags,i_experiment_tags},...
            'display','off','model','full');
    else % Without experiment tags 
        [twa.gf.active.p(1,c),~,twa.gf.active.stats{c}] = anovan...
            (p_assign_dist_mean{1,1}(:,c),{i_group_tags},...
            'display','off','model','full'); % Try without
    end
end

% Inactive 
for c = 1:numComp(2) % For each cluster
    if max(experiment_tags{2,1}) > 1 % With experiment tags 
        [twa.gf.inactive.p(1:3,c),~,twa.gf.inactive.stats{c}] = anovan...
            (p_assign_dist_mean{2,1}(:,c),{i_group_tags,i_experiment_tags},...
            'display','off','model','full');
    else % Without experiment tags 
        [twa.gf.inactive.p(1,c),~,twa.gf.inactive.stats{c}] = anovan...
            (p_assign_dist_mean{2,1}(:,c),{i_group_tags},...
            'display','off','model','full'); % Try without
    end
end

clear c e f spread_cols

%% Bout Proportions  
bout_proportions{1,1} = nan(max(fish_tags{1,1}),numComp(1),max(parameter_indicies{1,1}),...
    'single'); % fish x clusters x time windows 
bout_proportions{2,1} = nan(max(fish_tags{2,1}),numComp(2),max(parameter_indicies{2,1}),...
    'single'); % fish x clusters x time windows 

% Active
for f = 1:max(fish_tags{1,1}) % For each fish 
    for c = 1:numComp(1) % For each active bout type 
        for t = 1:max(parameter_indicies{1,1}) % For each time window 
            bout_proportions{1,1}(f,c,t) = size(find(fish_tags{1,1} == f & ...
                idx_numComp_sorted{1,1} == c & parameter_indicies{1,1} == t),1)/...
                    size(find(fish_tags{1,1} == f & parameter_indicies{1,1} == t),1); 
        end 
    end 
    disp(horzcat('Calculated Active Bout Proportions for fish ',num2str(f),...
        ' of ',num2str(max(fish_tags{1,1})))); 
end 

% Inactive 
for f = 1:max(fish_tags{2,1}) % For each fish 
    for c = 1:numComp(2) % For each inactive bout type 
        for t = 1:max(parameter_indicies{2,1}) % For each time window 
            bout_proportions{2,1}(f,c,t) = size(find(fish_tags{2,1} == f & ...
                idx_numComp_sorted{2,1} == c & parameter_indicies{2,1} == t),1)/...
                    size(find(fish_tags{2,1} == f & parameter_indicies{2,1} == t),1); 
        end 
    end 
    disp(horzcat('Calculated Inactive Bout Proportions for fish ',num2str(f),...
        ' of ',num2str(max(fish_tags{2,1})))); 
end 

clear f c t 

%% Bout Proportions Figures 
% Active 
figure; hold on;
for c = 1:numComp(1) % for each cluster
    subplot(3,5,c); title(horzcat('Active Cluster ',num2str(c))); hold on;
    counter = 1; clear scrap; % subplot
    for e = 1:max(i_experiment_tags) % for each experiment
        for g = 1:max(i_group_tags) % for each group
            legend_lines(g) = errorbar(nanmean(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2])),...
                nanstd(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]))/...
                sqrt(size(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]),1)),...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'linewidth',3);
            
            if e == 1 % For the first experiment
                legend_cell{g} = horzcat(geno_list.colheaders{g},', n = (',...
                    num2str(size(find(i_group_tags == g),1)),')');
            end
            
            scrap(1,counter) = max(nanmean(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2])) +...
                nanstd(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]))/...
                sqrt(size(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]),1)));
            scrap(2,counter) = min(nanmean(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2])) -...
                nanstd(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]))/...
                sqrt(size(permute(bout_proportions{1,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]),1)));
            
            counter = counter + 1; % add to counter
        end
    end
    
    % Add night patches
    y_lims = [(min(scrap(2,:)) - min(scrap(2,:))*0.05) ...
        (max(scrap(1,:)) + max(scrap(1,:))*0.05)]; % Add a bit of space either side
    
    a = 1; night_start = first_night; % Start counters
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
            1 (y_lims(2)-y_lims(1))],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1; night_start = night_start + 2; % Add to counters
    end
    
    % Figure Looks
    if c == 5 % For the 5th cluster
        [~,~,~,~] = legend(legend_cell,'Location','northeast'); % Generate axis
        legend('boxoff'); % Turn legend off
    end
    axis([0.5 size([days_crop(days) nights_crop(nights)],2)+0.5 ...
        y_lims]); % Set axis
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
    xlabel('Time (Days/Nights)','Fontsize',12); % X Labels
    set(gca, 'XTick', []); % Turn off X-Ticks
    ylabel('Proportion of Bouts','Fontsize',12); % Y Labels
    
end

% Inactive 
figure; hold on;
for c = 1:numComp(2) % for each cluster
    subplot(2,3,c); title(horzcat('Inactive Cluster ',num2str(c))); hold on;
    counter = 1; clear scrap; % subplot
    for e = 1:max(i_experiment_tags) % for each experiment
        for g = 1:max(i_group_tags) % for each group
            legend_lines(g) = errorbar(nanmean(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2])),...
                nanstd(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]))/...
                sqrt(size(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]),1)),...
                'color',cmap(g,:)+(1-cmap(g,:))*(1-(1/e^.5)),'linewidth',3);
            
            if e == 1 % For the first experiment
                legend_cell{g} = horzcat(geno_list.colheaders{g},', n = (',...
                    num2str(size(find(i_group_tags == g),1)),')');
            end
            
            scrap(1,counter) = max(nanmean(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2])) +...
                nanstd(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]))/...
                sqrt(size(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]),1)));
            scrap(2,counter) = min(nanmean(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2])) -...
                nanstd(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]))/...
                sqrt(size(permute(bout_proportions{2,1}(i_group_tags == g &...
                i_experiment_tags == e,c,time_window(1):time_window(2)),[1 3 2]),1)));
            
            counter = counter + 1; % add to counter
        end
    end
    
    % Add night patches
    y_lims = [(min(scrap(2,:)) - min(scrap(2,:))*0.05) ...
        (max(scrap(1,:)) + max(scrap(1,:))*0.05)]; % Add a bit of space either side
    
    a = 1; night_start = first_night; % Start counters
    for n = 1:size(nights,2) % For each night
        r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
            1 (y_lims(2)-y_lims(1))],...
            'FaceColor',night_color,'Edgecolor',[1 1 1]);
        uistack(r(a),'bottom'); % Send to back
        a = a + 1; night_start = night_start + 2; % Add to counters
    end
    
    % Figure Looks
    if c == 3 % For the 5th cluster
        [~,~,~,~] = legend(legend_cell,'Location','northeast'); % Generate axis
        legend('boxoff'); % Turn legend off
    end
    axis([0.5 size([days_crop(days) nights_crop(nights)],2)+0.5 ...
        y_lims]); % Set axis
    box off; set(gca, 'Layer','top'); set(gca,'Fontsize',12); % Format
    xlabel('Time (Days/Nights)','Fontsize',12); % X Labels
    set(gca, 'XTick', []); % Turn off X-Ticks
    ylabel('Proportion of Bouts','Fontsize',12); % Y Labels
    
end

clear c counter scrap e g a n y_lims r legend_cell legend_cols legend_lines

%% Bout Proportion Stats
% Grouping variables - note that these are the same for both states 
anova_group = repmat(i_group_tags,[size([days nights],2),1])';
anova_experiment = repmat(i_experiment_tags,[size([days nights],2),1])';
anova_time = []; 
for t = time_window(1):time_window(2) % For each time window 
    anova_time = [anova_time ; ones(size(i_experiment_tags,1),1)*mod(t,2)]; 
    % Allocate alternating zeros and ones to each time window 
end 
anova_time = anova_time(:)'; % vectorise 
if size(days_crop(days),2) == size(nights_crop(nights),2) % If there are an equal number of windows 
    anova_development = []; % development
    anova_development = zeros(1,size(anova_group,2)); % Pre-allocate 
    d = 1:size(anova_development,2)/(size(time_window(1):time_window(2),2)/2):...
        size(anova_development,2); % divide into "24h" windows 
    for t = 1:size(d,2)-1 
        anova_development(d(t):d(t+1)-1) = t; 
    end 
end 

% Calculation 
% Active 
for c = 1:numComp(1) % For each active cluster 
    clear scrap;
    scrap = permute(bout_proportions{1,1}(:,c,time_window(1):time_window(2)),[1 3 2]); 
    scrap = scrap(:)'; % Vectorise  
      
    if size(days_crop(days),2) == size(nights_crop(nights),2) % If comparing development 
        if max(experiment_tags{1,1}) > 1 % If comparing experiments 
            [twa.bp.active.p(:,c),~,twa.bp.active.stats{c}] = anovan(scrap,...
                {anova_group,anova_time,anova_development,anova_experiment},...
                'display','off','model','full');
        else % Development but no experiments 
            [twa.bp.active.p(1:7,c),~,twa.bp.active.stats{c}] = anovan(scrap,...
                {anova_group,anova_time,anova_development},...
                'display','off','model','full');
        end
    else % Without development 
        if max(experiment_tags{1,1}) > 1 % With experiments 
            [twa.bp.active.p(1:7,c),~,twa.bp.active.stats{c}] = anovan(scrap,...
                {anova_group,anova_time,anova_experiment},...
                'display','off','model','full'); % Try without
        else % Without experiments 
            [twa.bp.active.p(1:3,c),~,twa.bp.active.stats{c}] = anovan(scrap,...
                {anova_group,anova_time},...
                'display','off','model','full'); % Try without
        end     
    end
end

% Inactive 
for c = 1:numComp(2) % For each cluster 
    clear scrap;
    scrap = permute(bout_proportions{2,1}(:,c,time_window(1):time_window(2)),[1 3 2]); 
    scrap = scrap(:)'; % Vectorise  
      
    if size(days_crop(days),2) == size(nights_crop(nights),2) % If comparing development 
        if max(experiment_tags{2,1}) > 1 % If comparing experiments 
            [twa.bp.inactive.p(:,c),~,twa.bp.inactive.stats{c}] = anovan(scrap,...
                {anova_group,anova_time,anova_development,anova_experiment},...
                'display','off','model','full');
        else % Development but no experiments 
            [twa.bp.inactive.p(1:7,c),~,twa.bp.inactive.stats{c}] = anovan(scrap,...
                {anova_group,anova_time,anova_development},...
                'display','off','model','full');
        end
    else % Without development 
        if max(experiment_tags{2,1}) > 1 % With experiments 
            [twa.bp.inactive.p(1:7,c),~,twa.bp.inactive.stats{c}] = anovan(scrap,...
                {anova_group,anova_time,anova_experiment},...
                'display','off','model','full'); % Try without
        else % Without experiments 
            [twa.bp.inactive.p(1:3,c),~,twa.bp.inactive.stats{c}] = anovan(scrap,...
                {anova_group,anova_time},...
                'display','off','model','full'); % Try without
        end
        
    end
    
end

clear anova_experiment anova_group anova_time anova_development ...
    t d c scrap 

%% Cluster Parameters

% Active 
figure; hold on; 
for p = 1:6 % For each parameter 
    subplot(2,3,p); title(parameters{p}); 
    plotSpread(wake_cells(:,p+2)/unit_conversion(1,p),'distributionIdx',...
        idx_numComp_sorted{1,1},'distributionColors',cmap_cluster{1,1}); 
    set(gca,'xtick',[]); 
    xlabel('Clusters','Fontsize',12); 
    ylabel(units{p},'Fontsize',12); 
end 

% Inactive - scatter  
figure; hold on; 
for p = 1 % For each parameter 
    title(parameters{p}); 
    plotSpread(log10(sleep_cells(:,p+2)/unit_conversion(1,p)),'distributionIdx',...
        idx_numComp_sorted{2,1},'distributionColors',cmap_cluster{2,1}); 
    set(gca,'xtick',[]); 
    xlabel('Clusters','Fontsize',12); 
    ylabel(units{p},'Fontsize',12); 
end

%% Scattering bouts in PCA Space
% Active 
    % Should try coloring each cluster using posterior probability bins
figure; hold on; box off;  
for c = numComp(1):-1:1  % For each cluster 
    scatter((score(idx_numComp_sorted{1,1}==c,1)),score(idx_numComp_sorted{1,1}==c,2),...
        'markerfacecolor',cmap_cluster{1,1}(c,:),...
        'markeredgecolor',cmap_cluster{1,1}(c,:));
    xlabel('PC1','Fontsize',32); 
    ylabel('PC2','Fontsize',32); 
end 
set(gca,'Fontsize',18); 

% Add in ScatterHist to show the Gaussians 

%% Line Plots of each clusters feature 
scrap = zscore(wake_cells(:,3:end)); 
figure; hold on; 
for c = numComp(1):-1:1 % For each cluster 
    plot(nanmean(scrap(idx_numComp_sorted{1,1}==c,:)),'color',cmap_cluster{1,1}(c,:),'linewidth',3);
    set(gca,'Xtick',1:6);
    set(gca,'Fontsize',18); 
    set(gca,'XTickLabels',parameters(1:6),'Fontsize',12); 
    xlabel('Parameters','Fontsize',32); 
    ylabel('Z Score','Fontsize',32); 
end 
clear scrap; 