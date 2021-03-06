%% State_Space_5 

% V4 - Aims to use previous state space versions as a rough frame work to now 
    % work with data from Bed_frames (frame by frame data)
% V5 - Fit's a GMM to a sub-sample of the data, then assigns all data 
    % to cluster's based on posterior probabilities 
    
%% Required scripts 

% Lbmap 
%https://uk.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps

% Knee_pt 
%https://uk.mathworks.com/matlabcentral/fileexchange/35094-knee-point

%% Settings 
set(0,'DefaultFigureWindowStyle','docked'); % dock figures 
set(0,'defaultfigurecolor',[1 1 1]); % white background

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
    
    %delta_px_sq{1,f} = experiment.delta_px_sq; skipping this is easier on memory
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
% the mean from each column (as the means are not quite zero this seems
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

X{2,1} = []; % empty X
X{2,1} = sleep_cells(:,3); 

% tic
% X{2,1} = []; % empty X
% % Z-score each fishes data
% for f = 1:max(fish_tags{2,1}) % for each fish
%     X{2,1} = [X{2,1} ; zscore(sleep_cells(fish_tags{2,1} == f,3))];
%     if mod(f,100) == 0 % report every 100 fish 
%         disp(horzcat('Completed ',num2str(f),' fish of ',...
%             num2str(max(fish_tags{2,1}))));
%     end
% end
% toc

% Hard Coded Parameters
%ds = [0.01 0.1 1 5 10 15 25 50]; % set the number of data points to use for each clustering iteration
ds = 2000000; % two million
reps = 1; % set the number of repetitions
k_max = 1:20; % set values of k (clusters) to try
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
                
                tic % start timer
                GMModels{1} = fitgmdist(sample,k,...
                    'Options',options,'RegularizationValue',...
                    abs(rv),'Replicates',5); % Fit K gaussians
                % To the down-sampled data
                
                % Cluster the full data-set (using this mixture)
                [idx{1},~,P{1}] = cluster(GMModels{1},X{s,1});
                P{1} = max(P{1},[],2); % Keep only assigned probabilities (helps memory)
                
                time_sheet{s,1}(k,r,d) = toc; % record time for each iteration
                
                % Information Criteria
                AIC{s,1}(k,r,d)= GMModels{1}.AIC; % Extract AIC
                BIC{s,1}(k,r,d)= GMModels{1}.BIC; % Extract BIC
                
                % Store cluster centroids
                cluster_centroids{s,1}{k,r,d} = GMModels{1}.mu;                
                
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

clear sample sample_tags score_values score_zero rv ...
    GMModels idx P s d r k;

% Save Data
save(strcat(save_pathname,'\Full_Set.mat'),'-v7.3');

%% Choice of Sample Size 
    % Re-worked from Down_Sampling_Code_V4
    
% cmap_cluster{1,1} = lbmap(max(k_max),'RedBlue'); % generate a colormap for the clusters 
% 
% % DS Data (k,r,d)
% figure; 
% for d = 1:size(ds,2) % for each downsampling
%     subplot(2,4,d); hold on; set(gca,'Fontsize',12); 
%     title(horzcat(num2str(ds(d)*100),'% SubSamples'),'Fontsize',18); 
%     
%     for r = 1:reps % for each repetition
%         scatter(cluster_centroids{1,r,d}(:,1),cluster_centroids{1,r,d}(:,2),72,...
%             cmap_cluster{1,1},'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
%     end
%     
%     xlabel('PC 1','Fontsize',16); 
%     ylabel('PC 2','Fontsize',16); 
%     axis([-3 9 -3 9]);
%     
% end
% 
% % Comparison 
% cluster_centroids_comp = nan(size(ds,2),max(k_max),'single'); % ds x k 
% for d = 1:size(ds,2) % for each downsampling
%     for c = 1:max(k_max) % for each cluster
%         scrap = []; 
%         for r = 1:reps % for each repetition
%             scrap = [scrap ; cluster_centroids{1,r,d}(c,:)]; 
%         end
%         cluster_centroids_comp(d,c) = sum(pdist(scrap))/2;
%     end
% end
% 
% figure; subplot(2,1,1); set(gca,'Fontsize',12); hold on; 
% title('Variation in Centroid Position','Fontsize',18); 
% plot(1:size(ds,2),sum(cluster_centroids_comp,2),'color',cmap(1,:),'linewidth',3)
% ylabel('Total Distance (Euclidean)','Fontsize',16); 
% xticklabels({num2str(ds'*100)});
% xlabel('Percentage of Data Sampled','Fontsize',16); 
% 
% % Time Figure 
% ax = subplot(2,1,2); hold on; ax.FontSize = 12; 
% title('Time Taken for each Subsample','Fontsize',18);
% spread_cols = plotSpread(squeeze(time_sheet)/60,'showMM',4,'distributionColors',cmap_2(2,:)); 
% spread_cols{2}(1).LineWidth = 3; spread_cols{2}(2).LineWidth = 3; % Change marker width 
% spread_cols{2}(1).Color = cmap_2(1,:); spread_cols{2}(2).Color = cmap_2(1,:); 
% xticklabels({num2str(ds'*100)});
% xlabel('Percentage of Data Sampled','Fontsize',16); 
% ylabel('Time Taken (minutes)','Fontsize',16); 

%% Choice of Number of Clusters & Settings  

% Load data
load('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\New\Z_Score_Wake\Full_Set.mat'); 
strings{1,1} = 'Active'; strings{2,1} = 'Inactive'; % name tags 

% Info Criteria Figure 
figure; 
for s = 1:2 % for active & inactive
    ax(s) = subplot(1,2,s); set(gca,'FontName','Calibri'); hold on; title(strings{s,1});
    plot(BIC{s,1},'linewidth',3);
    box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
    if s == 2 % for inactive 
        [~,icons,plots,~] = legend('BIC','Location','northeast');
        legend('boxoff');
        set(icons(1),'Fontsize',26) ; set(plots,'LineWidth',3);
    end
    xlabel('Number of Clusters','Fontsize',26);
    ylabel('Criteria Value','Fontsize',26);
end

% Choice of number of clusters 
numComp = [10 6]; % Choose active & inactive numbers of clusters
scrap = lbmap(sum(numComp),'RedBlue'); % Color Scheme
cmap_cluster{1,1} = scrap(1:numComp(1),:); % generate a colormap for the clusters 
cmap_cluster{2,1} = scrap(numComp(1)+1:end,:); % generate a colormap for the clusters 

% Add to Figure 
for s = 1:2
    subplot(ax(s)); scatter(numComp(s),BIC{s,1}(numComp(s)),72,'k','filled'); 
end 

% Load Data 
Wake = load(horzcat('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\New\Z_Score_Wake\1s_2000000d_1r_',...
    num2str(numComp(1)),'k.mat')); 
Sleep = load(horzcat('D:\Behaviour\SleepWake\Re_Runs\Clustered_Data\New\Z_Score_Wake\2s_2000000d_1r_',...
    num2str(numComp(2)),'k.mat')); 

% Remove NaN Values from Sleep matracies  
sleep_cells(sleep_cells_nan_track,3) = NaN; % sleep cells 
Sleep.idx{1}(sleep_cells_nan_track,1) = NaN; % cluster assignments 
Sleep.P{1}(sleep_cells_nan_track,1) = NaN; % posterior probabilities 

% Merge Variables for ease of Looping 
idx{1,1} = Wake.idx{1,1}; idx{2,1} = Sleep.idx{1,1}; % clusters  
P{1,1} = Wake.P{1,1}; P{2,1} = Sleep.P{1,1}; % posterior probabilities 
cells{1,1} = wake_cells; cells{2,1} = sleep_cells; % cells 

% Grouping repeats of experiments (hard coded)
experiment_reps = [1 1 1 2 2 3 3 4 5]; % experiment groupings 
i_experiment_reps = i_experiment_tags;
for er = 1:max(experiment_reps) % for each repeat 
    found = find(experiment_reps == er); % find experiments
    
    for f = found % for each experiment in the repeat 
        i_experiment_reps(i_experiment_reps == f,1) = er; % tag with grouping variable 
    end 
    
end 

% Adjust colours 
for e = 1:size(cmap,2) % for each experiment
    if max(i_group_tags(i_experiment_tags==e)) == 1 % if theres just one group (e.g. WT experiments)
        cmap{e}(1,:) = [135 206 250]/255; % light sky blue
        cmap_2{e}(1,:) = cmap{e};
        cmap_2{e}(2,:) = [25 25 112]/255; % midnight blue
    else % for experiments with multiple groups 
        cmap_2{e} = flip(cmap_2{e}); % flip cmap (so it starts with blue)
        for c = 1:2:size(cmap_2{e},1) % for every other color 
            cmap_2{e}([c c+1],:) = cmap_2{e}([c+1 c],:); % swap pairs of colours around 
        end
        cmap{e} = cmap_2{e}(1:2:size(cmap_2{e},1),:); % Extract main colors
    end
end

% Adjust time windows (Hard coded) 
for e = 1:size(time_window,2) % for each experiment 
    if experiment_reps(e) < 4 
       time_window{e} = [3 6]; % take the middle two days/nights 
       days{e} = [2 3]; nights{e} = [2 3]; 
    else
       time_window{e} = [1 2]; % take the first day/night 
       days{e} = 1;
    end 
end 

clear s ax icons plots scrap Wake Sleep er f e c

%% Sorting Clusters by Mean Length 
             
for s = 1:2 % for active & inactive
    mean_cluster_length = nan(1,numComp(s),'single'); % pre-allocate
    for c = 1:numComp(s) % For each cluster
        mean_cluster_length(c) = nanmean(cells{s,1}(idx{s,1}==c,3));
        % Calculate mean bout length
    end
    
    [~,O] = sort(mean_cluster_length); % Sort by increasing bout length
    
    idx_numComp_sorted{s,1} = nan(size(cells{s,1},1),1,'single'); % Pre-allocate
    
    for c = 1:numComp(s)  % For each cluster
        idx_numComp_sorted{s,1}(idx{s,1} == O(c),:) = c;
        % Re-assign cluster numbers
    end
    
    clear s mean_cluster_length c O;
end

clear idx

%% Cluster Parameters - Fits 

tic
parameter_dists = cell(2,1); % structure {s,parameters}(cluster,min:max(parmater:))

for s = 1:2 % for active & inactive
    for p = 3:size(cells{s,1},2) % for each parameter
        for k = 1:numComp(s) % for each cluster
            clear pd;
            pd = fitdist(cells{s,1}(idx_numComp_sorted{s,1}==k,p),'kernel','Width',1); % Fit
            parameter_dists{s,p-2}(k,:) = pdf(pd,min(cells{s,1}(:,p)):max(cells{s,1}(:,p)));
        end
    end
end
toc

clear s p k pd 

%% Cluster Parameters - Fits Figure 

figure;
counter = 1; % start a counter 
for s = 1:2 % for active & inactive
    for p = 1:size(cells{s,1},2) - 2 % For each parameter
        % Figure Settings
        subplot(2,4,counter); hold on;
        box off; set(gca,'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
        if s ~= 2 % for active parameters 
            title(parameters{1}{p}); % Add title
        else % for inactive bout length 
            title(parameters{1}{10}); % Add title   
        end
        
        % Plot 
        crop = size(parameter_dists{s,p},2); 
        for k = 1:numComp(s) % for each cluster
            plot((1:crop)/unit_conversion{1}(s,p),parameter_dists{s,p}(k,:),'color',cmap_cluster{s,1}(k,:),'linewidth',3)
        end
        
        % Axes 
        set(gca,'XTick',...
            [1/unit_conversion{1}(s,p), 1/unit_conversion{1}(s,p)*10,...
            crop/unit_conversion{1}(s,p)]); % set x tick labels
        
        axis([1/unit_conversion{1}(s,p) crop/unit_conversion{1}(s,p) ...
            min(min(parameter_dists{s,p})) max(max(parameter_dists{s,p}))]); % Set axis limits
        
        % Set decimal places depending on units
        if unit_conversion{1}(s,p) > 1
            xtickformat('%.2f');
        else
            xtickformat('%.0f');
        end
        
        set(gca,'XScale','log'); % set log axis
        xlabel(units{1}(p),'Fontsize',12); % X labels
        ylabel('Probability Density','Fontsize',12); % Y label
        
        counter = counter + 1; % add to counter 
        
        % Legend Plot 
        if counter == 8 
            subplot(2,4,counter); hold on;
            box off; set(gca,'Layer','top'); set(gca,'Fontsize',12); set(gca,'FontName','Calibri'); % Set Font
            title('Legend'); 
            
            for s_2 = 1:2 % for active & inactive 
                for k = 1:numComp(s_2) % for each cluster 
                    scatter(k,mod(s_2,2),300,...
                        'markerfacecolor',cmap_cluster{s_2,1}(k,:),...
                        'markeredgecolor',cmap_cluster{s_2,1}(k,:));
                end 
            end 
            
            axis([0 max(numComp)+1 -1 2]);
            xlabel('Cluster Number','Fontsize',12); % X labels
            set(gca,'YTick',[0 1]); 
            set(gca,'YTickLabels',{'Inactive','Active'},'Fontsize',12);
        end 
        
    end
end

clear counter s p crop k s_2

%% Posterior Probabilities 

% Fit  
pd = fitdist(P{1,1},'Kernel','Width',0.05); % Fit kernel distribution 
p_dist{1,1} = pdf(pd,0:0.05:1); clear pd; % Save as pdf  
pd = fitdist(P{2,1},'Kernel','Width',0.05); % Fit kernel distribution
p_dist{2,1} = pdf(pd,0:0.05:1); % Save as pdf 

% Figure
figure; hold on; set(gca,'FontName','Calibri'); title('Posterior Probabilities'); 
plot(0:0.05:1,p_dist{1,1},'color',cmap_cluster{1,1}(1,:),'linewidth',3)
plot(0:0.05:1,p_dist{2,1},'color',cmap_cluster{2,1}(1,:),'linewidth',3)

box off; set(gca,'Layer','top'); set(gca,'Fontsize',32);
[~,icons,plots,~] = legend('Active Bouts','Inactive Bouts','Location','northwest');
legend('boxoff'); 
set(icons(1:2),'Fontsize',26) ; set(plots,'LineWidth',3);
xlabel('Posterior Probability','Fontsize',26);
ylabel('Probability Density','Fontsize',26);

clear pd icons plots 

%% Bout Proportions (slow...)

tic
bout_proportions{1,1} = nan(max(fish_tags{1,1}),numComp(1),max(parameter_indicies{1,1}),...
    'single'); % fish x clusters x time windows 
bout_proportions{2,1} = nan(max(fish_tags{2,1}),numComp(2),max(parameter_indicies{2,1}),...
    'single'); % fish x clusters x time windows 

for s = 1:2 % for active & inactive  
    % note that comms overhead makes this faster as a for rather than a parfor loop  
    for f = 1:max(fish_tags{s,1}) % For each fish
        for c = 1:numComp(s) % For each active bout type
            for t = 1:max(parameter_indicies{s,1}(fish_tags{s,1}==f)) % For each time window that fish uses 
                bout_proportions{s,1}(f,c,t) = sum(fish_tags{s,1}==f & idx_numComp_sorted{s,1}==c ...
                    & parameter_indicies{s,1}==t)/...
                    sum(fish_tags{s,1}==f & parameter_indicies{s,1}==t); 
                % the number of times fish (f) uses cluster (c) time (t) 
                % divide by the number of bouts fish (f) has that time (t) 
                
                % Note - will return zero's when a fish doesn't use a
                % particular bout type :-)
            end
        end
    end
end
toc

clear s f c t 

%% Bout Proportions Figure

for er = 1:max(experiment_reps) % for each group of experiments
    figure; counter = 1; % start counter to track clusters 
    for s = 1:2 % for active & inactive
        for c = 1:numComp(s) % for each cluster
            
            counter_2 = 1; % start a counter to track color shades 
            subplot(4,4,counter); set(gca,'FontName','Calibri'); ...
                title(horzcat(strings{s},' Cluster ',num2str(c))); hold on;
            
            for e = find(experiment_reps == er) % for each experiment in this group (er) 
                for g = 1:max(i_group_tags(i_experiment_tags == e)) % for each group
                    legend_lines(g) = errorbar(nanmean(permute(bout_proportions{s,1}(i_experiment_tags == e & ...
                        i_group_tags == g,c,time_window{e}(1):time_window{e}(2)),[1 3 2])),...
                        (nanstd(permute(bout_proportions{s,1}(i_experiment_tags == e & ...
                        i_group_tags == g,c,time_window{e}(1):time_window{e}(2)),[1 3 2]))/sqrt(sum(i_experiment_tags==e & i_group_tags==g))),...
                        'color',cmap{e}(g,:)+(1-cmap{e}(g,:))*(1-(1/counter_2^.5)),'linewidth',3);
                    
                    if e == find(experiment_reps == er,1,'first') % For the first experiment
                        legend_cell{g} = geno_list{e}.colheaders{g}; % add to legend 
                    end
                    
                    % determining good axis boundaries 
                    scrap(counter_2,g,1) = min(nanmean(permute(bout_proportions{s,1}(i_experiment_tags == e & ...
                        i_group_tags == g,c,time_window{e}(1):time_window{e}(2)),[1 3 2])) - ...
                        (nanstd(permute(bout_proportions{s,1}(i_experiment_tags == e & ...
                        i_group_tags == g,c,time_window{e}(1):time_window{e}(2)),[1 3 2]))/sqrt(sum(i_experiment_tags==e & i_group_tags==g)))); 
                    
                    scrap(counter_2,g,2) = max(nanmean(permute(bout_proportions{s,1}(i_experiment_tags == e & ...
                        i_group_tags == g,c,time_window{e}(1):time_window{e}(2)),[1 3 2])) + ...
                        (nanstd(permute(bout_proportions{s,1}(i_experiment_tags == e & ...
                        i_group_tags == g,c,time_window{e}(1):time_window{e}(2)),[1 3 2]))/sqrt(sum(i_experiment_tags==e & i_group_tags==g)))); 
                end
                counter_2 = counter_2 + 1;
            end
                        
            % Add night patches
            y_lims = [(min(min(scrap(:,:,1))) - (min(min(scrap(:,:,1)))*0.05)) ...
                (max(max(scrap(:,:,2))) + (max(max(scrap(:,:,2)))*0.05))]; % Add a bit of space either side
            
            a = 1; night_start = first_night{e}; % Start counters
            for n = 1:size(nights{e},2) % For each night
                r(a) = rectangle('Position',[(night_start-0.5) y_lims(1)...
                    1 (y_lims(2)-y_lims(1))],...
                    'FaceColor',night_color{e},'Edgecolor',[1 1 1]);
                uistack(r(a),'bottom'); % Send to back
                a = a + 1; night_start = night_start + 2; % Add to counters
            end
            
            % Figure Looks
            if c == 4 && s == 1 % For the 4th active cluster
                [~,~,~,~] = legend(legend_cell,'Location','best'); % Generate axis
                legend('boxoff'); % Turn legend off
            end
            
            axis([0.5 size([days_crop{e}(days{e}) nights_crop{e}(nights{e})],2)+0.5 ...
                y_lims]); % Set axis
            
            box off; set(gca,'Layer','top'); set(gca,'Fontsize',12); % Format
            xlabel('Time (Days/Nights)','Fontsize',12); % X Labels
            set(gca, 'XTick', []); % Turn off X-Ticks
            ylabel('Proportion of Bouts','Fontsize',12); % Y Labels
            
            counter = counter + 1; 
            
            clear scrap y_lims legend_lines legend_cell r 
        end
    end
end

clear er counter s c counter_2 e g legend_lines legend_cell scrap ...
    y_lims a night_start n r 

%% Bout Proportion Stats

for er = 1:max(experiment_reps) % for each group of experiments
    set_token = find(experiment_reps == er,1,'first'); % used for each experiments sets settings 
    
    for s = 1:2 % for active & inactive
        
        % Grouping Variables
        anova_group = repmat(i_group_tags(i_experiment_reps==er),...
            [size([days{set_token} nights{set_token}],2),1])'; % groups
        anova_experiment = repmat(i_experiment_tags(i_experiment_reps==er),...
            [size([days{set_token} nights{set_token}],2),1])'; % experiments
        
        anova_time = [];
        for t = time_window{set_token}(1):time_window{set_token}(2) % For each time window
            anova_time = [anova_time ; ones(sum(i_experiment_reps==er),1)*mod(t,2)];
            % Allocate alternating zeros and ones to each time window
        end
        anova_time = anova_time';
        
        % Development Grouping Variable
        if size(days_crop{set_token}(days{set_token}),2) == ...
                size(nights_crop{set_token}(nights{set_token}),2) ...
                && size(days_crop{set_token}(days{set_token}),2) > 1 % If there are an equal number of windows (>1)
            
            anova_development = []; % development
            anova_development = zeros(1,size(anova_group,2)); % Pre-allocate
            d = 1:size(anova_development,2)/(size(time_window{set_token}(1):...
                time_window{set_token}(2),2)/2):...
                size(anova_development,2); % divide into "24h" windows
            for t = 1:size(d,2)-1
                anova_development(d(t):d(t+1)-1) = t;
            end     
        else 
            anova_development = ones(size(anova_experiment)); % use all ones  
        end
        
        % Comparison
        for c = 1:numComp(s) % For each active cluster
            clear scrap;
            scrap = permute(bout_proportions{s,1}(i_experiment_reps==er,...
                c,time_window{set_token}(1):time_window{set_token}(2)),[1 3 2]);
            scrap = scrap(:)'; % Vectorise
            
            [twa.bp.p{s,er}(:,c),~,twa.bp.stats{s,er,c}] = anovan(scrap,...
                {anova_group,anova_time,anova_development,anova_experiment},...
                'display','off','model','full');
        end
        
        clear anova_development anova_experiment anova_group anova_time ... 
            scrap 
    end
end

clear er set_token s anova_group anova_experiment anova_time anova_development ...
    c scrap 

%% Mean "Goodness" of fit for each cluster

p_dist_mean{1,1} = nan(max(fish_tags{1,1}),numComp(1),'single'); % fish x clusters  
p_dist_mean{2,1} = nan(max(fish_tags{2,1}),numComp(2),'single'); % fish x clusters  
    % Note that using NaN's means that fish who don't have a particular
    % bout type won't contribute to the mean fit :-) 
tic
parfor s = 1:2 % active & inactive
    for f = 1:max(fish_tags{s,1}) % For each fish
        for c = 1:numComp(s) % For each cluster
            p_dist_mean{s,1}(f,c) = nanmean(P{s,1}...
                (idx_numComp_sorted{s,1} == c & fish_tags{s,1} == f));
            % the mean posterior probability for each cluster & fish 
        end
    end 
end
toc

clear s f c 

%% Mean "Goodness" of fit Figure 

for er = 1:max(experiment_reps) % for each group of experiments
    figure; counter = 1; % start subplot counter 
    for s = 1:2 % for active & inactive
        for c = 1:numComp(s) % for each cluster
            
            counter_2 = 1; % start color shade counter
            subplot(4,4,counter); set(gca,'FontName','Calibri'); ...
                title(horzcat(strings{s},' Cluster ',num2str(c))); hold on;
            
            for e = find(experiment_reps == er) % for each experiment
                spread_cols = plotSpread(p_dist_mean{s,1}(i_experiment_tags == e,c),...
                    'distributionIdx',i_group_tags(i_experiment_tags == e),...
                    'distributionColors',cmap{e}+(1-cmap{e})*(1-(1/counter_2^.5)),'showMM',2);
                set(findall(er,'type','line'),'markersize',15); % change marker sizes  
                spread_cols{2}.LineWidth = 3; spread_cols{2}.Color = 'k'; % Change Mean properties
                spread_cols{2}.MarkerSize = 12; 
                counter_2 = counter_2 + 1; % add to color shade counter 
            end

            ylabel('Posterior Probability','Fontsize',12);
            set(gca,'xticklabel',geno_list{e}.colheaders,'Fontsize',12); % Name each group
            counter = counter + 1; % add to subplot counter 
        end
    end
end

clear er counter s c counter_2 e spread_cols 

%% Mean "Goodness" of fit Stats 

for er = 1:max(experiment_reps) % for each group of experiments
    for s = 1:2 % for active & inactive
        for c = 1:numComp(s) % For each cluster
            [twa.gf.p{s,er}(1:3,c),~,twa.gf.stats{s,er,c}] = anovan...
                (p_dist_mean{s,1}(i_experiment_reps==er,c),...
                {i_group_tags(i_experiment_reps==er),i_experiment_tags(i_experiment_reps==er)},...
                'display','off','model','full');
        end      
    end
end

clear er s c 

%% OLD onwards 

% Scattering bouts in PCA Space
% Active 
    % Should try coloring each cluster using posterior probability bins
figure; hold on; box off;  
for c = numComp(1):-1:1  % For each cluster 
    scatter((score(idx_numComp_sorted{1,1}==c,1)),score(idx_numComp_sorted{1,1}==c,2),'.',...
        'markerfacecolor',cmap_cluster{1,1}(c,:),...
        'markeredgecolor',cmap_cluster{1,1}(c,:));
end 
xlabel('PC1','Fontsize',32);
ylabel('PC2','Fontsize',32);
set(gca,'Fontsize',18);

% Add in ScatterHist to show the Gaussians 
figure; 
scatterhist(score(:,1),score(:,2),'group',idx_numComp_sorted{1,1},...
    'kernel','on','color',cmap_cluster{1,1},'legend','on',...
    'linewidth',3); 

%% Line Plots of each clusters feature 
clear scrap; scrap = zscore(wake_cells(:,3:end)); 
figure; hold on; 
for c = numComp(1):-1:1 % For each cluster 
    plot(nanmean(scrap(idx_numComp_sorted{1,1}==c,:)),'color',cmap_cluster{1,1}(c,:),'linewidth',3);
    set(gca,'Xtick',1:6);
    set(gca,'Fontsize',18); 
    set(gca,'XTickLabels',parameters{1}(1:6),'Fontsize',12); 
    xlabel('Parameters','Fontsize',32); 
    ylabel('Z Score','Fontsize',32); 
end 
clear scrap; 

%% Bout Data Figures 
    % Should sort by something (e.g. length)

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