%% New Drug Analysis 

%% New FingerPrinting 

% Hard Coded Variables 
starting = 2; % Choose your first day (to use for the fingerprint) 
Cntrl = 1; % Define which group your WT controls are 

%Load Data - using multiselect
[filename, pathname] = uigetfile('*.mat', 'Select Geno files','MultiSelect','on'); %Select a geno file
if isequal(filename,0) %If no file is selected
    error('No File Selected') %Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) %Show selected filename
end;

for f = 1:size(filename,2) %For each file
    load(strcat(pathname,filename{f})); %Load the mat file
    multi_geno{f} = geno; %Store it in a strucutre called data
end;
clear pathname filename geno f;

% Finger Printing 
% Variables 1-24 
for e = 1:size(multi_geno,2) % For each experiment
    parameter_names = fieldnames(multi_geno{e}.summarytable.mean); % Find parameter Names
    a = 1; % Set Counter 
    
    for v = 1:size(parameter_names,1) % For each variable + the day/night means
        time_names = fieldnames(multi_geno{e}.summarytable.mean.(parameter_names{v}));
        
        for c = starting:starting+1 % For each period 
            for t = 2:-1:1 % For each night and day
                for g = 1:size(multi_geno{e}.data,2) % For each group
                    finger_prints_mean{e}(g,a) = ... % Take the mean 
                        multi_geno{e}.summarytable.mean.(parameter_names{v}).(time_names{t}){g}(c);
                    finger_prints_std{e}(g,a) = ...  % Take the Std 
                        multi_geno{e}.summarytable.std.(parameter_names{v}).(time_names{t}){g}(c);
                end
                    a = a + 1; % Add to Counter 
                                  
            end
        end
    end
    
end

% Variables 25/26 
for e = 1:size(multi_geno,2) % For each experiment
    for g = 1:size(multi_geno{e}.data,2) % For each group
        % Day 
        clear scrap; scrap = ... % Take the 
            multi_geno{e}.summarytable.averageWaking.day{g}(starting:starting+1,:); 
        finger_prints_mean{e}(g,25) = nanmean(scrap(:)');  
        finger_prints_std{e}(g,25) = nanstd(scrap(:)'); 
        
        % Night 
        clear scrap; scrap = ...
            multi_geno{e}.summarytable.averageWaking.night{g}(starting:starting+1,:); 
        finger_prints_mean{e}(g,26) = nanmean(scrap(:)'); 
        finger_prints_std{e}(g,26) = nanstd(scrap(:)'); 
    end
end 

%Z scores for every group and variable 
for e = 1:size(multi_geno,2) % For each experiment
    for g = 1:size(multi_geno{e}.data,2) % For each group
        for v = 1:size(finger_prints_mean{e},2) % For each parameter 
            z_scores{e}(g,v) = (finger_prints_mean{e}(g,v) - ...
                finger_prints_mean{e}(Cntrl,v))/finger_prints_std{e}(Cntrl,v); 
            % Calculate the difference divided by the WT Std
        end 
    end 
end

clear a c e g t v scrap parameter_names time_names Cntrl ...
    finger_prints_mean finger_prints_std multi_geno starting
%% Adding to the DataBase 
Het = []; Hom = []; 
for i = 1:size(z_scores,2) 
    Het = [Het ; z_scores{i}(2,:)]; 
    Hom = [Hom ; z_scores{i}(3,:)]; 
end 

figure; hold on
plot(nanmean(Het)); plot(nanmean(Hom));  

data = [data ; Het ; nanmean(Het) ; Hom ; nanmean(Hom)]; 

clear Het Hom i clear z_scores 

%% Clustering etc

load('C:\Users\Marcus\Documents\MATLAB\WorkSpaces\Drug_DataBase_Mutants.mat')

% Pca 
X = zscore(giantdata);  

[coeff,score,latent,tsquared,explained,mu] = pca(X); 

[knee_dim] = knee_pt(explained) % Choose this many dimensions 

% Clustering in PCA Space with Gaussian Mixed Models  

% Set-up 
AIC = zeros(1,size(score,1)-1); 
GMModels = cell(size(score,1)-1,1);
idx = cell(size(score,1)-1,1);
nlogl = cell(size(score,1)-1,1);
P = cell(size(score,1)-1,1);

options = statset('MaxIter',1000);

score_values = unique(score(:,1:knee_dim));
score_zero = knnsearch(score_values,0);

progressbar('Clusters') %Initialise progress bars
for k = 1:size(score,1)-1 % Try all possible no of clusters 
    
    GMModels{k} = fitgmdist(score(:,1:knee_dim),k,...
        'Options',options,'RegularizationValue',...
        abs(score_values(score_zero)),'Replicates',10); % Fit K gaussians
    
    AIC(k)= GMModels{k}.AIC; % Extract AIC 
    
    % Cluster using this mixing  
    [idx{k},nlogl{k},P{k}] = cluster(GMModels{k},score(:,1:knee_dim)); 
    
    progressbar(k/size(score,1)-1); %Update the Clustering progressbar
    k % Feedback 
end

save('C:\Users\Marcus\Documents\MATLAB\WorkSpaces\Drug_DataBase_Mutants_Clustering.mat')

%% Choice of Number of Clusters 
[minAIC,numComp] = min(AIC); % AIC 

for k = 1:size(GMModels,1) % BIC - note should add this above later 
    BIC(k)= GMModels{k}.BIC; % Extract BIC 
end 
[minAIC,numComp] = min(BIC); % BIC 

% Choose based on number of members 
cluster_size(1:size(GMModels,1),1:size(GMModels,1)) = NaN; % Pre-allocate
for m = 1:size(GMModels,1) % For each model  
    
    for c = 1:m % For each cluster 
        cluster_size(m,c) = size(find(idx{m}==c),1); % Find the number of 
            % Members 
    end 
    
end 
[minAIC,numComp] = knee_pt(nanstd(cluster_size')) % Choose a number of 
    % clusters with the lowest deviation in cluster size 

% Choose based on assigned probabilities 
p_assign_dist(1:size(P,1),1:size(P,1)) = NaN; % Pre-allocated

for k = 1:size(P,1) % For each clustering 
    for c = 1:size(score,1) % For each compound 
        p_assign_dist(c,k) = double(P{k}(c,idx{k}(c,1))); % Pull posterior probabilities
            % For each cluster assignment 
    end 
end 

%% Filtering 
% Filtering Out Compounds with poor negative posterior probabilities
    % Leave in all mutants 
p_assign(1:size(score,1)) = NaN; % Pre-allocated

for c = 1:size(score,1) % For each compound 
    p_assign(c) = P{numComp}(c,idx{numComp}(c,1)); % Pull posterior probabilities
        % For each cluster assignment 
end 

p_assign_cut_off = prctile(p_assign(1:5756),[75]); % Specify a probability cut off 
removed_compounds = find(p_assign(1:5756) < p_assign_cut_off); 
idx{numComp}(removed_compounds) = NaN; 

% Filter Out Cluster with less than 3 members 
for k = 1:max(idx{numComp}) % For each cluster 
    found = find(idx{numComp} == k); % Find it's members 
    cluster_size_filtered(k) = size(found,1); % Find the number of members
    
    if cluster_size_filtered(k) < 3 % If there are less than three members
       idx{numComp}(found) = NaN; % Remove them from the panel
       cluster_size_filtered(k) = NaN; 
       removed_compounds = [removed_compounds found']; 
    end 
    clear found 
end 

for r = removed_compounds
    Giant_Panel_Targets{r} = 'Filtered Out'; % Remove from targets
    Structural_Clusters(r) = NaN; % Remove From Structural Clusters 
end 

%% Re-assign Cluster Numbers 

cluster_numbers = unique(idx{numComp}); % Find remaining Clusters
cluster_numbers(isnan(cluster_numbers)) = []; % Remove NaN values 

for k = 1:size(cluster_numbers,1) % For each remaining cluster 
    found = find(idx{numComp} == cluster_numbers(k)); % Find cluster members 
    
    idx{numComp}(found) = k; % Adjust cluster numbering 
    
end 

cluster_size_filtered(isnan(cluster_size_filtered )) = []; % Remove NaNs 

%% Figures 
% Scatter 
cmap = jet(max(idx{numComp})); %For each group 
figure; hold on; 
% subplot(1,2,1); hold on
% scatter(1:size(objective_f,2),objective_f,'filled');
% scatter(knee_clusters,objective_f(knee_clusters),'r','filled');

for c = 1:max(idx{numComp}) % For each cluster 
    %subplot(1,2,2); hold on
    scatter(score(idx{numComp}==c,1),score(idx{numComp}==c,2),...
        'markerfacecolor',cmap(c,:),...
        'markeredgecolor',cmap(c,:));
end 

% Plot 
cmap = jet(max(idx{numComp})); %For each group 
figure; hold on; 
% subplot(1,2,1); hold on
% scatter(1:size(objective_f,2),objective_f,'filled');
% scatter(knee_clusters,objective_f(knee_clusters),'r','filled');

for c = 1:max(idx{numComp}) % For each cluster 
    subplot(ceil(max(idx{numComp})/10),10,c)
    plot(X(idx{numComp}==c,:)','color',cmap(c,:),'linewidth',1);
    axis([0 26 ylim]); 
    axis off
end 

% Organised Imagesc 
% organised_prints = []; 
% for c = 1:max(idx{numComp}) 
%     if c > 1
%     organised_prints = [organised_prints ; ...
%         nanmean(score(idx{numComp}==c,1:5))]; 
%     end
% end 

%% Tsne Figure 
% Settings 
    % Reduce to 2 dimensions (for visulaisation) 
    % Start with PCA to your previously determined number of dimensions
    % Try all possible perplexities (note that perplexity should be set to 
    % lower than your number of points 
    % Repeat 3 times for each perplexity 

    mappedX = cell(size(3:30:size(X,1),2));
    a = 1; 
    for per = 3:30:size(X,1) % For a range of possible perplexities
        mappedX{a} = tsne(X,[],2,knee_dim,per);
        a = a + 1 % Feedback 
    end
    
save('C:\Users\Marcus\Documents\MATLAB\WorkSpaces\Drug_DataBase_Mutants_Clustering.mat')

% Maximise distance between cluster centroids 
for per = 1:size(mappedX,2) % For a range of possible perplexities
        for c = 1:max(idx{numComp}) % For each cluster 
            mappedX_centroids_x(per,c) = nanmean(mappedX{per}...
                (idx{numComp}==c,1));
            mappedX_centroids_y(per,c) = nanmean(mappedX{per}...
                (idx{numComp}==c,2)); 
        end 
end 

for per = 1:size(mappedX,2) % For a range of possible perplexities
    a = 1;
    for m = 1:max(idx{numComp}) % For each cluster
        for n = m+1:max(idx{numComp}) % For each cluster comparison
            mappedX_of(per,a) = pdist2([mappedX_centroids_x(per,m),...
                mappedX_centroids_y(per,m)],...
                [mappedX_centroids_x(per,n),mappedX_centroids_y(per,n)]);
            a = a + 1;
        end
    end
end

% Alternatively Prioritise Tightness of Clustering 
for per = 1:size(mappedX,2) % For each perplexity 
    scrap = 0; clear found; 
    for c = 1:max(idx{numComp}) % For each cluster 
        found = mappedX{per}(idx{numComp}==c,:); 
        
        for m = 1:size(found,1) % For each member
            scrap = scrap + pdist2([mappedX_centroids_x(per,c),...
                mappedX_centroids_y(per,c)],...
                [found(m,1),found(m,2)]); 
        end 
        mappedX_dist(per,c) = scrap; 
        scrap = 0; clear found; 
    end 
    
end 

% Choice of Perplexity 
% [knee_perplexity] = knee_pt(nanmedian(mappedX_of,2)) % Distance
[knee_perplexity] = knee_pt(nanmedian(mappedX_dist')) % Tightness 

% Figure
cmap = jet(max(idx{numComp})); %For each cluster 

figure; 
scatter(mappedX{knee_perplexity}(:,1),mappedX{knee_perplexity}(:,2),144,...
    'markerfacecolor',[0.4392 0.5020 0.5647],...
    'markeredgecolor',[0.4392 0.5020 0.5647]); hold on; 

for c = 1:max(idx{numComp}) % For each cluster 
    scatter(mappedX{knee_perplexity}(idx{numComp}==c,1),...
        mappedX{knee_perplexity}(idx{numComp}==c,2),144,...
        'markerfacecolor',cmap(c,:),...
        'markeredgecolor',cmap(c,:)); hold on; 
end 

% for c = 1:max(idx{numComp}) % For each cluster 
%     text(mappedX_centroids_x(knee_perplexity,c),...
%         mappedX_centroids_y(knee_perplexity,c),num2str(c),...
%     'fontsize',22,'FontWeight','bold')
% end 
  
% Pharmacological Targets in Clusters 
for c = 1:max(idx{numComp}) % For each cluster 
    pharm_clusters{c} = Giant_Panel_Targets(idx{numComp}==c); 
end 

% Coloring tSne By Pharmacological target 
targets = unique(string(Giant_Panel_Targets)); % Find Unique Targets 
targets = rmmissing(targets); % Remove Missing targets 
targets = cellstr(targets); % Convert to Cell Array 
found = contains(string(targets),'Filtered Out'); % Find Filtered Compounds
targets(found) = []; % Remove these compounds
targets = targets(~cellfun('isempty',targets)); % Remove Empty Cells 

% Figure
cmap = jet(size(targets,1)); %For each pharmacological Target  

figure; hold on;
scatter(mappedX{knee_perplexity}(:,1),mappedX{knee_perplexity}(:,2),144,...
    'markerfacecolor',[0.4392 0.5020 0.5647],...
    'markeredgecolor',[0.4392 0.5020 0.5647]); 
for c = 1:size(targets,1) % For each target
    found = contains(string(Giant_Panel_Targets),targets{c}); % Find the indicies of each target
    % Note that this is a logical (unlike strfind)
    if sum(found) > 1
%         found_centroid = [nanmean(mappedX{knee_perplexity}(found,1)),...
%             nanmean(mappedX{knee_perplexity}(found,2))];
%         found_x = mappedX{knee_perplexity}(found,1);
%         found_y = mappedX{knee_perplexity}(found,2);
        
%         for f = 1:sum(found) % For each
%             plot([found_centroid(1),found_x(f)],...
%                 [found_centroid(2),found_y(f)],...
%                 'color',cmap(c,:));
%         end
%         
         scatter(mappedX{knee_perplexity}(found,1),...
        mappedX{knee_perplexity}(found,2),144,...
        'markerfacecolor',cmap(c,:),...
        'markeredgecolor',cmap(c,:)); hold on; 
        
%         scatter(found_centroid(1),found_centroid(2),360,'x',...
%             'markerfacecolor',cmap(c,:),'markeredgecolor',cmap(c,:),...
%             'linewidth',3);
        
      % Area Patches 
%     x = mappedX{knee_perplexity}(found,1); 
%     y = mappedX{knee_perplexity}(found,2); 
%     k = boundary(x,y,1); 
%     patch(x(k),y(k),cmap(c,:),'FaceAlpha',0.1); 
%     clear x y k 
    
   
    end
    
    %drawnow; pause(1);

    clear found_centroid found_x found_y 
    
end

% Coloring tsne By Pharmacological Cluster 
clear targets; 
targets = unique(Structural_Clusters); 
targets(isnan(targets)) = []; 

cmap = jet(size(targets,1)); %For each pharmacological Target  

figure; hold on
scatter(mappedX{knee_perplexity}(:,1),mappedX{knee_perplexity}(:,2),144,...
    'markerfacecolor',[0.4392 0.5020 0.5647],...
    'markeredgecolor',[0.4392 0.5020 0.5647]); 

for c = 1:size(targets,1) % For target  
    scatter(mappedX{knee_perplexity}(Structural_Clusters==c,1),...
        mappedX{knee_perplexity}(Structural_Clusters==c,2),144,...
        'markerfacecolor',cmap(c,:),...
        'markeredgecolor',cmap(c,:)); hold on; 
end 

%% Correlation - Hard Coding 
interest = 5813; % Use this to set the mutant/drug you're interested in 

num_return = 20; % Specify the number of compounds you'd like to see 

hom_mutants = [5757 5759 5777 5785 5805 5813]; 

numComp = 52; 
%% Methods Comparison 
% Raw data
cor_raw = corrcoef(giantdata'); 
[IX_raw,O_raw] = sort(cor_raw(interest,:)); 
sorted_names_raw = giantname(O_raw); 
sorted_targets_raw = Giant_Panel_Targets(O_raw); 
sorted_structures_raw = Structural_Clusters(O_raw); 
sorted_behavioural_clusters_raw = idx{numComp}(O_raw); 
sorted_Small_Panel_Logic = Small_Panel_Logic(O_raw); 

% PCA data 
% cor_score = corrcoef(score(:,1:knee_dim)'); 
% [IX_score,O_score] = sort(cor_score(interest,:));
% sorted_names_score = giantname(O_score); 
% sorted_targets_score = Giant_Panel_Targets(O_score); 
% sorted_structures_score = Structural_Clusters(O_score); 

% Weighted Correlation Data 
[IX_weight,O_weight] = sort(correlordergiantdata(interest,:));
sorted_names_weight = giantname(O_weight); 
sorted_targets_weight = Giant_Panel_Targets(O_weight); 
sorted_structures_weight = Structural_Clusters(O_weight); 

%% Figures

% Correlation Slopes 
a = 1; 
for m = hom_mutants % For each hom mutant
    [IX_hom_mutants_raw(a,:),O_hom_mutants_raw(a,:)] = sort(cor_raw(m,:)); 
    a = a + 1; 
end 

cmap = jet(size(hom_mutants,2));
figure; hold on
plot([0 size(giantdata,1)-1],[0 0],'--k','linewidth',1.5);
for m = 1:size(hom_mutants,2) % For each mutant 
    legend_lines(m) = plot(flip(IX_hom_mutants_raw(m,1:end-1)),'color',cmap(m,:),...
        'linewidth',3); 
end

% Nice Figure
axis([0 size(giantdata,1)-1 -1 1]);
[h,icons,plots,str] = legend(legend_lines,string(giantname(hom_mutants)),'location','southwest')
set(h,'Fontsize',17)
set(icons,'Linewidth',10)
set(icons(1:m),'Fontsize',17)
legend('boxoff')
box off; 
set(gca,'Fontsize',22);
xlabel('Compounds','Fontsize',34)
ylabel('Correlation','Fontsize',34)
set(gca, 'Layer','top')

%% Raw Data 
figure; 
for i = 1:num_return 
    subplot((num_return/5),5,i); hold on;  
    plot(giantdata(interest,:)*-1,'m','linewidth',1.5)
    plot(giantdata(interest,:),'r','linewidth',1.5)
    plot(giantdata(O_raw(i),:),'color',...
        [0.4392 0.5020 0.5647],'linewidth',1.5)
    %yl = ylim; 
    axis([1 26 ylim]) 
    try 
        title(string(sorted_names_raw{i}))
    catch 
    end 
end 

% Z-Scores 
% figure; 
% for i = 1:num_return 
%     subplot((num_return/5),5,i); hold on;  
%     plot(giantdata(interest,:)*-1,'m','linewidth',1.5)
%     plot(giantdata(interest,:),'r','linewidth',1.5)
%     plot(giantdata(O_score(i),:),'color',...
%         [0.4392 0.5020 0.5647],'linewidth',1.5)
%     %yl = ylim; 
%     axis([1 26 ylim]) 
%     try 
%         title(string(sorted_names_score{i}))
%     catch 
%     end 
% end 

% Weighted Analysis 
figure; 
for i = 1:num_return 
    subplot((num_return/5),5,i); hold on;  
    plot(giantdata(interest,:)*-1,'m','linewidth',1.5)
    plot(giantdata(interest,:),'r','linewidth',1.5)
    plot(giantdata(O_weight(i),:),'color',...
        [0.4392 0.5020 0.5647],'linewidth',1.5)
    %yl = ylim; 
    axis([1 26 ylim]) 
    try 
        title(string(sorted_names_weight{i}))
    catch 
    end 
end 

%% Representative Compounds

m = 6; % Set mutant 
compound = 5593; % Hard Code Compound
cmap = jet(size(hom_mutants,2)); % Generate colormap 
figure; hold on; clear scrap; clear legend_lines;
scrap = [giantdata(hom_mutants(m),:) giantdata(compound,:)]; 
for t = 1:2:24 
    rectangle('Position',[(t-0.5) min(scrap) 1 (max(scrap) - min(scrap))],...
        'FaceColor',[0.9608 0.9608 0.9608],'Edgecolor',[1 1 1]);
end 
rectangle('Position',[(26-0.5) min(scrap) 1 (max(scrap) - min(scrap))],...
        'FaceColor',[0.9608 0.9608 0.9608],'Edgecolor',[1 1 1])
plot([0 26.5],[0 0],'--k','linewidth',1.5);
legend_lines(1) = plot(giantdata(hom_mutants(m),:),'color',cmap(m,:),'linewidth',3);
legend_lines(2) = plot(giantdata(compound,:),'color',[0.4392 0.5020 0.5647],'linewidth',3);
axis([0.5 26.5 min(scrap) max(scrap)]); 

% Nice Figure
[h,icons,plots,str] = legend(legend_lines,[string(giantname(hom_mutants(m)))...
    string(giantname(compound))],'location','northwest');
set(h,'Fontsize',17)
set(icons,'Linewidth',10)
set(icons(1:2),'Fontsize',17)
legend('boxoff')
box off; 
set(gca,'Fontsize',32);
clear scrap; scrap = 2.5:4:27; scrap(end) = scrap(end) - 1; 
set(gca,'xtick',scrap)
set(gca,'xticklabels',{'Sleep','Sleep Bouts','Sleep Bout Length',...
    'Sleep Latency','Activity','Waking Activity','Total Activity'},...
    'Fontsize',20)
xlabel('Parameters','Fontsize',44)
ylabel('Z-Score','Fontsize',44)
set(gca, 'Layer','top')

%% Image Sc Versions  

figure; 
imagesc([giantdata(interest,:)*-1 ; giantdata(O_raw(1:num_return),:)])
colormap spring
colorbar

figure; 
imagesc([giantdata(interest,:)*-1 ; giantdata(O_score(1:num_return),:)])
colormap spring
colorbar

figure; 
imagesc([giantdata(interest,:)*-1 ; giantdata(O_weight(1:num_return),:)])
colormap spring
colorbar  

%% tSne 
figure; 
scatter(mappedX{knee_perplexity}(:,1),mappedX{knee_perplexity}(:,2),144,...
    'markerfacecolor',[0.4392 0.5020 0.5647],...
    'markeredgecolor',[0.4392 0.5020 0.5647]); hold on; 
scatter(mappedX{knee_perplexity}(interest,1),mappedX{knee_perplexity}(interest,2),144,...
        'markerfacecolor',cmap(1,:),...
    'markeredgecolor',cmap(1,:)); 

% Correlating 
% for i = 1:num_return
%     scatter(mappedX{knee_perplexity}(O_raw(end-i),1),mappedX{knee_perplexity}(O_raw(end-i),2),144,'g','filled')
% end 

% Anti-correlating
for i = 1:num_return
    scatter(mappedX{knee_perplexity}(O_raw(i),1),mappedX{knee_perplexity}(O_raw(i),2),144,'m','filled')
end 

% Coloring by Pharmacological Target 
figure; hold on; 
found = unique(sorted_targets_raw(1:num_return)); % Find unique targets
found = rmmissing(found); % Remove Missing targets 
cmap = jet(size(found,1)); % Assign colors to each target 
empty_cells = cellfun(@isempty,sorted_targets_raw(1:num_return)); 

scatter(mappedX{knee_perplexity}(:,1),mappedX{knee_perplexity}(:,2),...
    'markerfacecolor',[0.4392 0.5020 0.5647],...
    'markeredgecolor',[0.4392 0.5020 0.5647]); hold on; 
scatter(mappedX{knee_perplexity}(interest,1),mappedX{knee_perplexity}(interest,2),'g','filled')

for i = 1:num_return % For each pharm target
    if empty_cells(i) == 0 
    scatter(mappedX{knee_perplexity}(O_raw(i),1),mappedX{knee_perplexity}(O_raw(i),2),...
    'markerfacecolor',cmap(find(contains(string(found),string(sorted_targets_raw(i))) == 1,1,'first'),:),...
    'markeredgecolor',cmap(find(contains(string(found),string(sorted_targets_raw(i))) == 1,1,'first'),:));        
    end 
end 

% Coloring by Pharmacological Structure 
figure; hold on; 
found = unique(sorted_structures_raw(1:num_return)); % Find unique targets
found(isnan(found)) = []; % Remove Structural Clusters 
cmap = jet(size(found,1)); % Assign colors to each cluster 

scatter(mappedX{knee_perplexity}(:,1),mappedX{knee_perplexity}(:,2),...
    'markerfacecolor',[0.4392 0.5020 0.5647],...
    'markeredgecolor',[0.4392 0.5020 0.5647]); hold on; 
scatter(mappedX{knee_perplexity}(interest,1),mappedX{knee_perplexity}(interest,2),'g','filled')

for i = 1:num_return % For each pharm target
    if isnan(sorted_structures_raw(i)) == 0 % If it has a strucutral cluster assigned
    scatter(mappedX{knee_perplexity}(O_raw(i),1),mappedX{knee_perplexity}(O_raw(i),2),...
    'markerfacecolor',cmap(find(found == sorted_structures_raw(i),1,'first'),:),...
    'markeredgecolor',cmap(find(found == sorted_structures_raw(i),1,'first'),:));   
    end 
end 

% Coloring by correlation 
cmap = cool(size(IX_raw,2)); % Set colormap 
figure; hold on;
for i = 1:size(IX_raw,2) % For each compound  
    scatter(mappedX{knee_perplexity}(O_raw(i),1),mappedX{knee_perplexity}(O_raw(i),2),144,...
    'markerfacecolor',cmap(i,:),...
    'markeredgecolor',cmap(i,:)); 
end 

%% Weighted Correlation Analysis 

correlordergiantdata(1:size(giantdata,1),1:size(giantdata,1)) = NaN; 

for k = 1:size(giantdata,1) % For each compound
    for s = 1:size(giantdata,1) % For each comparison
        k % Report K
        s; % Report S
        x = 0; % Set variable
        y = 0; % Set variable
        r = 0; % Set variable
        xbase = 0; % Set variable
        ybase = 0; % Set variable
        
        for i = 1:size(giantdata,2) % For each parameter
            x = x + (Weights(i)*(giantdata(k,i)))^2;
            y = y + (Weights(i)*(giantdata(s,i)))^2;
        end
        
        xbase = sqrt(x/sum(Weights));
        ybase = sqrt(y/sum(Weights));
        
        for j = 1:size(giantdata,2) % For each parameter 
            r = r + Weights(j)^2*(giantdata(k,j)/xbase)*(giantdata(s,j)/ybase);
        end
        
        correlordergiantdata(k,s) = r/sum(Weights);
            
    end
end
