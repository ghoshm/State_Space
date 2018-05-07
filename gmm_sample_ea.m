function [ea, idx, idx_cts, ea_dist, ...
    ea_links, ea_idx, lifetimes, th, sample_a,sample_a_n] = ...
    gmm_sample_ea(X,reps,k_vals,a_size,s_vals,rv,GMM_reps,max_its,method,nn)

%% Info 
    % marcus.ghosh.11@ucl.ac.uk 
    
% Details 
    % Clustering using Gaussian Mixture Models, as both K, 
    % and the data the model is fit to is varied over multiple iterations. 
    % Each iteration evidence (pairwise co-occurances in the same cluster) 
    % is accumulated on a set of probe points. After this process the evidence 
    % accumulation matrix is hierarchically clustered and the number of
    % clusters is determined by a maximum lifetime cut. Finally clusters
    % are normalised for size and all data is assigned using KNN. 
    
    % This process is suited to datasets where there are too many points 
    % to cluster in a single pass 
    
    % Based upon (Fred & Jain 2005 - 10.1109/ICPR.2002.1047450) 
    
% Steps 
    % Sample (a_size) probe points for evidence accumulation 
    % Unifromly sample values of K and sample sizes from set ranges (k_vals and s_vals) 
    % Iteratively (for reps iterations)  
        % Fit a GMM to the sampled data with k components 
        % Assign each probe point to the component with the highest
            % posterior probability for that point 
        % Accumulate evidence (pairwise co-occurances in the same cluster) 
            % on the probe points 
    % Hierarchically cluster the evidence accumulation matrix and determine
        % the final number of clusters using a maximum lifetime cut 
    % Normalise the clusters for size by sampling an equal number of points 
        % from each cluster   
    % Assign all data points to final size normalised clusters using KNN 
    
% Dependencies 
  % Statistics and Machine Learning Toolbox

% Inputs 
    % X - input data; n (points) x d (dimensions/parameters) matrix 
    % reps - number of iterations, for sampling/clustering/evidence accumulation   
    % k_vals - range of component numbers to fit with GMM; vector 
    % a_size - number of probe points to sample for evidence accumulation
    % s_vals - range of sample sizes to use; [min,max] 
    % rv - regularization values for GMM 
    % GMM_reps - number of reps for each GMM model 
    % max_its - maximum number of iterations for each GMM model 
    % method - linkage method for hierarchical clustering (eg. 'average')  
    % nn - number of nearest neighbours for KNN assignment 
    
% Outputs 
    % ea - evidence accumulation matrix; a_size x a_size 
    % idx - final cluster assignments; n x 1 
    % idx_cts - cluster centroids from every iteration; centroids x d 
    % ea_dist - evidence accumulation distance matrix; a_size x a_size 
    % ea_links - hierarchical cluster tree; (a_size - 1) x 3
    % ea_idx - probe points cluster assignments; a_size x 1 
    % lifetimes - hierarchical partition lifetimes; (a_size - 2) x 3 
    % th - maximum lifetime threshold
    % sample_a - indicies in X of probe points sampled; a_size x 1 
    % sample_a_n - indicies of points in size normalised clusters, used for
        % KNN; (smallest cluster size x number of clusters) x 2 (1 = index in X,
            % 2 = final cluster assignment) 
    
% % Example
% % Example Data
% X = [normrnd(1,1,100000,2) ; normrnd(10,1,100000,2)];
% 
% % Settings
% reps = 100; % set the number of repetitions
% k_vals = 2:10; % set values of k (clusters) to try
% a_size = 1000; % number of probe points
% s_vals = [100,1000]; % min & max points to sample (uniformly)
% rv = 0.0001; % regularization value for GMM
% GMM_reps = 5; % number of GMM Models to fit per iteration
% max_its = 1000; % number of GMM iterations (per repetition)
% method = 'average'; % linkage measure
% nn = 50; % number of nearest neighbours
% 
% % Cluster 
% [ea, idx, idx_cts, ea_dist, ...
%     ea_links, ea_idx, lifetimes, th, sample_a,sample_a_n] = ...
%     gmm_sample_ea(X,reps,k_vals,a_size,s_vals,rv,GMM_reps,max_its,method,nn);
% 
% % Scatter Figure 
% figure; box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
% set(gca,'FontName','Calibri');
% 
% gscatter(X(:,1),X(:,2),idx); % scatter points 
% 
% % Evidence Accumulation Figure  
% figure; box off; set(gca, 'Layer','top'); set(gca,'Fontsize',32); % Format
% set(gca,'FontName','Calibri');
% 
% % Dendrogram
% subplot('position',[0.0500    0.8178    0.9000    0.1322]);
% [H,~,perm] = dendrogram(ea_links,size(ea,1),'colorthreshold',th); % dendrogram
% axis off;
% 
% % Maximum lifetime threshold line
% line(get(gca,'xlim'), [th th],'LineStyle',':','LineWidth',1.5);
% text(1,double(th),'Maximum Lifetime Cut','verticalalignment','bottom',...
%     'FontName','Calibri','FontSize',16);
% 
% % EA Matrix
% subplot('position', [0.0500    0.0500    0.9000    0.7478]);
% imagesc(ea(perm,perm));
% axis off
% c = colorbar;
% c.Label.String = 'E.A. Index';
% c.Location = 'southoutside';
% c.FontSize = 16;

%% Function 

% Setup
options = statset('MaxIter',max_its); % Max number of GMM iterations
ea = zeros(a_size,a_size,'single'); % allocate ea matrix
[~,sample_a(:,1)] = datasample(X,a_size,'Replace',false); % sample answers (indicies in X) 
idx_cts = []; % allocate
s_sizes(:,1) = datasample(s_vals(1):s_vals(2),reps); % uniform sampling of size samples  
ks(:,1) = datasample(k_vals,reps); % uniform sampling of k values  

% Clustering & Evidence Accumulation 
for r = 1:reps % for each iteration
    
    clear sample k GMModels idx;
    
    sample = datasample(X,s_sizes(r,1),'Replace',false); % sample points
    k = ks(r,1); % choose a value for k
    
    GMModels = fitgmdist(sample,k,...
        'Options',options,'RegularizationValue',...
        rv,'Replicates',GMM_reps); % Fit k gaussians to sample, GMM_reps times
    
    idx = cluster(GMModels,X(sample_a(:,1),:)); % cluster the answer points
    
    idx_cts = [idx_cts ; grpstats(X(sample_a(:,1),:),idx,'mean')]; % answer centroids
    
    ea = ea + ((repmat(idx,1,size(idx,1)) - repmat(idx',size(idx,1),1)) ==0); % add to the EA Matrix
    
    % Report Progress 
    disp(horzcat('Iteration ',num2str(r),' of ',num2str(reps)));
    
end

% construct the dendrogram
disp('Constructing dendrogram'); % report progress 
ea = ea./reps; % scale (0-1)
ea_dist = 1 - ea; % invert
ea_dist(eye(size(ea_dist))==1) = 0; % make sure diagonal values are zero 
ea_links = linkage(squareform(ea_dist),method); % compute linkage

% compute the maximum lifetime cut
lifetimes = diff(ea_links(:,3)); % diff linkage distances
[~, Ith] = max(lifetimes); % longest lifetime index

th = max(ea_links(Ith,3) + lifetimes(Ith)*0.5); % lifetime threshold 

% apply the cut to the dendrogram
ea_idx = cluster(ea_links,'cutoff',th,'criterion','distance');

% Normalise clusters for size 
sample_a_n = []; % allocate 
for c = 1:max(ea_idx) % for each cluster 
    sample_a_n = [sample_a_n ; ...
        datasample(sample_a(ea_idx == c,1),min(histcounts(ea_idx,...
        'binmethod','integer')),'Replace',false)];
    % sample indicies refering to X,
    % sample as many points as the smallest cluster 
end 
sample_a_n(:,2) = reshape(repmat(1:max(ea_idx),min(histcounts(ea_idx,...
    'binmethod','integer')),1),[],1); 
ea_idx_n = sample_a_n(:,2); % normalised cluster indicies   

% Fit points to clusters
disp('Fitting data into clusters'); % report progress
try
    kn = knnsearch(X(sample_a_n(:,1),:),X,'K',nn); % find nn nearest neighbours
    idx = mode(reshape(ea_idx_n(kn(:)),size(kn)),2); % assign to mode neighbour cluster

catch % if the X by nn array is too large for memory 
    disp('Caught large matrix, chunking data'); % report mode 
    idx = zeros(size(X,1),1,'single'); % allocate 
    cks = [1:round(size(X,1)/1000):size(X,1) (size(X,1)+1)]; % break into 1000 chunks 
    for i = 1:(length(cks) - 1) % for each chunk
        kn = knnsearch(X(sample_a_n(:,1),:),X(cks(i):(cks(i+1)-1),:),...
            'K',nn); % find nn nearest neighbours
        idx(cks(i):(cks(i+1)-1),1) = ...
            mode(reshape(ea_idx_n(kn(:)),size(kn)),2); % assign to mode neighbour cluster
    end 
end

end