function [ea, idx, idx_cts, ea_dist, ...
    ea_links, ea_idx, lifetimes, th, sample_a,sample_a_n] = ...
    gmm_sample_ea(X,reps,k_vals,a_size,s_vals,rv,GMM_reps,max_its,method,nn)

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
        datasample(sample_a(ea_idx == c,1),min(histcounts(ea_idx)),...
        'Replace',false)];
    % sample indicies refering to X,
    % sample as many points as the smallest cluster 
end 
sample_a_n(:,2) = reshape(repmat(1:max(ea_idx),min(histcounts(ea_idx)),1),[],1); 
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