function [ea, idx, idx_cts, ea_dist, ...
    ea_links, ea_idx, lifetimes, th] = ...
    gmm_sample_ea(X,reps,k_vals,a_size,s_size,rv,GMM_reps,max_its,method)

% Setup
options = statset('MaxIter',max_its); % Max number of GMM iterations
ea = zeros(a_size,a_size,'single'); % allocate ea matrix
[~,sample_a(:,1)] = datasample(X,a_size); % sample answers
idx_cts = []; % allocate

% Clustering & Evidence Accumulation 
for r = 1:reps % for each iteration
    
    clear sample k GMModels idx;
    
    sample = datasample(X,s_size); % sample points
    k = datasample(k_vals,1); % choose a value for k
    
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
ea = ea./reps; % scale (0-1)
ea_dist = 1 - ea; % invert
ea_dist(eye(size(ea_dist))==1) = 0; % remove diagonal values
ea_links = linkage(squareform(ea_dist),method); % compute linkage

% compute the maximum lifetime cut
lifetimes = diff(ea_links(:,3)); % diff linkage distances
[~, Ith] = max(lifetimes); % longest lifetime index

th = max(ea_links(Ith,3) + lifetimes(Ith)*0.5);

% apply the cut to the dendrogram
ea_idx = cluster(ea_links,'cutoff',th,'criterion','distance');

% Fit points to clusters 
kn = knnsearch(X(sample_a(:,1),:),X,'K',9);
idx = mode(reshape(ea_idx(kn(:)),size(kn)),2); 

end