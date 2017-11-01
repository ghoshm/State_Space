function [GMModels,AIC,BIC,idx,nlogl,P] = GhoshMarcusModels(data,num_clusters)

options = statset('MaxIter',100); % Hard coded number of iterations 

AIC = zeros(1,num_clusters); 
BIC = zeros(1,num_clusters); 
GMModels = cell(num_clusters,1);
idx = cell(num_clusters,1);
nlogl = cell(num_clusters,1);
P = cell(num_clusters,1);
score_values = unique(data); % Find unique scores 
score_zero = knnsearch(score_values,0); % Find the closest to zero
rv = score_values(score_zero); % Regularization value 

parfor k = 1:num_clusters % Try this many groups
    
    GMModels{k} = fitgmdist(data,k,...
        'Options',options,'RegularizationValue',...
        abs(rv),'Replicates',5); % Fit K gaussians
    
    AIC(k)= GMModels{k}.AIC; % Extract AIC
    BIC(k)= GMModels{k}.BIC; % Extract BIC
    
    % Cluster using this mixing
    [idx{k},nlogl{k},P{k}] = cluster(GMModels{k},data);
end

end
 