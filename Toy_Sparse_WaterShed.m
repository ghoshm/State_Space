% Toy Sparse WaterShed Algorithm 

% Based upon 
%https://uk.mathworks.com/help/stats/cluster-data-from-mixture-of-gaussian-distributions.html

% Generate two gaussians 
    % + a 3rd = 1st, slightly offset 
rng default;  % For reproducibility
mu1 = [1 2];
sigma1 = [3 .2; .2 2];
mu2 = [-1 -2];
sigma2 = [2 0; 0 1];
X = [mvnrnd(mu1,sigma1,200); mvnrnd(mu2,sigma2,100) ;...
    mvnrnd(mu2-.5,sigma2,200)];
n = size(X,1);

figure;
scatter(X(:,1),X(:,2),10,'ko')

% Fit a 3 component GMM 
options = statset('Display','final');
gmm = fitgmdist(X,14,'Options',options);

% Cluster 
idx = gmm.cluster(X); 

scatterhist(X(:,1),X(:,2),'Group',idx,'Kernel','on')

% Sparse WaterShedding 
results=struct();
options=optimoptions('fminunc','Algorithm','quasi-newton','ObjectiveLimit',-1e100,'HessUpdate','steepdesc','Display','final-detailed',...
    'MaxFunEvals',inf,'MaxIter',5000);
for iCluster=1:size(gmm.mu,1)
    % Grab this cluster's mean coords
    x0=gmm.mu(iCluster,:);
    
    % Run gradient descent on each data point, we invert the PDF to find its maximum and scale by 1e-30 to bring the PDFs values into
    % a more typical range for fminunc
%     [results(iCluster).x,results(iCluster).fval,results(iCluster).exitflag,results(iCluster).output,results(iCluster).grad]=...
%         fminunc(@(x)-gmm.pdf(x)*1e-30,x0,options);
    [results(iCluster).x,results(iCluster).fval,results(iCluster).exitflag,results(iCluster).output,results(iCluster).grad]=...
        fminunc(@(x)-gmm.pdf(x),x0,options);
end

%fminunc outputs: 
    % x - the local minimum 
    % fval - the value of the objective function @ the solution 
    % Exitflag - exit condition  
    % Output - information about the process 
    % Grad - gradient of the solution @ x 

% Options: 
    % Quasi-Newton - computationally easy optimization algorithm 
    % Objective-Limit - a stopping criteria 
    % Hess Update (Steepdesc) - 
    % Display 
    % MaxFunEvals 
    % MaxIter 
    
% Notes: 
    % Optimization Algorithms 
        % Noting that the search for a minimum or maximum of a scalar-valued
        % function is nothing else than the search for the zeroes of the 
        % gradient of that function
    
    
   options=optimoptions('fminunc','Algorithm','quasi-newton','ObjectiveLimit',-1e100,'HessUpdate','steepdesc','Display','final-detailed',...
    'MaxFunEvals',inf,'MaxIter',5000); 
% Figure: 
figure; hold on; axis([-3 5  -1 1.5]); 
for i = 1:14
scatter(gmm.mu(i,1),gmm.mu(i,2),'k','filled');
scatter(results(i).x(1),results(i).x(2),'r','filled');
plot([results(i).x(1),gmm.mu(i,1)],[results(i).x(2),gmm.mu(i,2)],...
    '--k','linewidth',1.5); 
pause(5); 
end 

% Full WaterShed 
[~,~,~,~,clusters]=wshedProcess(embeddingValues,desiredK); 

