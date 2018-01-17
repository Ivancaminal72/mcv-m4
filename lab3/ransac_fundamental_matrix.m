function [F, idx_inliers] = ransac_fundamental_matrix(p1, p2, th)
[~, Npoints] = size(p1);

% ransac
it = 0;
max_it = 1000;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
while it < max_it
    
    points = randomsample(Npoints, 8);
    F = fundamental_matrix(p1(:,points), p2(:,points));
    inliers = compute_inliers(F, p1, p2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^8;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    max_it = min(max_it, log(1-p)/log(pNoOutliers));
    
    it = it + 1;
end

% compute F from all the inliers
F = fundamental_matrix(p1(:,best_inliers), p2(:,best_inliers));
idx_inliers = best_inliers;


function idx_inliers = compute_inliers(F, p1, p2, th)
    % transformed points (in both directions)
    Fp1 = F * p1;
    Ftp2 = F' * p2;
    
    % compute the numerator of the Sampson error
    p2tFp1= diag(p2'*Fp1)';
    
    % compute the Sampson error
    d2 = p2tFp1.^2 ./(Fp1(1,:).^2 + Fp1(2,:).^2 ... 
                + Ftp2(1,:).^2 + Ftp2(2,:).^2);
    
    idx_inliers = find(d2 < th.^2);
    
    function idx_inliers = compute_inlierss(F, x1, x2, th)
    x1=normalise(x1);
    x2=normalise(x2);
    x2tFx1 = zeros(1,length(x1));
	for n = 1:length(x1)
	    x2tFx1(n) = x2(:,n)'*F*x1(:,n);
    end
    Fx1 = F*x1;
    Ftx2 = F'*x2;
    
    d2 =  x2tFx1.^2 ./ (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);
    
    idx_inliers = find(d2 < th.^2);
    


function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat