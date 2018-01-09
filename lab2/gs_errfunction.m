function E = gs_errfunction(P0, Xobs)
    x  = reshape(Xobs(1:length(Xobs)/2),2,[]);
    xp = reshape(Xobs(length(Xobs)/2+1:end),2,[]); 
    xhat = reshape(P0(9+1:end),2,[]);
    H  = reshape(P0(1:9),3,3);
    
    x = [x ; ones(1,length(x))];
    xp = [xp ; ones(1,length(xp))];
    xhat = [xhat ; ones(1,length(xhat))];
    xhatp = H*xhat;
    
    E = zeros(2, size(x,2));
    for i=1:size(x,2)
        % compute the symmetric geometric error
        E(1,i) = pdist([reshape(h2c(x(:,i)),1,[]); reshape(h2c(xhat(:,i)),1,[])],'euclidean');
        E(2,i) = pdist([reshape(h2c(xp(:,i)),1,[]); reshape(h2c(xhatp(:,i)),1,[])],'euclidean');
    end
end