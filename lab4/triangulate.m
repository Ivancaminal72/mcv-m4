function X = triangulate(x1, x2, P1, P2, imsize)
    x = homog(x1);
    xp = homog(x2);
    
    H = [2/imsize(1), 0          , -1;
         0          , 2/imsize(2), -1;
         0          , 0          , 1];
     
     x = H*x;
     xp = H*xp;
     P1 = H*P1;
     P2 = H*P2;
     
     A = [x(1)*P1(3,:) - P1(1,:);
          x(2)'*P1(3,:) - P1(2,:);
          xp(1)*P2(3, :) - P2(1,:);
          xp(2)*P2(3,:) - P2(2,:)];
      
    [~,~, V] = svd(A);
    
    X= V(:,size(V,1));
    X = X ./ X(end);
          
end