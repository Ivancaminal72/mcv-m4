function H = homography2d(x1, x2)
    A = zeros(8,9);
    Ai = zeros(2,9);
    
    for i=1:size(x1,2)
        Ai(1,4:6) = -(x1(3,i)*(x2(:,i)));
        Ai(1,7:9) = (x1(2,i)*(x2(:,i)));
        Ai(2,1:3) = (x1(3,i)*(x2(:,i)));
        Ai(2,7:9) = -(x1(1,i)*(x2(:,i)));
        
        A(i:i+1,:) = Ai;
    end
    [~,~,V] = svd(A);
    H=reshape(V(:,size(V,2)),[3,3]);
end