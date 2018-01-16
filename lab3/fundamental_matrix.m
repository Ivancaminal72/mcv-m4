function F_es = fundamental_matrix(x1_test, x2_test)
    p1 = x1_test./x1_test(3,:);
    p2 = x2_test./x2_test(3,:);
    
    W = [(p1(1,:).*p2(1,:))', (p1(1,:).*p2(2,:))', (p1(1,:))', ...
         (p1(2,:).*p2(1,:))', (p1(2,:).*p2(2,:))', (p1(2,:))', ...
         (p2(1,:))'         , (p2(2,:))'         , ones(8,1)];
    
    [~,~,V] = svd(W);
    F = reshape(V(:,size(V,2)),[3,3]);
    [U,D,V] = svd(F);
    D(3,3)=0;
    F_es = U*D*V';
end