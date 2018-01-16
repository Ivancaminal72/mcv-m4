function F_es = fundamental_matrix(x1_test, x2_test)
    [p1, T1] = normalise2dpts(x1_test);
    [p2, T2] = normalise2dpts(x2_test);
    
    W = [(p1(1,:).*p2(1,:))', (p1(1,:).*p2(2,:))', (p1(1,:))', ...
         (p1(2,:).*p2(1,:))', (p1(2,:).*p2(2,:))', (p1(2,:))', ...
         (p2(1,:))'         , (p2(2,:))'         , ones(8,1)];
    
    [~,~,V] = svd(W);
    F = reshape(V(:,size(V,2)),[3,3])';
    [U,D,V] = svd(F);
    D(3,3)=0;
    F_es = U*D*V';
    F_es = inv(T2)*F_es*T1;
end
