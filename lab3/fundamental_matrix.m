function F_rank2 = fundamental_matrix(x1_test, x2_test)
    [p1, T1] = normalise2dpts(x1_test);
    [p2, T2] = normalise2dpts(x2_test);
    
%     W = [(p1(1,:).*p2(1,:))', (p1(1,:).*p2(2,:))', (p1(1,:))', ...
%          (p1(2,:).*p2(1,:))', (p1(2,:).*p2(2,:))', (p1(2,:))', ...
%          (p2(1,:))'         , (p2(2,:))'         , ones(8,1)];
    
    u1 = p1(1,:); 
    v1 = p1(2,:);
    u2 = p2(1,:);
    v2 = p2(2,:);
    W = [ (u1.*u2)', (v1.*u2)', u2', (u1.*v2)', (v1.*v2)', v2', u1', v1', ones(length(p1),1) ];

    
    [~,~,V] = svd(W);
    F_rank3 = reshape(V(:,size(V,2)),[3,3])';
    [U,D,V] = svd(F_rank3);
    D(3,3)=0;
    F_rank2 = U*D*V';
    F_rank2 = T2'*F_rank2*T1;
end
