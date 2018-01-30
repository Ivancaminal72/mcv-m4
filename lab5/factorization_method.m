function [Pproj, Xproj] = factorization_method(x1, x2, initialization, th1, th2)
    [x1,T1]=normalise2dpts(x1);
    [x2,T2]=normalise2dpts(x2);
    M = [x1; x2];

    
    if (initialization == "ones")
        lambda = ones(size(M));
    elseif (initialization == "strum")
        lambda1 = ones(3,size(M,2));
        F21 = fundamental_matrix(x2, x1);
        [~, ~, V] = svd(F21);
        e = V(:,length(V))./V(length(V),length(V));
        lambda2 = lambda1;
        for j = 1:length(x2)
            lambda2(:,j) =  (x1(:,j)'*F21*(cross(e,x2(:,j))/norm(cross(e,x2(:,j)))^2))*lambda1(:,j);
        end
        lambda = [lambda1;lambda2];
    else
        error("Not valid lambda initialization");
    end
    
    dx_old = Inf;
    while true
        new_lambda = lambda;
        d_old = Inf;
        d = 1;
        cols = true;
        while (abs(d - d_old)/d) < th1
            if cols
                for i=1:size(lambda,2)
                    new_lambda(:,i) = lambda(:,i)./norm(lambda(:,i));
                end            
            else
                for i=1:size(lambda,1)
                    new_lambda(i,:) = lambda(i,:)./norm(lambda(i,:));
                end
            end
            d_old = d;
            d = sum(sum(sqrt((new_lambda-lambda).^2))); 
            cols = ~cols;
        end

        M = M.*lambda;
        [U, D, V] = svd(M);        
        Pproj=U*D(:,1:4);
        V = V';
        Xproj=V(1:4,:);        
        x12 = [x1; x2];
        dx = sum(sum(sqrt((x12-Pproj*Xproj).^2)));        
        disp((dx-dx_old)/dx);
        
        if (abs(dx - dx_old)/dx) < th2
            break;
        else            
            xp = Pproj*Xproj;
            lambda(1,:) = xp(3,:);
            lambda(2,:) = xp(6,:);
            dx_old = dx;
        end
    end
    Pproj(1:3,:) = inv(T1)*Pproj(1:3,:);
    Pproj(4:6,:) = inv(T2)*Pproj(4:6,:);    
end