function H = homography2d(x1, x2)
    [x1, T1] = DLT_normalization(x1);
    [x2, T2] = DLT_normalization(x2);
    
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
    H=inv(T2')*reshape(V(:,size(V,2)),[3,3])*T1;
end

function [ x_normalized, T ] = DLT_normalization(x)
    x(1,1:4)=x(1,1:4)./x(3,1:4);
    x(2,1:4)=x(2,1:4)./x(3,1:4);
    x(3,1:4)=1;
    
    x_mean = mean(x(1:2,1:4));
    x_new(1,1:4) =  x(1,1:4)-x_mean(1);
    x_new(2,1:4) = x(2,1:4)-x_mean(2);
    dist = sqrt(x_new(1,1:4).^2 + x_new(2,1:4).^2);
    meandist = mean(dist(:).^2);  
    s = sqrt(2)/meandist;
    T = [s   0   -s*x_mean(1)
         0   s   -s*x_mean(2)
         0   0      1      ];
    x_normalized = T*x;
end
