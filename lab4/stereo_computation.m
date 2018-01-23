function disparity = stereo_computation(imL, imR, dmin, dmax, win_size, match_cost)
    imL = double(imL);
    imR = double(imR);
    
    if(mod((win_size-1),2) == 0)
        padding = (win_size-1)/2;
    else
        padding = win_size/2;
        win_size = win_size +1;
    end
    disparity = zeros(size(imL));
    [height,width] = size(imL);
    if(match_cost == 'BW')%Define bilateral weigths
        WL = zeros(size(win_size));
        for i = 1:length(win_size)
            for j = 1:length(win_size)                
                WL(i,j) = exp(-(pdist([[i,j];[padding,padding]],'euclidean'))/win_size);
            end
        end
        WR = WL;
    else %Define standard weights        
        W = (1/(win_size^2))*ones(win_size,win_size);
    end
    
    imL = padarray(imL, [padding, padding]);
    imR = padarray(imR, [padding, padding]);
    for i = 1:height
        for j = 1:width
            if((j+dmax) < width)
                x_max = j+dmax; 
            else
                x_max = width; 
            end
            if((j-dmax) > 1)
                x_min = j-dmax;
            else
                x_min = 1;
            end
            best_x=Inf;                        
            if match_cost == "SSD"
                ssd_min=Inf;
                for x = x_min:x_max
                    patchL=imL(i:i+win_size-1,j:j+win_size-1);
                    patchR=imR(i:i+win_size-1,x:x+win_size-1);
                    ssd=sum(sum(W.*((patchL-patchR).^2)));
                    if(ssd < ssd_min)
                        best_x = x;
                        ssd_min = ssd;
                    end
                end
            elseif match_cost == "NCC"         
                ncc_max=-Inf;
                for x = x_min:x_max
                    patchL=imL(i:i+win_size-1,j:j+win_size-1);
                    patchR=imR(i:i+win_size-1,x:x+win_size-1);
                    L_norm = sum(sum(W.*patchL));
                    R_norm = sum(sum(W.*patchR));
                    patchL_norm =(patchL-L_norm);
                    patchR_norm =(patchR-R_norm);
                    sigmaL = sqrt(sum(sum(W.*patchL_norm.^2)));
                    sigmaR = sqrt(sum(sum(W.*patchR_norm.^2)));
                    ncc=sum(sum(W.*patchL_norm.*patchR_norm))/(sigmaL*sigmaR);                    
                    if(ncc > ncc_max)
                        best_x = x;
                        ncc_max = ncc;
                    end
                end
            elseif match_cost == "BW"   
                E_min=Inf;
                for x = x_min:x_max
                    patchL=imL(i:i+win_size-1,j:j+win_size-1);
                    patchR=imR(i:i+win_size-1,x:x+win_size-1);           
                    e=sum(sum(abs(patchL-patchR)));
                    E=sum(sum(WL.*WR.*e))/sum(sum(WL.*WR));
                    if(E < E_min)
                        best_x = x;
                        E_min = E;
                    end
                end
            end
            disparity(i,j)= max(dmin,abs(best_x-j))/dmax;
        end
    end
end