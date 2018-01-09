function pc = h2c(ph)
    [s1, s2] = size(ph);
    if s1==1
        pc = zeros(1,s2-1);
        for i = 1:s2-1
            pc(i) = ph(i)/ ph(s2);
        end
    elseif s2==1
        pc = zeros(s1-1,1);
        for i = 1:s1-1
            pc(i) = ph(i)/ ph(s1);
        end
    else
        error("not supported")
    end
end