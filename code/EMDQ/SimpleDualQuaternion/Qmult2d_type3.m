function s = Qmult2d_type3(q,r)  % q = [0 x x 0], r = [x 0 0 x]
    if (size(q,2)==2)
        s1 = q(:,1).*r(:,1) - q(:,2).*r(:,2);
        s2 = q(:,2).*r(:,1) + q(:,1).*r(:,2);
       
        s = [s1 s2];
    end
    
end

