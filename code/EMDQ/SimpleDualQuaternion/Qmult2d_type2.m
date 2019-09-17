function s = Qmult2d_type2(q,r) % q = [x 0 0 x], r = [0 x x 0]
    if (size(q,2)==2)
        s1 = q(:,1).*r(:,1) + q(:,2).*r(:,2);
        s2 = q(:,1).*r(:,2) - q(:,2).*r(:,1);
      
        s = [s1 s2];
    end
    
end

