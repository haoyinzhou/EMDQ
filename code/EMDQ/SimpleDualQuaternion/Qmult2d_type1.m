function s = Qmult2d_type1(q,r) % q = [x 0 0 x], r = [x 0 0 x]
    if (size(q,2)==2)
        s0 = q(:,1).*r(:,1) - q(:,2).*r(:,2);
        s3 = q(:,2).*r(:,1) + q(:,1).*r(:,2);

        s = [s0 s3];
    end
    
end

