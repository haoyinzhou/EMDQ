function [ normDQ ] = DQnorm( dq )
    N = size(dq,1);
    if (size(dq,2) == 8)
        normDQ = zeros(N, 2);
        normDQ(:,1) = sqrt(dq(:,1).^2 + dq(:,2).^2 + dq(:,3).^2 + dq(:,4).^2); 
        q0_conj = [dq(:,1), -dq(:,2), -dq(:,3), -dq(:,4)];
        q1_conj = [dq(:,5), -dq(:,6), -dq(:,7), -dq(:,8)];
       
        temp1 = Qmult(q0_conj, dq(:,5:8));
        temp2 = Qmult(q1_conj, dq(:,1:4)); 
        normDQ(:,2) = (temp1(:,1) + temp2(:,1)) ./ (2 * normDQ(:,1));
    elseif(size(dq,2) == 4)

    end

end