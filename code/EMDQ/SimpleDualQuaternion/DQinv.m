function [ dq_inv ] = DQinv( dq )
    N = size(dq,1);
    if (size(dq,2) == 8)
%         dq_inv = zeros(N,8);
%         normdq = DQnorm(dq);
%         normdq_2 = zeros(N,2);
%         normdq_2(:,1) = normdq(:,1) .* normdq(:,1);
%         normdq_2(:,2) = 2 * normdq(:,1) * normdq(:,2);
%         dq_conj = [dq(:,1),-dq(:,2:4),dq(:,5),-dq(6:8)];        
%         for l=1:8
%             if l <= 4
%                 temp = [dq_conj(:,l), zeros(N,1)];
%                 res = Ddiv(temp,normdq_2);
%                 dq_inv(:,l) = dq_inv(:,l)+res(:,1);
%                 dq_inv(:,l+4) = dq_inv(:,l+4)+res(:,2);
%             else
%                 temp = [zeros(N,1), dq_conj(:,l)];
%                 res = Ddiv(temp,normdq_2);
%                 dq_inv(:,l) = dq_inv(:,l)+res(:,2);
%                 dq_inv(:,l-4) = dq_inv(:,l-4)+res(:,1);
%             end
%         end
        dq_inv = [dq(:,1) -dq(:,2:4) dq(:,5) -dq(:,6:8)];
        
    elseif(size(dq,2) == 4)
        dq_inv = [dq(:,1) -dq(:,2:4)];

    end

end


