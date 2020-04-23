function [ angle ] = dquat2angle( dq ) % works when angle is 'XYZ'
   N = size(dq,1);
    if (size(dq,2) == 8)
        angle = zeros(N,3);
        angle(:,1) = atan2(2*dq(:,3).*dq(:,1)-2*dq(:,2).*dq(:,4), 1-2*dq(:,3).*dq(:,3)-2*dq(:,4).*dq(:,4));
        angle(:,2) = asin(2*dq(:,2).*dq(:,3) + 2*dq(:,4).*dq(:,1));
        angle(:,3) = atan2(2*dq(:,2).*dq(:,1)-2*dq(:,3)*dq(:,4) , 1-2*dq(:,2).*dq(:,2)-2*dq(:,4).*dq(:,4));
    elseif (size(dq,2) == 4)
        angle = zeros(N,1);
      
    end

end

