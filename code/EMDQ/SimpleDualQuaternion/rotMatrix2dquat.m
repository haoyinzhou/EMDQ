function dq = rotMatrix2dquat(R)
    N = size(R,3);       
    if (size(R,1) == 3)
        dq = zeros(N,8);    
        dq(:,1) = 0.5 * sqrt(1 + R(1,1,:) + R(2,2,:) + R(3,3,:));
        temp = 0.25 / dq(:,1);
        dq(:,2) = temp .* ( R(3,2,:) - R(2,3,:));
        dq(:,3) = temp .* ( R(1,3,:) - R(3,1,:));
        dq(:,4) = temp .* ( R(2,1,:) - R(1,2,:));
    elseif (size(R,1) == 2)
%         dq = zeros(N,4);    
%         dq(:,1) = R(1,1,:);
%         dq(:,2) = R(2,1,:);
        dq = zeros(N,4);    
        dq(:,1) = 0.5 * sqrt(2.0 + R(1,1,:) + R(2,2,:));
        temp = 0.25 / dq(:,1);
        dq(:,2) = temp .* ( R(2,1,:) - R(1,2,:));
end
