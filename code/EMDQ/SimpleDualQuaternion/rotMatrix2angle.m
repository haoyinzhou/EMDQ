function angle = rotMatrix2angle(R)
    N = size(R,3);       
    if (size(R,1) == 3)
        angle = zeros(N,3);    
        angle(:,1) = asin(-R(1,3,:));
        angle(:,2) = atan2(R(2,3,:), R(3,3,:));
        angle(:,3) = atan2(R(1,2,:), R(1,1,:));
    elseif (size(R,1) == 2)
        angle = atan2(R(1,2,:), R(1,1,:));
    end
end
