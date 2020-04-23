function [ R ] = angle2rotMatrix( angle )
    N = size(angle,1);
    if (size(angle,2) == 3)
        R = zeros(3,3,N);
        sx = sin(angle(:,3));
        sy = sin(angle(:,1));
        sz = sin(angle(:,2));
        cx = cos(angle(:,3));
        cy = cos(angle(:,1));
        cz = cos(angle(:,2));
        R(1,1,:) = cy.*cx;
        R(1,2,:) = cy.*sx;
        R(1,3,:) = -sy;
        R(2,1,:) = sz.*sy.*cx - cz.*sx;
        R(2,2,:) = sz.*sy.*sx + cz.*cx;
        R(2,3,:) = sz.*cy;
        R(3,1,:) = cz.*sy.*cx + sz.*sx;
        R(3,2,:) = cz.*sy.*sx - sz.*cx;
        R(3,3,:) = cz.*cy;        
    elseif (size(angle,2) == 1)
        R = zeros(2,2,N);
        R(1,1,:) = cos(angle(:));
        R(1,2,:) = sin(angle(:));
        R(2,1,:) = -R(1,2,:);
        R(2,2,:) = R(1,1,:);
    end


end
