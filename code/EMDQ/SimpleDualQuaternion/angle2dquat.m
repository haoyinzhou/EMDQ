function DQ = angle2dquat(angle) % works when angle is 'XYZ'
    N = size(angle,1);
    if (size(angle,2) == 3)
        DQ = zeros(N,8);
%         c1 = cos(angle(:,1) * 0.5);
%         s1 = sin(angle(:,1) * 0.5);
%         c2 = cos(angle(:,2) * 0.5);
%         s2 = sin(angle(:,2) * 0.5);
%         c3 = cos(angle(:,3) * 0.5);
%         s3 = sin(angle(:,3) * 0.5);
%         DQ(:,1) = c1 .* c2 .* c3 - s1 .* s2 .* s3; 
%         DQ(:,2) = -s1 .* c2 .* c3 - c1 .* s2 .* s3; 
%         DQ(:,3) = -c1 .* s2 .* c3 + s1 .* c2 .* s3; 
%         DQ(:,4) = -s1 .* s2 .* c3 - c1 .* c2 .* s3;         
%         DQ(:,1:4) = repmat(2.0 * (-(DQ(:,1) < 0) + 0.5),[1,4]) .* DQ(:,1:4); % this aims to make DQ(:,1 > 0)
        t0 = cos(angle(:,3) * 0.5);
        t1 = sin(angle(:,3) * 0.5);
        t2 = cos(angle(:,2) * 0.5);
        t3 = sin(angle(:,2) * 0.5);
        t4 = cos(angle(:,1) * 0.5);
        t5 = sin(angle(:,1) * 0.5);
        DQ(:,1) = t0 .* t2 .* t4 + t1 .* t3 .* t5;
        DQ(:,2) = -t0 .* t3 .* t4 + t1 .* t2 .* t5;
        DQ(:,3) = -t0 .* t2 .* t5 - t1 .* t3 .* t4;
        DQ(:,4) = -t1 .* t2 .* t4 + t0 .* t3 .* t5;
    elseif (size(angle,2) == 1)
        DQ = zeros(N,4);
        DQ(:,1) = cos(angle(:,1) * 0.5);
        DQ(:,2) = -sin(angle(:,1) * 0.5);
    end
end
        




