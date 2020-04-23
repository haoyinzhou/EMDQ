function T = dquat2trans(dq)
    N = size(dq,1);
    if (size(dq,2) == 8)
        T = 2.0 * dq(:,6:8);
    elseif (size(dq,2) == 4)
        T = 2.0 * dq(:,3:4);
    end
end
        




