function dq = pos2dquat(v)
    N = size(v,1);
    if (size(v,2) == 3)
        dq = [ones(N,1), zeros(N,4), v];
    elseif (size(v,2) == 2)
        dq = [ones(N,1), zeros(N,1), v];
    end
end