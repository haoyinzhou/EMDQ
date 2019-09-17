function dq = trans2dquat(T)
    N = size(T,1);
    if (size(T,2) == 3)
        dq = [ones(N,1), zeros(N,4), 0.5*T];
    elseif (size(T,2) == 2)
        dq = [ones(N,1), zeros(N,1), 0.5*T];
    end
end
        




