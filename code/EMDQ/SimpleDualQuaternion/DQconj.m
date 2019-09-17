function dq = DQconj(dv)
    if (size(dv,2) == 8)
        dq = [dv(:,1),-dv(:,2:5),dv(:,6:end)];
    elseif (size(dv,2) == 4)
        dq = [dv(:,1),-dv(:,2),dv(:,3:end)];
    end
end
