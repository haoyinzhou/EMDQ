function ds = DQmult_dqxdt(dq,dt)
    if (size(dq,2) == 8)
%         ds0 = Qmult(dq(:,1:4),dr(:,1:4));
        ds1 = Qmult(dq(:,1:4),dt(:,5:8))+(dq(:,5:8));
        ds = [dq(:,1:4) , ds1];
    elseif (size(dq,2) == 4)
%         ds0 = Qmult(dq(:,1:2),dr(:,1:2));
        ds1 = Qmult(dq(:,1:2),dt(:,3:4))+(dq(:,3:4));
        ds = [dq(:,1:2), ds1];
    end
end
    

