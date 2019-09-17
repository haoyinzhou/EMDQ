function ds = DQmult(dq,dr)
    if (size(dq,2) == 8)
        ds0 = Qmult(dq(:,1:4),dr(:,1:4));
        ds1 = Qmult(dq(:,1:4),dr(:,5:8))+Qmult(dq(:,5:8),dr(:,1:4));
        
%         [Qmult(dq(:,1:4),dr(:,1:4));
%         Qmult(dq(:,1:4),dr(:,5:8));
%         Qmult(dq(:,5:8),dr(:,1:4))]        
        ds = [ds0 , ds1];
    elseif (size(dq,2) == 4)
        ds0 = Qmult2d_type1(dq(:,1:2),dr(:,1:2));
        ds1 = Qmult2d_type2(dq(:,1:2),dr(:,3:4))+Qmult2d_type3(dq(:,3:4),dr(:,1:2));
%         
%         [Qmult(dq(:,1:2),dr(:,1:2));
%          Qmult(dq(:,1:2),dr(:,3:4));
%          Qmult(dq(:,3:4),dr(:,1:2))]
        
        ds = [ds0, ds1];
    end
end
    

