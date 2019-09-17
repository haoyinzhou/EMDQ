function [ outDualQ ] = GenerateDQfromRt( R, T)

    if (size(R,1) == 3)
        QR = rotMatrix2dquat(R);
%         QT = trans2dquat(T);
        t2 = T * R;
        QT = trans2dquat(t2);
        outDualQ = DQmult(QT,QR);
    elseif (size(R,1)==2)
        QR = rotMatrix2dquat(R);
%         QT = trans2dquat(T);
%         outDualQ = [QR(:,1:2) QT(:,3:4)];
        t2 = T * R;
        QT = trans2dquat(t2);
        outDualQ = DQmult(QT,QR);
    end
end
