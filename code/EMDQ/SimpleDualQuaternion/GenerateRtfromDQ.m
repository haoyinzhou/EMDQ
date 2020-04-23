function [ R, T ] = GenerateRtfromDQ(dq)

    N = size(dq,1);
    if (size(dq,2) == 8)
        R = dquat2rotMatrix(dq);
        QR = [dq(:,1:4), zeros(N,4)];
        QT2 = DQmult(dq, DQinv(QR));
        t2 = dquat2trans(QT2);
        T = t2 * R';
    elseif (size(dq,2)==4)
        R = dquat2rotMatrix(dq);
        QR = [dq(:,1:2), zeros(N,2)];
        QT2 = DQmult(dq, DQinv(QR));
        t2 = dquat2trans(QT2);
        T = t2 * R';
    end
end
