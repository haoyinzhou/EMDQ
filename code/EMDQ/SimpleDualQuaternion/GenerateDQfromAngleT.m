function [ outDualQ ] = GenerateDQfromAngleT( angle, T)
    QR = angle2dquat(angle);
    QT = trans2dquat(T);
    outDualQ = DQmult(QR,QT);
end
