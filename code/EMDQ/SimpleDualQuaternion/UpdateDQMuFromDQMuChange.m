function [ dq_new, mu_new ] = UpdateDQMuFromDQMuChange( dq_old, mu_old, dq_change, mu_change )
%UPDATEDQMUFROMDQMUCHANGE Summary of this function goes here
%   Detailed explanation goes here
    [R2_change, t2_change] = GenerateRtfromDQ(dq_change);
    dq_change_new = DQmult(dq_change, trans2dquat((1-mu_old)/mu_old*t2_change));
    dq_new = DQmult(dq_old,dq_change_new);
    mu_new = mu_old * mu_change;
end

