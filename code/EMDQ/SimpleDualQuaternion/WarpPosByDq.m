function [ pos_out ] = WarpPosByDq( dq_input, pos_in)
    if (size(dq_input,2) == 8 && size(pos_in,2) == 3)
        dq_conj = DQconj(dq_input);
        Qpos = pos2dquat(pos_in);
        Q = DQmult(dq_conj,Qpos);
        Q = DQmult(Q,dq_input);   
        pos_out = Q(:,6:8) ./ repmat(sqrt(sum(Q(:,1:4).^2,2)),[1,3]);
    elseif (size(dq_input,2) == 4 && size(pos_in,2) == 2)
%         Qpos = pos2dquat(pos_in);
%         Q = DQmult_dqxdt(dq_input, Qpos);
%         pos_out = Q(:,3:4) ./ repmat(sqrt(sum(Q(:,1:2).^2,2)),[1,2]);
        dq_conj = DQconj(dq_input);
        Qpos = pos2dquat(pos_in);
        Q = DQmult(dq_conj,Qpos);
        Q = DQmult(Q,dq_input);    
        pos_out = Q(:,3:4) ./ repmat(sqrt(sum(Q(:,1:2).^2,2)),[1,2]);
    end
end

