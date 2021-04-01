% this is the test code of this simplified dual quaternion box.
% maybe it is helpful for beginners to better understand dual quaternion
% I give some simple examples, please uncomment them and run the code.
% Haoyin Zhou, Surgical Planning Lab, Brigham and Women's Hospital, Harvard
% Medical School. 
% email: zhouhaoyin@bwh.harvard.edu or zhouhaoyin@gmail.com

clear
clc
close all

%% 2D test 1: dq-based transform is equal to the [R t] transform
% angle = randn(1,1);
% R = angle2rotMatrix(angle);
% t = 100 * randn(2,1);
% 
% X2 = [5 10];
% X1_Rt = R * X2' + t;
% 
% dq = GenerateDQfromRt( R, t');
% X1_dq = WarpPosByDq(dq, X2);
% X1_dq = X1_dq';
% 
% [X1_Rt X1_dq] % we have X1_Rt = X1_dq
% X3 = 100 * randn(2,1);
% Xdiff = X3 - X1_dq;
% QTdiff = trans2dquat(Xdiff');
% dqnew = DQmult(dq,QTdiff);
% X2_update = WarpPosByDq(dqnew, X2);
% [X3, X2_update']
% 
%% 2D test 2: dq-based interplation is smooth 
% N = 100;
% x = 1:N; x = x';
% y = zeros(length(x),1);
% C = [x y];
% 
% t = zeros(N,2);
% angle = zeros(N,1);
% R = zeros(2,2,N);
% warp = zeros(N, 4);
% 
% t(1,:) = 100 * randn(1,2)+ 2;%[0 10];
% t(N,:) = 100 * randn(1,2);%[0 5];
% angle(1) = randn(1,1);
% angle(N) = randn(1,1);
% 
% R(:,:,1) = [cos(angle(1)) -sin(angle(1)); sin(angle(1)) cos(angle(1))];
% R(:,:,N) = [cos(angle(N)) -sin(angle(N)); sin(angle(N)) cos(angle(N))];
% 
% warp(1,:) = GenerateDQfromRt(R(:,:,1),t(1,:));
% warp(N,:) = GenerateDQfromRt(R(:,:,N),t(N,:));
% for i = 1:N
%     weight1 = N - i;
%     weightN = i - 1;    
%     sumweight = weight1 + weightN;
%     weight1 = weight1 / sumweight;
%     weightN = weightN / sumweight;    
%     warp(i,:) = weight1 * warp(1,:) + weightN * warp(N,:);
% end
% P1 = WarpPosByDq(warp, C);
% 
% figure
% set(gcf, 'color', 'w');
% hold on
% plot(C(:,1),C(:,2),'ro');
% plot(P1(:,1),P1(:,2),'b*');
% axis equal
% grid on

%% 3D test 1: dq-based transform is equal to the [R t] transform
% angle = randn(1,3);
% 
% R = angle2rotMatrix(angle); %angle2dcm(angle(1),angle(2),angle(3));
% t = 100 * randn(3,1);
% 
% X2 = [5 10 3];
% X1_Rt = R * X2' + t;
% 
% dq = GenerateDQfromRt( R, t');
% X1_dq = WarpPosByDq(dq, X2);
% X1_dq = X1_dq';
% 
% [X1_Rt X1_dq] % we have X1_Rt = X1_dq
% 
% % we test how to update a dq from a new translation
% X3 = 100 * randn(3,1);
% Xdiff = X3 - X1_dq;
% QTdiff = trans2dquat(Xdiff');
% dqnew = DQmult(dq,QTdiff);
% X2_update = WarpPosByDq(dqnew, X2);
% 
% [X3, X2_update']

%% 3D test 2: dq-based interplation is smooth 
% N = 100;
% x = 1:N; x = x';
% y = 1:N; y = y';
% z = zeros(N,1);
% C = [x y z];
% 
% mu = 0.5;
% t = zeros(N,3);
% angle = zeros(N,3);
% R = zeros(3,3,N);
% warp = zeros(N, 8);
% 
% t(1,:) = 100 * randn(1,3);%[0 10];
% t(N,:) = 100 * randn(1,3);%[0 5];
% angle(1,:) = randn(1,3);
% angle(N,:) = randn(1,3);
% 
% R(:,:,1) = angle2rotMatrix(angle(1,:));%angle2dcm(angle(1,1),angle(1,2),angle(1,3));
% R(:,:,N) = angle2rotMatrix(angle(N,:));%angle2dcm(angle(N,1),angle(N,2),angle(N,3));
% 
% warp(1,:) = GenerateDQfromRt(R(:,:,1),t(1,:));
% warp(N,:) = GenerateDQfromRt(R(:,:,N),t(N,:));
% for i = 1:N
%     weight1 = N - i;
%     weightN = i - 1;    
%     sumweight = weight1 + weightN;
%     weight1 = weight1 / sumweight;
%     weightN = weightN / sumweight;    
%     warp(i,:) = weight1 * warp(1,:) + weightN * warp(N,:);
% end
% P1 = WarpPosByDq(warp, C);
% 
% figure
% hold on
% plot3(C(:,1),C(:,2),C(:,3),'ro');
% plot3(P1(:,1),P1(:,2),P1(:,3),'b*');
% axis equal
% grid on

%% test conversion among dual quaternion, rotMatrix, angles and translations.
% angle = [-0.919695072404640,-0.706421954056373,-0.642052446655979];%randn(1,3);
% t = [-151.973772137986,-174.539386354885,-312.890661533416];%100 * randn(1,3);
% % 
% % angle = -0.779440665402206;%randn(1,1);
% % t = [52.0463114879734,-104.080086146455];%100 * randn(1,2);
% % 
% R = angle2rotMatrix(angle);
% 
% dq = rotMatrix2dquat(R);
% R2 = dquat2rotMatrix(dq);
% [R R2]
% 
% dq3 = angle2dquat(angle); 
% [dq;dq3]
% 
% angle2 = rotMatrix2angle(R);
% [angle; angle2]
% % 
% % R3 = angle2rotMatrix(angle2);
% % [R R2 R3]
% 
% % dq = trans2dquat(t);
% % t2 = dquat2trans(dq);
% % [t; t2]
% 
% % dq = GenerateDQfromRt(R, t);
% % [R2, t2] = GenerateRtfromDQ(dq);
% % [R R2]
% % [t; t2]


%% DQ inv test
% % angle = randn(1,3);
% % t = 100 * randn(3,1);
% % X2 = 1000*randn(1,3);
% angle = randn(1,1);
% t = 100 * randn(2,1);
% X2 = 100*randn(1,2);
% 
% R = angle2rotMatrix(angle);
% 
% X1_Rt = R * X2' + t;
% 
% dq = GenerateDQfromRt( R, t');
% X1_dq = WarpPosByDq(dq, X2);
% X1_dq = X1_dq';
% [X1_Rt X1_dq]
% 
% dq_inv = DQinv(dq);
% X1 = X1_Rt';
% 
% X2_dqinv = WarpPosByDq(dq_inv, X1);
% [X2; X2_dqinv]


%% dq propogation test

% X0 = [1;2;3];
% angle1 = randn(1,3);
% t1 = 100 * randn(3,1);
% R1 = angle2rotMatrix(angle1);
% dq1 = GenerateDQfromRt( R1, t1');
% mu1 = 0.5;
% 
% angle2 = randn(1,3);
% t2 = 100 * randn(3,1);
% R2 = angle2rotMatrix(angle2);
% dq2 = GenerateDQfromRt( R2, t2');
% mu2 = 1.2;
% 
% X1_rt = mu1 * (R1 * X0 + t1);
% X1_dq = mu1 * WarpPosByDq(dq1, X0');
% 
% X2_rt = mu2 * (R2 * X1_rt + t2);
% X2_dq = mu2 * WarpPosByDq(dq2, X1_dq);
% 
% % method 2
% % R_all = R2 * R1;
% % t_all = R2 * t1 + t2 / mu1;
% % dq_all = GenerateDQfromRt( R_all, t_all');
% 
% %method 3
% % [R2_dq, t2_dq] = GenerateRtfromDQ(dq2);
% % dq2_new = DQmult(dq2, trans2dquat((1-mu1)/mu1*t2_dq));
% % dq_all = DQmult(dq1,dq2_new);
% % mu_all = mu1 * mu2;
% 
% [ dq_all, mu_all ] = UpdateDQMuFromDQMuChange( dq1, mu1, dq2, mu2 );
% 
% X_dq_all = mu_all * WarpPosByDq(dq_all, X0');
% 
% [X2_rt'; X2_dq; X_dq_all]

%% dq inv propogation test
% X0 = [1;2];
% angle1 = randn(1,1);
% t1 = 100 * randn(2,1);
% R1 = angle2rotMatrix(angle1);
% dq1 = GenerateDQfromRt( R1, t1');
% mu1 = 0.5;
% 
% angle2 = randn(1,1);
% t2 = 100 * randn(2,1);
% R2 = angle2rotMatrix(angle2);
% dq2 = GenerateDQfromRt( R2, t2');
% mu2 = 1.2;
% 
% delta_mu = mu2 / mu1;
% delta_R = R2 * R1';
% delta_t = mu1 * (t2 - delta_R * t1);
% delta_dq = GenerateDQfromRt(delta_R, delta_t');
% 
% X2_rt = mu2 * (R2 * X0 + t2);
% X2_dq = mu2 * WarpPosByDq(dq2, X0');
% 
% [ dq_all, mu_all ] = UpdateDQMuFromDQMuChange( dq1, mu1, delta_dq, delta_mu );
% X2_all = mu_all * WarpPosByDq(dq_all, X0');
% 
% [X2_rt'; X2_dq; X2_all]


%% 2D and 3D relationship

angle2d = 0.1;%randn(1,1);
R2d = angle2rotMatrix(angle2d);
t2d = [1.0 2.0];%]100 * randn(1,2);

X2d = [5.0 8.0];
X2d_Rt = R2d * X2d' + t2d';

dq2d = GenerateDQfromRt( R2d, t2d);
X2d_dq = WarpPosByDq(dq2d, X2d);

[X2d_Rt X2d_dq'] 


angle3d = [0 0 angle2d];
R3d = angle2rotMatrix(angle3d); % = [R2d 0; 0 1]
t3d = [t2d, 0];
X3d = [X2d 0];
dq3d = GenerateDQfromRt( R3d, t3d);

X3d_Rt = R3d * X3d' + t3d';
X3d_dq = WarpPosByDq(dq3d, X3d);
[X3d_Rt X3d_dq'] 



%% 2D
% N = 101;
% x2d = (1:N); x2d = x2d';
% % y = 1:N; y = y';
% y2d = zeros(N,1);
% 
% C2d = [x2d y2d];
% 
% figure
% set(gcf, 'color', 'w');
% hold on
% plot(C2d(:,1),C2d(:,2),'ro');
% plot(C2d(1:5:end,1),C2d(1:5:end,2),'k*');
% 
% % a = 80;
% for k = 1:5
%     a = 12 * (k+1); 
%     t2d = zeros(N,2);
%     angle2d = zeros(N,1);
%     R2d = zeros(2,2,N);
%     warp2d = zeros(N, 4);
% 
%     angle2d(1) =  (a) * pi / 180;
%     angle2d(N) =  (-a) * pi / 180;
% 
%     R2d(:,:,1) = [cos(angle2d(1)) -sin(angle2d(1)); sin(angle2d(1)) cos(angle2d(1))];
%     R2d(:,:,N) = [cos(angle2d(N)) -sin(angle2d(N)); sin(angle2d(N)) cos(angle2d(N))];
%     t2d(1,:) = C2d(1,:)' -R2d(:,:,1) * C2d(1,:)' + [-10 10]';
%     t2d(N,:) = C2d(N,:)' -R2d(:,:,N) * C2d(N,:)' + [ 10 10]';
% 
%     warp2d(1,:) = GenerateDQfromRt(R2d(:,:,1),t2d(1,:));
%     warp2d(N,:) = GenerateDQfromRt(R2d(:,:,N),t2d(N,:));
%     for i = 1:N
%         weight1 = N - i;
%         weightN = i - 1;    
%         sumweight = weight1 + weightN;
%         weight1 = weight1 / sumweight;
%         weightN = weightN / sumweight;    
%         warp2d(i,:) = weight1 * warp2d(1,:) + weightN * warp2d(N,:);
%     end
%     P12d = WarpPosByDq(warp2d, C2d);
%     
%     plot(P12d(:,1),P12d(:,2),'b*-');
%     
%     if (k == 1)
%         for j = 1:5:N
%             diff = P12d(j,:) - C2d(j,:);
%             diff = 1.1 * diff;
%             h1 = quiver(C2d(j,1), C2d(j,2), diff(1), diff(2),'k-','LineWidth', 2.0,'MaxHeadSize',1.5);
%             temp = [P12d(j,:); C2d(j,:)];
%             plot(temp(:,1),temp(:,2),'k-','LineWidth', 1.5); 
%         end
%     end
% end
% 
% axis equal
% grid on
% 
%% 3D
% N = 101;
% x3d = (1:N); x3d = x3d';
% y3d = zeros(N,1);%0.1*(1:N); y3d = y3d';
% z3d = zeros(N,1);
% 
% C3d = [x3d y3d z3d];
% 
% figure
% set(gcf, 'color', 'w');
% hold on
% plot3(C3d(:,1),C3d(:,2),C3d(:,3),'ro');
% plot3(C3d(1:5:end,1),C3d(1:5:end,2),C3d(1:5:end,3),'k*');
% 
% for k = 1:7
%     a = 10 * (k-4); 
%     t3d = zeros(N,3);
%     angle3d = zeros(N,3);
%     R3d = zeros(3,3,N);
%     warp3d = zeros(N, 8);
% 
%     angle3d(1,:) =  [ 50  0   a] * pi / 180;
%     angle3d(N,:) =  [-50 -0  -a] * pi / 180;
% 
%     R3d(:,:,1) = angle2rotMatrix(angle3d(1,:));% angle2dcm(angle3d(1,1),angle3d(1,2),angle3d(1,3));
%     R3d(:,:,N) = angle2rotMatrix(angle3d(N,:));% angle2dcm(angle3d(N,1),angle3d(N,2),angle3d(N,3));
% 
%     % t(1,:) = randn(1,3);
%     % t(N,:) = randn(1,3);
%     t3d(1,:) = C3d(1,:)' -R3d(:,:,1) * C3d(1,:)' + [-0 00 10]';
%     t3d(N,:) = C3d(N,:)' -R3d(:,:,N) * C3d(N,:)' + [ 0 00 10]';
% %     t(1,:) = [0 0]';
% %     t(N,:) = [0 0]';
% 
%     warp3d(1,:) = GenerateDQfromRt(R3d(:,:,1),t3d(1,:));
%     warp3d(N,:) = GenerateDQfromRt(R3d(:,:,N),t3d(N,:));
%     for i = 1:N
%         weight1 = N - i;
%         weightN = i - 1;    
%         sumweight = weight1 + weightN;
%         weight1 = weight1 / sumweight;
%         weightN = weightN / sumweight;    
%         warp3d(i,:) = weight1 * warp3d(1,:) + weightN * warp3d(N,:);
%     end
%     P13d = WarpPosByDq(warp3d, C3d);
% 
%     plot3(P13d(:,1),P13d(:,2),P13d(:,3),'b*-');
%     
%     if (k == 7)
%     for j = 1:5:N
%         diff = P13d(j,:) - C3d(j,:);
%         diff = 1.1 * diff;
%         h1 = quiver3(C3d(j,1), C3d(j,2), C3d(j,3), diff(1), diff(2), diff(3),'k-','LineWidth', 2.0,'MaxHeadSize',0.5);
%         temp = [P13d(j,:); C3d(j,:)];
%         plot3(temp(:,1),temp(:,2),temp(:,3),'k-','LineWidth', 1.5); 
%     end
%     end
% 
% end
% 
% axis equal
% grid on



