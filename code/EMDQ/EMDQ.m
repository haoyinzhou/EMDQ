%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full code of EMDQ algorithm
% Zhou, Haoyin, and Jayender Jagadeesan. "Smooth Deformation Field-based Mismatch Removal in Real-time." 
% <CopyRight 2019> Haoyin Zhou, Jayender Jagadeesan
% Surgical Planning Laboratory
% Brigham and Women's Hospital
% Harvard Medical School
%
% This algorithm includes two steps: 
% (1) reweighting and 1-point RANSAC based outliers removal, which provide the initial value for the next step.
% (2) EM algorithm + dual quaternion-based deformation field generation, and feature matching outliers are also removed.
%
% This algorithm works well with large number of feature matches
%
% Input Pamameters:
% X1 and X2: the corresponding feature points. (D*N matrix, D = 2 or 3 is the dimension, N is the number of points)
% params: please find the file: EMDQ_Initialization_2D(or 3D).m
%
% Output Pamameters:
% inliersMask_RANSAC: outliers removal results of R1P-RNSC
% inliersMask_final: outliers removal results of EMDQ
% dq_points and mu_points: the dual-quertion and scale factors of the feature matches 
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inliersMask_R1P_RNSC, inliersMask_EMDQ, dq_points, mu_points] = EMDQ( X1, X2, params)
%% Iniitalize the outputs 
    N = size(X1,2); D = size(X1,1);
    if (D < 2 || D > 3) return; end
    if D == 2
        dqsize = 4;
    else
        dqsize = 8;
    end
    
    X1t = X1';    X2t = X2';   

    [inliersMask_R1P_RNSC,inliersMask_EMDQ] = deal(logical(zeros(1,N)));
    dq_points = [ones(N,1) zeros(N,dqsize-1)];
%     dq_points = [ones(N,1) zeros(N,1) (X1t-X2t)];
    mu_points = ones(N,1);
    
    %% R1P_RNSC   
    distance_all = 1e6 * ones(1,N);   
    SumRNSCWeight = ones(1,N);
    inliersRatio = 0.0;
    punish =  1e+6;
    MaxTry = min(100,N);
    ransactryid = randperm(N);
%     SuccessfulRANSACCount = 0;
    
    for ransactry = 1:MaxTry
        benchPid = ransactryid(ransactry);  % randomly select a point as the benchmark.
        if inliersMask_R1P_RNSC(benchPid) > 0
            continue;
        end
        if (punish < ransactry)
            break;
        end        
        X1_nobench = [X1(:,1:benchPid-1), X1(:,benchPid+1:end)];
        X2_nobench = [X2(:,1:benchPid-1), X2(:,benchPid+1:end)];
        S1 = X1_nobench - repmat(X1(:,benchPid),1,N-1);
        S2 = X2_nobench - repmat(X2(:,benchPid),1,N-1);

        R = eye(D); % initialiaztion
        mu = 1.0;   % initialiaztion      
        for iter = 1:3
            RS2 = mu * R * S2;
            distance_this = (sum((S1 - RS2).^2)).^0.5;
            p_this = params.inliersThreshold ./ distance_this; %exp(-(pixeldistance.^2) / (2 * sigma2));
            p_this(p_this > 1.0) = 1.0;
            p_this(p_this < 1e-8) = 1e-8;

            weightMatrix = repmat(p_this,D,1);
            wS1 = weightMatrix .* S1;
            wS2 = weightMatrix .* S2;
            [UR,SR,VR] = svd(wS1 * wS2');
            R = UR * VR';
%             mu = norm([wS1(1,:) wS1(2,:)]) / norm([wS2(1,:) wS2(2,:)]);
            mu = norm(reshape(wS1,[(N-1)*D,1])) / norm(reshape(wS2,[(N-1)*D,1]));
        end  
        
        distance_this = [distance_this(1:benchPid-1) 0.0 distance_this(benchPid:end)];
        inliersMask_this = distance_this < params.inliersThreshold;
        
        if sum(inliersMask_this) < params.LeastNumberOfInlierRANSACTrial 
            continue;
        end    
        if (D == 2)
            temp = unique(round(X2(:,inliersMask_this))','rows'); % this step aims to avoid the situation when multiple X1 match the same X2
            if size(temp,1) < params.LeastNumberOfInlierRANSACTrial
                continue;
            end
        end
%         SuccessfulRANSACCount = SuccessfulRANSACCount + 1;
        t = X1(:,benchPid) / mu - R * X2(:,benchPid);
%         fprintf('benchPid = %d, iter = %d, sum(p) = %.1f (%.1f)\n', benchPid, iter,sum(p_this), LeastSumP);
        dq_this = GenerateDQfromRt(R,t');
        idx = distance_this < distance_all & inliersMask_this;
        dq_points(idx,:) = repmat(dq_this, [sum(idx) 1]);
        mu_points(idx) = mu;
        distance_all(idx) = min(distance_this(idx), distance_all(idx));
        inliersMask_R1P_RNSC(idx) = 1;%

        p_this = exp(-distance_this.^2 / (2 * params.inliersThreshold^2)); % 
        SumRNSCWeight(idx) = sum(p_this);
        
        inliersRatio = sum(inliersMask_R1P_RNSC > 0) / N;        
        temp = 1-params.LeastNumberOfInlierRANSACTrial/((1-inliersRatio)*N+1e-8);
        temp = max(0.01,temp);
        punish = log(1-0.95) / log(temp); 
    end
 
    if (inliersRatio < 0.05)
         fprintf('inliersRatio = %.1f, too small number of inliers!\n', 100 * inliersRatio);          
         return;
    end    

    %% EMDQ algorithm, find the neighboring connections between points (KD tree needed)
    % generate Kraw (the weight matrix) based on the distances between points
    NeighborCount = min(params.NeighborCount, N-1);
    [Kraw1,neighborX1] = GenerateKDTreeWeightMatrix(X1, X1, NeighborCount, 1.0/(params.r2));
    [Kraw2,neighborX2] = GenerateKDTreeWeightMatrix(X2, X2, NeighborCount, 1.0/(params.r2));
%     Kraw1 = Kraw1(2:end,:); neighborX1 = neighborX1(2:end,:);
%     Kraw2 = Kraw2(2:end,:); neighborX2 = neighborX2(2:end,:);
    Kraw = [Kraw1; Kraw2];
    neighborX = [neighborX1;neighborX2];    
    Kraw = Kraw .* reshape(inliersMask_R1P_RNSC(neighborX)+0.1,size(neighborX)) + 1e-8;
    
    [Kraw,vid] = sort(Kraw,'descend');
    for i = 1:N
        neighborX(:,i) = neighborX(vid(:,i),i);
    end    
%     M = min(NeighborCount, median(sum(Kraw > 0.1 * median(Kraw(1,:)))) );
    Kraw = Kraw(2:NeighborCount+1,:);
    neighborX = neighborX(2:NeighborCount+1,:); 
    Kraw = Kraw';
    neighborX = neighborX';

    %% Initial EM from 1P-RANSAC
    sigma2 = params.inliersThreshold^2; % initialization
    gamma = inliersRatio; % initialization
    gamma = max(min(gamma,0.95),0.05);
    
    distance_all = distance_all';
    p = exp(-distance_all.^2 / (2 * sigma2)); % initialization
    p = p ./ (p + params.a * sigma2 *(1-gamma)/gamma );
    p = p .* SumRNSCWeight';    
    
    %% EM iterations
    p(p < 1e-5) = 1e-5;
    p_old = p;
    for iter = 1:params.EMMaxIterNumber
        % M-step   
        p_neighbor = p(neighborX);% .* p_RNSCConsensus(neighborX);    
        Kraw_p = p_neighbor .* Kraw;
        K = Kraw_p ./ (repmat(sum(Kraw_p,2), [1,NeighborCount]) + 1e-9);

        % E-step
        for i=1:dqsize
            temp = dq_points(:,i);
            dq_points(:,i) = sum(temp(neighborX) .* K, 2);
        end
        mu_points = sum(mu_points(neighborX) .* K, 2);

        X1_dq = WarpPosByDq(dq_points, X2t);
        X1_dq = X1_dq .* repmat(mu_points, [1,D]);
        Xdiff = X1t - X1_dq;
        distance2 = sum((Xdiff).^2,2);
        p = exp(-(distance2) / (2 * sigma2));
%         p = p ./ (p + params.a * sigma2 *(1-gamma)/gamma );
        p = p ./ (p + params.a * (2 * pi * sigma2)^(D/2) *(1-gamma)/gamma );
        p(p<1e-8) = 1e-8;

        sigma2 = sum(p .* distance2) / (sum(p) + 1e-8);        
        gamma = sum( (p>params.pThreshold) & (distance2 < params.inliersThreshold^2) ) / N;
        gamma = max(min(gamma,0.95),0.05);

        QTdiff = trans2dquat(Xdiff ./ repmat(mu_points,[1,D]) );
        dq_points = DQmult(dq_points,QTdiff);
        
        % termination
        pdiff = mean(abs(p - p_old));
        if iter > 1 && pdiff < params.pdiffthreshold
            break;
        end
        p_old = p;        
%         fprintf('%d: gamma = %f, sigma = %f, pdiff = %f\n', iter, gamma, sqrt(sigma2), pdiff);
    end
% toc
    inliersMask_EMDQ = (p > params.pThreshold)' & (distance2 < params.inliersThreshold^2)';
%     inliersMask_final = (p > params.pThreshold)';
end


function [K,neighborX] = GenerateKDTreeWeightMatrix(X, Y, NeighborCount, beta)
    [D,N] = size(Y);
    M = NeighborCount;
    kdtreeX = vl_kdtreebuild(X);
    [neighborX, ~] = vl_kdtreequery(kdtreeX, X, Y, 'NumNeighbors', NeighborCount);    
    X_neighbor = reshape(X(:,neighborX),[D*NeighborCount,N]);
   
    distance2 = repmat(Y,[M 1]) - X_neighbor;
    distance2 = distance2.^2;
    
    dis2 = zeros(M,N);
    for i = 1:D
        dis2 = dis2 + distance2(i:D:end,:);
    end
    K = exp(-beta * dis2) + 1e-8;
end




