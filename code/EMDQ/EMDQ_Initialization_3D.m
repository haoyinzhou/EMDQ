function [ params ] = EMDQ_Initialization_3D(scale)

    %% RANSAC-related parameters
    % default value: inliersThreshold = 0.1 of the scale.  
    % used in both RANSAC and EM process to determine inliers
    params.inliersThreshold = 0.1  * scale;
    
    % default value: LeastNumberOfInlierRANSACTrial = 5
    % In a RANSAC trial, only when the detected number of inliers is larger
    % than this threshold, then the results of this trial will be further considered.
    params.LeastNumberOfInlierRANSACTrial = 5; 
    
    %% The number of neighboring points found by KD-tree 
    % Please try to use a larger number if you think the algorithm is not stable.
    params.NeighborCount = 50; % this is because in 3D space
                                % more neighboring points may be outliers
                                %, so we use a larger number
                                
    %% To build connections between feature matches
    % default value: r = 0.3 of the scale.
    params.r2 = (0.3 * scale)^2; 
    
    %% The parameter for computing p of points in the EM algorithm
    % default value: a = 20.0 / scale^3 of the scale
    % p = exp() / (exp() + a * (2 * pi * sigma2)^(D/2) * (1-gamma)/gamma)
    params.a = 20.0 / scale^3;%
    %% The termination conditions of the EM algorithm
    % default value: EMMaxIterNumber = 15, which is the maximum number of
    % EM iterations
    params.EMMaxIterNumber = 15; 

    % default value: pdiffthreshold = 0.005,
    % when mean(abs(p - p_old)) < pdiffthreshold, the EM iteration stops.     
    params.pdiffthreshold = 0.005;  
    
    %% The parameter to determine whethe a point is an inlier or not
    % default value: pThreshold = 0.5, which means this feature match
    % is more likely to be an inlier
    params.pThreshold = 0.5; 
    
    %%
    params.N_sparse = 100;
    
end


