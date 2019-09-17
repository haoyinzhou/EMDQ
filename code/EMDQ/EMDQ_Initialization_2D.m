function [ params ] = EMDQ_Initialization_2D(scale)

    %% RANSAC-related parameters
    % default value: inliersThreshold = 20 pixels for 800x600 reslution,
    % "scale" is used to adapt to differet image resolutions.  
    % inliersThreshold is used in both RANSAC and EM process to determine inliers
    params.inliersThreshold = 20.0 * scale;
    
    % default value: LeastNumberOfInlierRANSACTrial = 5
    % In a RANSAC trial, only when the detected number of inliers is larger
    % than this threshold, then the results of this trial will be further considered.
    params.LeastNumberOfInlierRANSACTrial = 5; 
    
    %% The number of neighboring points found by KD-tree 
    params.NeighborCount = 16;
    
    %% To build connections between feature matches
    % default value: r = 80 pixels for 800x600 reslution.
    % The pixel distance between points for computing weights =
    % 1.0/(params.r2) between points.
    params.r2 = (80 * scale)^2; 
    
    %% The parameter for computing p of points in the EM algorithm
    % default value: a = 1e+5 
    params.a =  1e-5 / scale^2; %2*pi/(1e+5) / scale;
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
    
    %% Number of points used in sparse RANSAC
    params.N_sparse = 100;
    
end


