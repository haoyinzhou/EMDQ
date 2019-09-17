% To generate the grid coord and the dual quaternion
%
% Output: 
% gridcoord and gridcoord_dq: the two coordinates of the grid displacments.
% dq_grid and mu_grid:        the dual quaternion and scale related to the grid points.
% gridmask:                   the mask.
%
% Input:
% coord_points:                coordinates of the feature matches
% dq_points and mu_points:     dual quaternion and scale of the feature matches
% gridmindis:                  minimum distance in pixel between the generated grid points
% params:                      the same params used in the main EMDQ algorithm, we only use r2
%                              to generate the weight matrix
% W and H:                     width and height of the image
% TH:                          a threshold


function [gridcoord, gridcoord_dq, dq_grid, mu_grid, gridmask] = GenerateGridDeformationField( coord_points, dq_points, mu_points, gridmindis, params, W, H, TH)
    [gridcoord] = GenerateGrid(W,H,gridmindis);

    N = size(coord_points,2);
    WeightMatrix_raw = GenerateDistanceWeightMatrix(gridcoord, coord_points, 1.0/(params.r2));
    gridmask = sum(WeightMatrix_raw,2) > TH;
%     gridmask = max(WeightMatrix_raw') > TH; %0.8; %1.5 * exp(-1.0/(params.r2) * 50^2);
    
    WeightMatrix_Point2Grid = WeightMatrix_raw ./ (repmat(sum(WeightMatrix_raw,2), [1,N]) + 1e-9);

    dq_grid = WeightMatrix_Point2Grid * dq_points;
    mu_grid = WeightMatrix_Point2Grid * mu_points;

    gridcoord = gridcoord';
    gridcoord_dq = WarpPosByDq(dq_grid, gridcoord);
    gridcoord_dq = gridcoord_dq .* repmat(mu_grid, [1,2]);

end

function [Y] = GenerateGrid(W,H,mindis)
    numberofsteps = (W - 1) / mindis;
    mindisnew = ceil(numberofsteps);
    stepnew = (W - 1) / mindisnew;
    c1 = round(1:stepnew:W);    
    numberofsteps = (H - 1) / mindis;
    mindisnew = ceil(numberofsteps);
    stepnew = (H - 1) / mindisnew;
    c2 = round(1:stepnew:H);    
    [x,y] = meshgrid(c1,c2);
    x = reshape(x, [1,size(x,1) * size(x,2)]);
    y = reshape(y, [1,size(y,1) * size(y,2)]);    
    Y = [x;y];
end


function [K] = GenerateDistanceWeightMatrix(x, y, beta)
    N = size(x,2); M = size(y,2);
    distance2 = repmat(x',[1 1 M]) - permute(repmat(y',[1 1 N]),[3 2 1]);
    distance2 = squeeze(sum(distance2.^2,2));
    K = exp(-beta * distance2) + 1e-8;
end
