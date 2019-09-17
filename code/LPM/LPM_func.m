%% This is the function version of Ma's original code

function [ CorrectIndex ] = LPM_func( X, Y )

    lambda1   = 0.8;   lambda2   = 0.5;   
    numNeigh1 = 6;     numNeigh2 = 6; 
    tau1      = 0.2;   tau2      = 0.2;   

    x1 = X; y1 = Y;
    [numx1,~] = size(x1);
    p1 = ones(1,numx1);
    Xt = X';Yt = Y';
    vec=Yt-Xt;
    d2=vec(1,:).^2+vec(2,:).^2;


    %%  iteration 1
    % % % constructe K-NN by kdtree
    kdtreeX = vl_kdtreebuild(Xt);
    kdtreeY = vl_kdtreebuild(Yt);  
    [neighborX, ~] = vl_kdtreequery(kdtreeX, Xt, Xt, 'NumNeighbors', numNeigh1+3) ;
    [neighborY, ~] = vl_kdtreequery(kdtreeY, Yt, Yt, 'NumNeighbors', numNeigh1+3) ;

    % % % calculate the locality costs C and return binary vector p
    [p2, C] = LPM_cosF(neighborX, neighborY, lambda1, vec, d2, tau1, numNeigh1);

    %%  iteration 2
    idx = find( p2 == 1 );
    if length(idx)>= numNeigh2+4
        kdtreeX = vl_kdtreebuild(Xt( :, idx ));
        kdtreeY = vl_kdtreebuild(Yt( :, idx ));
        [neighborX, ~] = vl_kdtreequery(kdtreeX, Xt(:,idx), Xt, 'NumNeighbors', numNeigh2+3) ;
        [neighborY, ~] = vl_kdtreequery(kdtreeY, Yt(:,idx), Yt, 'NumNeighbors', numNeigh2+3) ;
        neighborX = idx(neighborX);
        neighborY = idx(neighborY);
        [p2, C] = LPM_cosF(neighborX, neighborY, lambda2, vec, d2, tau2, numNeigh2);
    end

    CorrectIndex = find(p2 == 1);

end

