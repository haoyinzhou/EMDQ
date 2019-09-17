function [ scale ] = GetDataScaleForEMDQ_3D( X1,X2 )

    [D,N] = size(X1);
    if D ~= 3
        fprintf('the data dimension should be 3.\n');
    end
    
    x1 = double(X1);
    x1_mean = mean(x1,2);    
    x1 = x1 - repmat(x1_mean,[1,N]);
    scale_x1 = sqrt(sum(sum(x1.^2))/N);

    x2 = double(X2);
    x2_mean = mean(x2,2);    
    x2 = x2 - repmat(x2_mean,[1,N]);
    scale_x2 = sqrt(sum(sum(x2.^2))/N);

    scale = (0.5 * (scale_x1 + scale_x2));        

end

