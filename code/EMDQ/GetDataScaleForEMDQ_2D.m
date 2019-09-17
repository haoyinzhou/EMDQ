function [ scale ] = GetDataScaleForEMDQ_2D( imagesize_W,imagesize_H )
    % for 2D image feature matching, we use 800x600 reslution images as the reference resolution.
    scalex = imagesize_W / 800;
    scaley = imagesize_H / 600;
    scale = 0.5 *(scalex + scaley);
end

