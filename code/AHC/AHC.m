%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AHC - Augmented Homogeneous Coordinate Matrix based Mismatch Removal:
%   Detect the mismatches between two point sets.
%   Input:
%       X,Y: The coordinate of the two point sets.
%       et: The end_threshold for terminating the iterations.
%   Output:
%       TrueInd: The indices of the detected true matches.
%
%   Authors: 
%       Yan Zheng (yanzheng96117@163.com)
%   Date:
%       07/29/2017
%   Reference:
%       The Augmented Homogeneous Coordinates Matrix Based Projective Mismatch Removal for Partial-Duplicate Image Search
%       Yan Zheng and Zhouchen Lin



function  result=AHC(X,Y,et)  
     result=[];      
%% Remove the given matched pairs with repeated matched points in the target image and leave one pair randomly.
     [~,S_tmp,~]=unique(sum(Y));                                              
%% Iterative procedure for selecting anchor matches
     sigma=1.98; 
     MaxErr=100;                                                             
     while MaxErr>et                                                        
       [S_tmp,Max_Err]=AHCM(X,Y,S_tmp,sigma);
       if length(S_tmp)>4;Anchors=S_tmp;MaxErr=Max_Err;else;break;end
       sigma=sigma*0.98;                                                    
     end 
   %% get inliers
      if MaxErr<=et      
         result=AHCM(X,Y,Anchors,et);
      end
%       try
%         result = result;
%       catch
%         result = 0
%       end
end


%% The subfunction implements each step in the iteration. 
%   Input:
%       X,Y: The coordinate of the two point sets.
%       S: The index set of anchor matches for this iteration:
%       F: The threshold set to judge the reprojection error.
%   Output:
%       Anchors: The filtered anchor sets after each iteration.
%       MaxErr: The maximum reprojection error among the anchor matches.

function  [Anchors,MaxErr]=AHCM(X,Y,S,F)
    warning('off');
    
%% Estimating the Matched Points Based on Existing Anchor Matches
    X=double(X); Y=double(Y);
    MaxErr=100;
    AT=X(:,S); BT=Y(:,S); 
    
    u=[repmat(BT(1,:),3,1).*AT;AT];                                        
    L=inv(u*u');B=[X;X].*(L(1:3,:)'*X);
    a=-sum(B(4:6,:))./sum(B(1:3,:));

    v=[repmat(BT(2,:),3,1).*AT;AT]; 
    L=inv(v*v'); B=[X;X].*(L(1:3,:)'*X);
    b=-sum(B(4:6,:))./sum(B(1:3,:));
  
 %%  Filter the anchor matches by romoving pairs with reprojection errors below the iterated threshold
    r=Y(1:2,:)-[a;b];
    rr=sum(r.*r);
    if F>=2                                                       
        Anchors=find(rr<=F^2);
        return;
    else
      [~,mu,sigma]=zscore(r(:,S),[],2);                                   
       Z=[abs((r(1,:)-mu(1))/sigma(1));abs((r(2,:)-mu(2))/sigma(2))];                                       
       Anchors=find(Z(1,:)<=F & Z(2,:)<=F );
       MaxErr=sqrt(max(rr(Anchors)));if isempty( MaxErr);MaxErr=100;end    
     end
end