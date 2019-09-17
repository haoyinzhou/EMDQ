

function [ ind ] = LMR_func( X, Y, net )

    [N,~] = size(X);
    Xt = X';Yt = Y';
    
    Klist=[2,3,4, 5,6,7,...
           8,9,10,11,12,...
           13,15,18,20,25,...
           ];
    x=[];y=[];
    
    [feature]=MPC(Xt,Yt,Klist);
    x=[x;feature];

%     load('Net.mat')
    testVal0 = net(x.').';
    testVal = testVal0(:,1) >testVal0(:,2);
    ind = find( testVal==1 );

end

