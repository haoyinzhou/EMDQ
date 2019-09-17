function [ ] = ShowMatchOutliersRemovalResults_3D(Pts3D1, Pts3D2, X1, X2, mask, Text )

    [d1,d2] = size(X1);
    if (d1 == 2 || d1 == 3)
        X1_all = X1';
        X2_all = X2';
    else
        X1_all = X1;
        X2_all = X2;        
    end
    
    figure
    set(gcf, 'color', 'w');
    hold on
    plot3(X1_all(:,1),X1_all(:,2),X1_all(:,3),'.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.5);
    plot3(X2_all(:,1),X2_all(:,2),X2_all(:,3),'.', 'Color', [0.2 0.5 0.8], 'MarkerSize', 0.5);
    plot3(Pts3D1(1:5:end,1),Pts3D1(1:5:end,2),Pts3D1(1:5:end,3),'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
    plot3(Pts3D2(1:5:end,1),Pts3D2(1:5:end,2),Pts3D2(1:5:end,3),'o', 'Color', [0.2 0.5 0.8], 'MarkerSize', 0.1);
    for i = 1:size(X1_all,1)    
        temp = [X1_all(i,:);X2_all(i,:)];
        if (mask(i) > 0)
            plot3(temp(:,1),temp(:,2),temp(:,3),'b-','LineWidth', 2); 
        else
            plot3(temp(:,1),temp(:,2),temp(:,3),'k-','LineWidth', 0.01); 
        end
    end
%     grid on
    axis equal
    title(Text);
end

