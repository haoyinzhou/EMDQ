% Visualization of the results.

function [ ] = ShowMatchOutliersRemovalResults( I, X1, X2, mask, Text )

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
    imshow(I);
    hold on
    for i = 1:size(X1_all,1)    
        temp = [X1_all(i,:);X2_all(i,:)];
        if (mask(i) > 0)
            plot(temp(:,1),temp(:,2),'y-','LineWidth', 1.5); 
        else
            plot(temp(:,1),temp(:,2),'k-','LineWidth', 1.5); 
        end
    end
    title(Text);
    
end

