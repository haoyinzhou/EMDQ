function [] = ShowGridDeformationField(I, gridcoord, gridcoord_dq, gridmask)
    diff = gridcoord_dq - gridcoord;
    diff = 1.1 * diff; % only for matlab display, because the arrow is shorter
                       
    figure;
    set(gcf, 'color', 'w');
    imshow(I);
    hold on
    for i = 1:size(gridcoord_dq,1)
%          temp = [gridcoord(i,:);gridcoord_dq(i,:)];
         if (gridmask(i) > 0)
%              plot(temp(:,1),temp(:,2),'y-','LineWidth', 0.5); 
            h1 = quiver(gridcoord(i,1), gridcoord(i,2), diff(i,1), diff(i,2),'y-','LineWidth', 2.0,'MaxHeadSize',200);
%             quiver(x, y, dxda, dyda)
%             set(h1,'AutoScale','on', 'AutoScaleFactor', 5);
         else
%             quiver(gridcoord(i,1), gridcoord(i,2), diff(i,1), diff(i,2),'r-', 'MarkerSize', 10);
         end
    end
%     plot(gridcoord(:,1),gridcoord(:,2),'ro','MarkerSize', 3,'LineWidth', 2); 
    plot(gridcoord(gridmask,1),gridcoord(gridmask,2),'bo','MarkerSize', 3,'LineWidth', 1.5); 
    
end


