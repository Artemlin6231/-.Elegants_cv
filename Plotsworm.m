location_x_read = 'vCoordv1_x.xls';
location_y_read = 'vCoordv1_y.xls';
X = readtable(location_x_read);
Y = readtable(location_y_read);

figure(4);
hold on;
for i=1:length(STATS)
    x = cell2mat(table2cell(X(:,i)));
    y = cell2mat(table2cell(Y(:,i)));
    x(x==0) = [];
    y(y==0) = [];
    plot(x,y);
    xlabel('X coordinates')
    ylabel('Y coordinates')
    
end

hold off;    
