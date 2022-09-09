clear; clc;
vidObj = VideoReader('outputdata1.avi');
startFrame = 1;
endFrame = 248;
nFrames = 248;
CorrectFactor = 0.04;

vidObj.CurrentTime = 0;
erode_value = 0;
dilate_value = 5;
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);

v_new = VideoWriter('newfile1.avi','Motion JPEG AVI');
v_new.Quality = 100;
v_new.FrameRate = 10;
open(v_new);
for k = 1:nFrames
    Frame = k;
    result_numb = k;
    s(k).cdata = readFrame(vidObj);
    BW1 = s(k).cdata;
    BW2 = im2bw(BW1);
    se = strel('disk',5);
    L = imclose(BW2,se);
    BW3 = L;
    BW5 = bwareaopen(BW3,350);
    LS = bwlabeln(BW5);
    S = regionprops(LS,'Area');
    [L,NUM] = bwlabel(LS);
    STATS = regionprops(L, {'Area', 'Centroid', 'FilledArea', 'Eccentricity','Orientation','MinorAxisLength','MajorAxisLength','Image','BoundingBox'});
    FigureName = ['Tracking Results for Frame ', num2str(Frame)];
    RGB = label2rgb(L, @jet, 'k');
    figure(1)
    set(1, 'Name', FigureName);
    imshow(s(k).cdata);
    hold on
    for index = 1:length(STATS)
        xbar = STATS(index).Centroid(1);
        ybar = STATS(index).Centroid(2);
        a = STATS(index).MajorAxisLength/2;
        b = STATS(index).MinorAxisLength/2;
        theta = pi*STATS(index).Orientation/180;
        R = [cos(theta) sin(theta)
               -sin(theta) cos(theta)];
        xy = [a*cosphi; b*sinphi];
        xy = R*xy;
        x = xy(1,:)+xbar;
        y = xy(2,:)+ybar;
        wormImage = STATS(index).Image;
        BW8 = bwmorph(wormImage,'thin',Inf);

        v0 = find(bwperim(BW8),1);
        D0 = bwdistgeodesic(BW8,v0,'quasi-euclidean');
        [v1y,v1x]=find(D0 == max(max(D0)), 1);
        D1 = bwdistgeodesic(BW8,v1x,v1y,'quasi-euclidean');
        gdPX = max(max(D1));
        [v2y,v2x] = find(D1 == gdPX,1);
        vCoordv1x = v1x + STATS(index).BoundingBox(1);
        vCoordv1y = v1y + STATS(index).BoundingBox(2);
        vCoordv2x = v2x + STATS(index).BoundingBox(1);
        vCoordv2y = v2y + STATS(index).BoundingBox(2);
           
           % Record the path (between v1 and v2) of the graph diameter
        D2 = bwdistgeodesic(BW8, v2x, v2y, 'quasi-euclidean');
        D = round((D1 + D2)) / 8;
        D(isnan(D)) = Inf;
        path = bwmorph(imregionalmin(D),'thin',Inf);
        [Py,Px] = find(path);

        [B] = bwtraceboundary(imregionalmin(D),[v1y v1x],'E');
        cent_num = round(numel(B(:,1))/4);
        wormcenter_x = B(cent_num,2)+STATS(index).BoundingBox(1);
        wormcenter_y = B(cent_num,1)+STATS(index).BoundingBox(2); 
        location_x(result_numb,index) = wormcenter_x;
        location_y(result_numb,index) = wormcenter_y;
        vCoordv1_x(result_numb,index) = vCoordv1x;
        vCoordv1_y(result_numb,index) = vCoordv1y;
        vCoordv2_x(result_numb,index) = vCoordv2x;
        vCoordv2_y(result_numb,index) = vCoordv2y;
        wormskel_x = B(:,2)+STATS(index).BoundingBox(1);
        wormskel_y = B(:,1)+STATS(index).BoundingBox(2);
        
        plot(wormskel_x,wormskel_y,'r',wormcenter_x,wormcenter_y,'ow',vCoordv1x,vCoordv1y,'ow',vCoordv2x,vCoordv2y,'ow');
        hold on
     end

     data_matrix = zeros(size(location_x(result_numb,:)));
     temp_data_x = zeros(size(location_x(result_numb,:)));
     temp_data_y = temp_data_x;
     temp_data_v1x = temp_data_x;
     temp_data_v1y = temp_data_x;
     temp_data_v2x = temp_data_x;
     temp_data_v2y = temp_data_x;
     zero_i=[]; 
   
        if k > startFrame
                lonely_point_count = 1;
                for i = 1:numel(location_x(result_numb,:))
                    count_zeros = 0;
                    for j = 1:numel(location_x(result_numb,:))
                        if abs(location_x(result_numb,i) - location_x(result_numb-1,j)) < 50 && abs(location_y(result_numb,i)-location_y(result_numb-1,j)) < 50
                            temp_data_x(j) = location_x(result_numb,i);
                            temp_data_y(j) = location_y(result_numb,i);
                            temp_data_v1x(j) = vCoordv1_x(result_numb,i);
                            temp_data_v1y(j) = vCoordv1_y(result_numb,i);
                            temp_data_v2x(j) = vCoordv2_x(result_numb,i);
                            temp_data_v2y(j) = vCoordv2_y(result_numb,i);
                        else
                            count_zeros=count_zeros+1;
                        end
                    end
                    if count_zeros==numel(location_x(result_numb,:))
                        zero_i(lonely_point_count) = i;
                        lonely_point_count = lonely_point_count + 1;
                    end
                end
                clear zero_i
            location_x(result_numb,:) = temp_data_x;
            location_y(result_numb,:) = temp_data_y;
            vCoordv1_x(result_numb,:) = temp_data_v1x;
            vCoordv1_y(result_numb,:) = temp_data_v1y;
            vCoordv2_x(result_numb,:) = temp_data_v2x;
            vCoordv2_y(result_numb,:) = temp_data_v2y;
        end
        plot(location_x(result_numb,:),location_y(result_numb,:),'+r')

        writeVideo(v_new,getframe(gcf));
        hold off
end    % END for Frame
close(gcf);
close(v_new);
xlswrite('location_x',location_x(startFrame:endFrame,:));
xlswrite('location_y',location_y(startFrame:endFrame,:));
xlswrite('vCoordv1_x',vCoordv1_x(startFrame:endFrame,:));
xlswrite('vCoordv1_y',vCoordv1_y(startFrame:endFrame,:));
xlswrite('vCoordv2_x',vCoordv2_x(startFrame:endFrame,:));
xlswrite('vCoordv2_y',vCoordv2_y(startFrame:endFrame,:));