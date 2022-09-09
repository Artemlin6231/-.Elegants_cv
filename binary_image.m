I = imread('outputdata\binary_images_manual070011.tiff');
imshow(I);
I1 = max(I(:))-I;
g1=double(I1)./(double(I1)+10);

%imshow(g1);

se = strel('line',3,5)
background = imopen(I1,se);
%imshow(background);
I2 = I - background;
%imshow(I2)
I3 = imadjust(I2);
J = imnoise(I1,'gaussian',0,0.0001);
imshow(J);
K = wiener2(J,[5 5]);
figure(1);
%imshow(K);
BW=im2bw(I1, 0.1);
%imshow(BW)
imshowpair(I1,BW, 'montage');
Ifilled=imfill(BW,'holes');
