%temp_test
im1 = double(im)/255;
im1 = COLORSPACE('Luv<-RGB',im1);

[x y c] = size(im);
im2 = zeros(x,y,c);

range = max(max(im1(:,:,1))) - min(min(im1(:,:,1)));
range = double(range);
im2(:,:,1) = (im1(:,:,1)- min(min(im1(:,:,1))));
im2(:,:,1) = (im2(:,:,1)/range);


range = max(max(im1(:,:,2))) - min(min(im1(:,:,2)));
range = double(range);
im2(:,:,2) = (im1(:,:,2)- min(min(im1(:,:,2))));
im2 (:,:,2) = (im2(:,:,2)/range); 

range = max(max(im1(:,:,3))) - min(min(im1(:,:,3)));
range = double(range);
im2(:,:,3) = (im1(:,:,3)- min(min(im1(:,:,3))));
im2(:,:,3) = (im2(:,:,3)/range); 

im2 = double(im2*255);
%im2 = im2;

im2 = uint8(im2);

%figure(1), imshow(im2(:,:,1), []), title('LUV space : L');

color_mask = segmentImgOpt( 0.5, 100, 50, im, 'output_rgb.ppm', 1);

color_mask = segmentImgOpt( 0.3, 100, 20, im2, 'output_luv.ppm', 1);


