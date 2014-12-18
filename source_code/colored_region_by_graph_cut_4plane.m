function [im, color_mask] = colored_region_by_graph_cut_4plane(InFileName,start_frame,blksize,problist)

bi_lat_filtering = 0;

end_frame = start_frame; 
frame_index = start_frame:end_frame; %1:50; % number of frames % 300
InSequence=aviread(InFileName, frame_index); %read in
FileInfo=aviinfo(InFileName);

width = FileInfo.Width;
height = FileInfo.Height;

im = frame2im(InSequence(1));
 
%resize image
im1 = im;
%pause


%im1 = imread('under_seg.ppm');  %test
%im1 = imresize(im, [height/blksize width/blksize]);
%problist = imresize(problist, blksize); 

if (bi_lat_filtering == 0)
    %%%%% convert to LUV color space
    im2 = conv_LUV(im1);
else
    %%%% bilateral filtering does 
    im1 = double(im1/255.0);
    % Set bilateral filter parameters.
    w     = 5;       % bilateral filter half-width
    sigma = [3 0.1]; % bilateral filter standard deviations
    im2 = bfilter2(im1);
    im2 = double(im2*255.0);
    im2 = uint8(im2);
end

% im3 = zeros(size(problist,1), size(problist,2), 4);
% im3 = uint8(im3);
% im3(:,:,1:3) = im2;
% scaled_problist = double(problist)*128; %256
% im3(:,:,4) = uint8(scaled_problist);  %problist

% sigma: to smooth the image.
% c: constant for threshold function.
% min_size: minimum component size (enforced by post-processing stage).
% im: image to segment.
if ( bi_lat_filtering == 0)
    %color_mask = segmentImgOpt_4plane( 0.4, 80, 100, im2, 'ouput.ppm', 1);  %(0.5, 100, 100)  %im2
    %color_mask = segmentImgOpt( 0.5, 100, 100, im2, 'ouput.ppm', 1);  %(0.5, 200, 100)
    %color_mask = segmentImgOpt( 0.4, 3, 10, im2, 'ouput.ppm', 1);
    color_mask = segmentImgOpt( 0.4, 100, 100, im2, 'ouput.ppm', 1); %%(0.5, 200, 100)
else
    color_mask = segmentImgOpt_4plane( 0, 100, 50, im2, 'ouput.ppm', 1);  %(0.5, 200, 100)  %im2
    %color_mask = segmentImgOpt( 0.5, 100, 100, im2, 'ouput.ppm', 1);  %(0.5, 200, 100)
    %color_mask = segmentImgOpt( 0.4, 3, 10, im2, 'ouput.ppm', 1);
end


%color_mask = imresize(color_mask, blksize); 

end





function im2 = conv_LUV(im)
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
end
