% PROJ : depth from motion
% simple motion estimation/segmentation by phase correlation method
% Author : Subarna Tripathi
% Date : May 2010
% example use : 
% [blklist, problist] = motion_segmentation_new ('test.avi', 4, 3, 4, edge_mask, color_mask);

function [blklist, problist, edge_mask] = motion_segmentation_new (InFileName, blksize, start_frame, end_frame, edge_mask)

%InFileName = 'akio_cif.avi'; %'films_CIF_part2_fr16.avi'; %'akio_cif.avi'; %'foreman_cif.avi'; %'renata_qcif.avi'; %'racingcar.avi'; 
% -------- initialize variables --------
flag = 0;
if( start_frame > end_frame)
    flag = 1;
    temp = start_frame;
    start_frame = end_frame;
    end_frame = temp;
end

frame_index = start_frame:end_frame; %1:50; % number of frames % 300
last_index = size(frame_index);
last_index = last_index(2);

InSequence=aviread(InFileName, frame_index); %read in
%InSequence=aviread(InFileName); %read in

FileInfo=aviinfo(InFileName);

width = FileInfo.Width;
height = FileInfo.Height;

blklist = zeros(height/blksize, width/blksize);
problist = zeros(height/blksize, width/blksize);

%-------- process frame by frame -------
%bck =frame2im(InSequence(start_frame_num));
if ( flag == 0 )
    im1_o = frame2im(InSequence(1));
    im2_o = frame2im(InSequence(last_index));
else
    im1_o = frame2im(InSequence(last_index));
    im2_o = frame2im(InSequence(1));
end


figure(100), imshow(im1_o), title('img');
figure(101), imshow(im2_o), title('img');
disp('func : motion segmentation_new');
%pause

% pre-align two images in frequecy domain 
if(size(im1_o,3)==3) 
    im1=rgb2gray(im1_o); 
    im2=rgb2gray(im2_o); 
else
    im1 = im1_o;
    im2 = im2_o;
end 

%total_level = 10;

% x = zeros(total_level);
% y = zeros(total_level);

xtran = 0;
ytran = 0;

level_index = 0;
bck_x_motion = 0;
bck_y_motion = 0;

%%%% subarna : new test
[delta, flag, bck_x_motion, bck_y_motion] = computedelta_new(im1, im2, width, height, blksize);

% here we need to calculate the edge mask (moving edge after camera motion detection 
edge_mask = var_moving_edge (InFileName, blksize, start_frame, end_frame, bck_x_motion, bck_y_motion);%%end_frame, %new_moving_edge

[blklist, problist, total_levels] = compute_delta_label_new(delta, im1, im2, blksize, edge_mask, flag, bck_x_motion, bck_y_motion);
%total_levels

max_val = max(max(problist));
if ( 0 ~= max_val )
    problist = problist/max_val;
end

%pause

%%%%% color segments %%%%%
%figure(50);imshow(rgb2gray(im1_o));title('orig image');%hold on;

%figure(4); imshow(blklist, []); title('segmented blocks')
%figure(5); imshow(problist, []); title('relative motion feature value of blocks')

return;
