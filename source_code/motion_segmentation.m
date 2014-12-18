% PROJ : depth from motion
% simple motion estimation/segmentation by phase correlation method
% Author : Subarna Tripathi
% Date : March 2010
% example use : 
% [blklist, problist] = motion_segmentation ('test.avi', 4, 3, 4);

function [blklist, problist] = motion_segmentation (InFileName, blksize, start_frame, end_frame)

InFileName = 'akio_cif.avi'; %'films_CIF_part2_fr16.avi'; %'akio_cif.avi'; %'foreman_cif.avi'; %'renata_qcif.avi'; %'racingcar.avi'; 
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
while (1)
    level_index = level_index+1;    
    %for level_index = 1: total_level
    %%%%% input for this level of motion
    %figure(2), imshow(uint8(im1)), title('im1');
    %figure(3), imshow(uint8(im2)), title('im2');
        
    if (level_index == 1 )
        color = [1 0 0];
    elseif (level_index == 2 )
        color = [0 1 0];
    elseif (level_index == 3 )
        color = [1 1 0];
    elseif (level_index == 4 )
        color = [0 1 1];
    elseif (level_index == 5 )
        color = [1 0 1];
    elseif (level_index == 6 )
        color = [0 0 0];
    else
        color = [1 1 1];
    end
    
    %[delta, im1, im2, xtran, ytran, blklist, problist, bck_x_motion, bck_y_motion] = computedelta(im1, im2, width, height, level_index, color,blksize, blklist, problist, xtran, ytran, bck_x_motion, bck_y_motion);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ( level_index == 1 )
       [delta, im1, im2, xtran, ytran, blklist, problist, bck_x_motion, bck_y_motion] = computedelta(im1, im2, width, height, level_index, color,blksize, blklist, problist, xtran, ytran, bck_x_motion, bck_y_motion);
    end
    maxvalue = max(max(delta));
    [ytran,xtran] = find( delta == maxvalue); 
    for temp_index = 2 : level_index 
       delta(ytran, xtran) = 0;
       maxvalue = max(max(delta));
       [ytran,xtran]=find(delta == maxvalue);
    end
          
    %%%% care for negative motion vector
    if ( ytran > height/2)
        ytran = -(height - ytran);
    end
    if ( xtran > width/2)
        xtran = -(width - xtran);
    end
    
    %%%% for debugging
    %s = sprintf('level = %d  [xtran, ytran]= [%d, %d]',level_index, xtran, ytran)
     
    [blklist, problist] = compute_delta_label(im1, im2, width, height, level_index, color,blksize, blklist, problist, xtran, ytran, bck_x_motion, bck_y_motion);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %x(level_index) = xtran;
    %y(level_index) = ytran;
    
    [numy, numx] = find(blklist);
    per = size(numy);
    %%%% most of the pixels covered
    if( per(1) > 0.99*(height*width/(blksize*blksize)) )
        %level_index
        %pause
        break;
    end
end

max_val = max(max(problist));
if ( 0 ~= max_val )
    problist = problist/max_val;
end

total_level = level_index

%%%%% color segments %%%%%
figure(10);imshow(rgb2gray(im1_o));title('orig image');hold on;

horz_blk = width/blksize;
vert_blk = height/blksize;
for i = 1:vert_blk
    for j = 1: horz_blk
        
        yrange = (i-1)*blksize+1 : i*blksize;
        xrange = (j-1)*blksize+1 : j*blksize;
            
%         if ( blklist(i,j) == 1 )
%             color = [1 0 0];
%         elseif (blklist(i,j) == 2 )
%             color = [0 1 0];
%         elseif (blklist(i,j) == 3 )
%             color = [1 1 0];
%         elseif (blklist(i,j) == 4 )
%             color = [0 1 1];
%         elseif (blklist(i,j) == 5 )
%             color = [1 0 1];
%         elseif (blklist(i,j) == 6 )
%             color = [0 0 0];
%         else
%             color = [1 1 1];  
%         end    

        if ( blklist(i,j) > 1 ) 
            color = [1/blklist(i,j)  1/blklist(i,j)  1/blklist(i,j) ];
        else
            color = [1 1 1];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %draw block boundaries
%         line([(j-1)*blksize j*blksize],[(i-1)*blksize (i-1)*blksize],'color',color); %% line(X,Y) 
%         hold on
%         line([(j-1)*blksize j*blksize],[i*blksize i*blksize],'color',color); %% line(X,Y)
%         hold on
%         line([(j-1)*blksize (j-1)*blksize],[(i-1)*blksize i*blksize],'color',color); %% line(X,Y)
%         hold on
%         line([j*blksize j*blksize],[(i-1)*blksize i*blksize],'color',color); %% line(X,Y)
%         hold on
        %%%%%%%%%%%%%%%%%%%%%%%%%%           
        end
end


%end


figure(2); imshow(blklist, []); title('segmented blocks')
figure(3); imshow(problist, []); title('relative motion feature value of blocks')


return;
