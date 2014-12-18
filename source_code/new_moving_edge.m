% Subarna Tripathi
% determines edge feature of blocks, 
% two cases, if moving(connected comp) and not-moving (hough)

function blk_texture_list = new_moving_edge (InFileName, blksize, start_frame, end_frame)


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
FileInfo=aviinfo(InFileName);
width = FileInfo.Width;
height = FileInfo.Height;

blk_texture_list = zeros(height/blksize, width/blksize);

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

% test = edge(im1, 'canny'); %'robert'
% figure(5); imshow(test, []), title('canny');
% pause

Z = imsubtract(im1,im2);
%Z = double(Z);

horz_blk = width/blksize;
vert_blk = height/blksize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check for stationary image pairs
stationarity_pair = 0;

[R C] = find(Z > 50);
non_stationarity_size = size(R);
if ( non_stationarity_size(1) < 50 )
    %%%%% make Z in different way
    stationarity_pair = 1;
    Z = im1;    
    Z = medfilt2(Z,[3 3]);
    %%%% replace one edge detection to 4 edge detection responses %%%%
    P = edge(Z, 'robert');
    %P = edge(Z, 'canny');
    Z = P;
    %Z = double(Z);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%% replace one edge detection to 4 edge detection responses %%%%
%P = edge(Z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);imshow(Z);title('image'), hold on

color = [1 0 0];
for i = 1:vert_blk
    for j = 1: horz_blk        
        yrange = (i-1)*blksize+1 : i*blksize;
        xrange = (j-1)*blksize+1 : j*blksize;       
        X  = Z(yrange, xrange); 
        
        %%%% replace one edge detection to 4 edge detection responses %%%%           
        if ( stationarity_pair == 0 )
             BW = (X >= 40 );
            %BW = (new_Y <= 10 & new_Y > 1 );
            [L, num] = bwlabeln(BW,8);

            edge_dir = num;      
            for class_num = 1:num
                [r,c] = find(L == class_num );
                a = size(r, 1);
                sum1 = sum(sum(X(r,c)))/(a);            

                if ( a < 4 | (sum1 < 50 )) 
                    edge_dir = edge_dir-1;
                end           
            end
            
            if (edge_dir == 0 )
                blk_texture_list(i,j) = 0.1;
            elseif (edge_dir == 1 )
                blk_texture_list(i,j) = 0.8;
            elseif (edge_dir == 2 )
                blk_texture_list(i,j) = 0.4;
            elseif (edge_dir == 3 )
                blk_texture_list(i,j) = 0.25; 
            else
                %edge_dir
                %pause
                blk_texture_list(i,j) = 1/edge_dir;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            [I J] = find(X > 0);
            edge_dir = hough1(X, blksize);
            
            %%%%%%%%%%%%%%%%%% edge strength
            %local_edge_pixel = size(I);
            %edge_pixel = local_edge_pixel(1);
            %blk_texture_list(i,j) = 2*edge_pixel/(blksize*blksize); %0.1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            if (edge_dir == 0 )
                blk_texture_list(i,j)= 0;
            elseif (edge_dir == 1) 
                blk_texture_list(i,j) = 1;
            elseif (edge_dir == 2) 
                blk_texture_list(i,j) = 0.8;
            elseif (edge_dir == 3) 
                blk_texture_list(i,j) = 0.6;
            elseif (edge_dir == 4) 
                blk_texture_list(i,j) = 0.4;
            else
                blk_texture_list(i,j) = 1/edge_dir;   
            end
        end        
    end
end

%%%% post-processing block texture list for smoothing %%%%%
%blk_texture_list = imfilter(blk_texture_list, [ 0 1 0; 1 2 1; 0 1 0],'replicate');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





