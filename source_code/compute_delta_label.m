function [blklist, problist] = compute_delta_label(im1, im2, width, height, level, color, blksize, blklist, problist, xtran, ytran, bck_x_motion, bck_y_motion)

if ( level == 1 )
    bck_x_motion = xtran;
    bck_y_motion = ytran;
end

im3 = zeros(height, width);
im4 = zeros(height, width);

if ( xtran > 0 )
    xtran = xtran - 1;
elseif (xtran < 0 )
    xtran = xtran+1;
end

if ( ytran > 0 )
    ytran = ytran - 1;
elseif ( ytran < 0 )
    ytran = ytran + 1;
end

xindex1 = xtran + 1: width+xtran;
yindex1 = ytran + 1: height+ytran;

new_yindex = yindex1 - ytran;
new_xindex = xindex1 - xtran;


%%%%% masking negative indices
for i = 1:height
    if(new_yindex(i) > height)
        new_yindex(i) = height;
    end
    if(new_yindex(i) < 1)
        new_yindex(i) = 1;
    end
end

for i = 1:width
    if(new_xindex(i) > width)
        new_xindex(i) = width;
    end
    if(new_xindex(i) < 1)
        new_xindex(i) = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

im3(1:height, 1:width ) = im1(new_yindex, new_xindex);
%figure(4), imshow(uint8(im3)), title('im1 with -xtrans, -ytrans ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xindex2 = 1: width+xtran;
% yindex2 = ytran + 1: height+ytran;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xindex1 = xtran/2 + 1: width+xtran/2;
yindex1 = ytran/2 + 1: height+ytran/2;

xindex2 = -ytran/2 + 1: width-xtran/2;
yindex2 = -xtran/2 + 1: height-ytran/2;

im4 = im2 - im3;

%figure(2), imshow(im1), title('im1');hold on
% figure(3), imshow(im2), title('im2');
%figure(5), imshow(uint8(im4)), title('intermediate representation');hold on


horz_blk = width/blksize;
vert_blk = height/blksize;


range = 1:blksize;
black_blk = zeros(blksize, blksize);

for i = 1:vert_blk
    for j = 1: horz_blk
        
        yrange = (i-1)*blksize+1 : i*blksize;
        xrange = (j-1)*blksize+1 : j*blksize;
      
        extracted_block = im4(yrange, xrange);
        
        var = sum(sum(abs(extracted_block)));
        if ( (var/(blksize*blksize)) < 2 * level ) %2*level
           
            if ( blklist(i,j) == 0 )
                %%%% prepare input for next dominant motion
                im1(yrange, xrange) = black_blk;
                im2(yrange, xrange) = black_blk;
                
                blklist(i,j) = (level-1);
                
                problist(i,j) = (abs(xtran - bck_x_motion)/width) + (abs(ytran - bck_y_motion)/height);
            end
            
            %blklist(i,j) = level;
        end
    end
end


end