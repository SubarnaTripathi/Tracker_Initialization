function [delta, im1, im2, xtran, ytran, blklist, problist, bck_x_motion, bck_y_motion] = computedelta(im1, im2, width, height, level, color,blksize, blklist, problist, xtran_old, ytran_old, bck_x_motion, bck_y_motion)

flag = 0;
zero_pad = zeros(height, width);

im1=double(im1); 
im2=double(im2);

F1 = fft2(im1); 
F2 = fft2(im2); 
fz=F1.*conj(F2); 
fm=abs(F1.*F2); 

if ( fm ~= 0 )
    div=fz./fm;
    delta = ifft2(div);
else
    flag = 1;
    delta = zero_pad;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(1), mesh(1:width,1:height,abs(delta))
%xlabel('F_x'), ylabel('F_y'), zlabel('Magnitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( flag == 0 )
    maxvalue = max(max(delta));
    %delta=fliplr(flipud(delta)); 
    [ytran,xtran] = find( delta == maxvalue); 

    for level_index = 2: level
        delta(ytran, xtran) = 0;
        maxvalue = max(max(delta));
        [ytran,xtran]=find(delta == maxvalue);
    end

else %flag    
    ytran = 1;
    xtran = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% care for negative motion vector

if ( ytran > height/2)
    ytran = -(height - ytran);
end
if ( xtran > width/2)
    xtran = -(width - xtran);
end

%sprintf('xtran = %d, ytran = %d', xtran, ytran);

if ( level_index == 1 )
    bck_x_motion = xtran;
    bck_y_motion = ytran;
end

%fprintf('xtrans :  %d \n ytrans :  %d \n',xtran,ytran);

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
%im4(1:height, 1:width) = [im1((uint16(yindex1 - ytran/2)), (uint16(xindex1 - xtran/2)))];
%im4 = im4 - double(im2((uint16(yindex2 + ytran/2)), (uint16(xindex2 + xtran/2)))); % instead of taking average, take difference

%im4 = im4./2;

% figure(2), imshow(im1), title('im1');
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
        if ( (var/(blksize*blksize)) < 3 * level ) %2*level
            if ( blklist(i,j) == 0 )
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %draw block boundaries
%                 line([(j-1)*blksize j*blksize],[(i-1)*blksize (i-1)*blksize],'color',color); %% line(X,Y) 
%                 hold on
%                 line([(j-1)*blksize j*blksize],[i*blksize i*blksize],'color',color); %% line(X,Y)
%                 hold on
%                 line([(j-1)*blksize (j-1)*blksize],[(i-1)*blksize i*blksize],'color',color); %% line(X,Y)
%                 hold on
%                 line([j*blksize j*blksize],[(i-1)*blksize i*blksize],'color',color); %% line(X,Y)
%                 hold on
                %%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%% prepare input for next dominant motion
                im1(yrange, xrange) = black_blk;
                im2(yrange, xrange) = black_blk;
                
                blklist(i,j) = (level-1);
                
                problist(i,j) = abs(xtran - bck_x_motion)/width + abs(ytran - bck_y_motion)/height;
            end
            
            %blklist(i,j) = level;
        end
    end
end

end