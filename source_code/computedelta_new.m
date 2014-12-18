function [delta, flag, bck_x_motion, bck_y_motion] = computedelta_new(im1, im2, width, height, blksize)

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

if ( flag == 0 )
    maxvalue = max(max(delta));
    [ytran,xtran] = find( delta == maxvalue); 
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

bck_x_motion = xtran;
bck_y_motion = ytran;