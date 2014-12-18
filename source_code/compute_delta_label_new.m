function [blklist, problist, total_levels] = compute_delta_label_new(delta, im1, im2, blksize, edge_mask, flag, bck_x_motion, bck_y_motion)

EDGE_THRESOLD = 0.6;

[height width] = size(delta);
blklist = zeros(height/blksize, width/blksize);
problist = zeros(height/blksize, width/blksize);
im3 = zeros(height, width);
im4 = zeros(height, width);

w = width/blksize;
h = height/blksize;

range = 1:blksize;
black_blk = zeros(blksize, blksize);

tot = w*h;
heter = 0;


for i = 1:h
    for j = 1: w        
        yrange = (i-1)*blksize+1 : i*blksize;
        xrange = (j-1)*blksize+1 : j*blksize;        
        indices = j+(i-1)*w;
        
        if ( j == 1 )
            adjBlocks = [ indices-w+1;    % NE corner
                indices-w;      % top
                indices-1; ];   % left
        else
             adjBlocks = [ indices-w-1;    % NW corner
                indices-w;      % top
                indices-1; ];   % left
        end
        % mark the invalid indices    
        if(j == 1) adjBlocks(3) = 0; end   % first column means no left, NW
        if(i == 1) adjBlocks(2) = 0; end   % first row means no top, NE/NW 
              
        % pick out only the valid indices
        validOnes = find((adjBlocks<(h*w)) & (adjBlocks>0));
        adjGroups = adjBlocks(validOnes);
        
        need_check = check_motion_needed(i,j,edge_mask,EDGE_THRESOLD,indices,adjGroups);
        delta = double(delta);
                                  
        %%% add loop for mv finding here for edge/non-uniform block only
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ( need_check == 1 )
            [xtran, ytran, level, total_levels] = calc_best_warp(delta, im1, im2, bck_x_motion, bck_y_motion, xrange, yrange, blksize);
            blklist(i,j) = (level*10-1);  
            motion_value = (abs(xtran - bck_x_motion)/width) + (abs(ytran - bck_y_motion)/height);
            problist(i,j) =  motion_value;
            heter = heter + 1;
        else        
            r_u = i;
            c_u = j;
        end
              
    end
end

figure(11), imshow(edge_mask, []), title('edge mask');
[r1 c1] = find(edge_mask>EDGE_THRESOLD);
edge_blocks = size(r1, 1)
tot_blocks = tot;
non_uniform_block = heter;
end



function need_check = check_motion_needed(i,j,edge_mask,EDGE_THRESOLD,indices,adjGroups)
need_check = 0;
      if ( (edge_mask(i,j) >= EDGE_THRESOLD) )
          %%% to remove noise, take edge blocks whose any of the direct
          %%% neighbour is also edge block
          state = spurious_edge(i,j,edge_mask, EDGE_THRESOLD);
          if ( 0 == state)
            need_check = 1;
          end
      end       
      return;
end      
 

function state = spurious_edge(i,j,edge_mask, EDGE_THRESOLD)
    state = 1;
    top = 0;
    left= 0;
    bottom=0;
    right=0;
    NW=0;
    NE=0;
    SW=0;
    SE=0;

    [height width]=size(edge_mask);

    if ((i-1 >= 1) & (edge_mask(i-1,j) > EDGE_THRESOLD)) top = 1;end
    if ((j-1 >= 1) & (edge_mask(i,j-1) > EDGE_THRESOLD)) left = 1; end
    if ((i+1 <= height) & (edge_mask(i+1,j) > EDGE_THRESOLD)) bottom = 1; end
    if ((j+1 <= width) & (edge_mask(i,j+1) > EDGE_THRESOLD)) right = 1; end
    if ((i-1 >0 & j-1 >0) & (edge_mask(i-1,j-1) > EDGE_THRESOLD)) NW = 1; end
    if ((i-1 >0 & j+1 <= width) & (edge_mask(i-1,j+1) > EDGE_THRESOLD)) NE = 1; end
    if ((i+1 <= height & j+1 <= width) & (edge_mask(i+1,j+1) > EDGE_THRESOLD)) SE = 1; end
    if ((i+1 <= height & j-1 >0) & (edge_mask(i+1,j-1) > EDGE_THRESOLD)) SW = 1; end

 
    if (( top == 1)|(left == 1)|(bottom == 1)|(right == 1)|(NW == 1)|(NE == 1)|(SE == 1)|(SW == 1))        
        state = 0;
    end

end



 
function [xtran, ytran, level, total_levels] = calc_best_warp(delta, im1, im2, bck_x_motion, bck_y_motion, xrange, yrange, blksize)

[height width] = size(delta);
maxvalue = max(max(delta));
[ytran,xtran] = find( delta == maxvalue); 
level_index = 0;

minvar = 9999;
new_xtran = 0;
new_ytran = 0;

% allowable maximum segments
total_levels = 4;

level_index = 1;

level = 1;
%for level_index = 1: total_levels
while (1)
    maxvalue = max(max(delta)); %%% level_index-th highest peak value
    
    if ( level_index > 1 )
        if ( maxvalue < 0.05 ) % 0.03
            break
        end
    end
    
    [ytran,xtran]=find(delta == maxvalue);
       
    old_ytran = ytran;
    old_xtran = xtran;

    %%%% care for negative motion vector
    if ( ytran > height/2) ytran = -(height - ytran + 1);  end
    if ( xtran > width/2)  xtran = -(width - xtran + 1);   end
    
    if ( xtran > 0 )  xtran = xtran - 1;   end
    if ( ytran > 0 )  ytran = ytran - 1;   end
    
    %%%%% new add : subarna
    new_yindex = yrange + ytran;
    new_xindex = xrange + xtran;
    
    %%%%% masking out-of-range values   
    invalidOnes = find((new_yindex> height));  new_yindex(invalidOnes) = height;
    invalidOnes = find((new_yindex < 1));  new_yindex(invalidOnes) = 1;
    
    invalidOnes = find((new_xindex> width));  new_xindex(invalidOnes) = width;
    invalidOnes = find((new_xindex < 1));  new_xindex(invalidOnes) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im3 = im1(new_yindex, new_xindex);
    extracted_block = im2(yrange, xrange) - im3;
        
    var = sum(sum(abs(extracted_block)));
    
    if ( var < minvar )
        minvar = var;
        new_xtran = xtran;
        new_ytran = ytran;
        level = level_index;
    end
    
%     if ( (var/(blksize*blksize)) < 2 * level_index ) %2*level
%         break; 
%     end
    
    delta(old_ytran, old_xtran) = 0;
    level_index = level_index + 1;
 end

 xtran = new_xtran;
 ytran = new_ytran;
 total_levels = level_index;
end

      
      
      