function [obj_map,level, problist_pp] = assign_obj_label(color_mask, problist_orig, edge_mask_orig, blksize, rgb_image)

problist = problist_orig;

% resize to frame resolution
problist = imresize(problist_orig, blksize );
edge_mask = imresize(edge_mask_orig, blksize);
obj_map = zeros(size(color_mask));
problist_pp = zeros(size(color_mask));


obj_map_temp = obj_map;
obj_map_temp1 = obj_map_temp;

raw_max_motion_val = max(max(problist));

U = unique(color_mask);
num_obj_layer = size(U,1);

level = 0;
for i = 1:num_obj_layer
      [I J] = find(color_mask == U(i));
      
      raw_h = (max(I) - min(I));
      raw_w = (max(J) - min(J));
      
      region_size = size(I,1);
      max_motion_val = max(problist(color_mask == U(i)));
     
      temp_arr = problist(color_mask == U(i));
      [i1 j1 v] = find(temp_arr);      
      median_motion_val = median(v); % median of non-zero entries  
           
      %median_motion_val = (median_motion_val/raw_max_motion_val);  %%%%
      %normalize value, no need max is always 1
            
      m_threshold1 = 0.2;  %0.3, % 0.6
      m_threshold2 = 0.8; %0.6;  % 0.5, %0.3, % 0.6
      m_threshold3 = 0.9;
      e_threshold1 = 0.2; % 0.4
      e_threshold2 = 0.7;
      e_threshold3 = 0.003;
      e_threshold4 = 0.25;
      
      %[I1 J1] = find(problist(color_mask == U(i)) == max_motion_val);
      [I1 J1] = find(problist(color_mask == U(i)));
      motion_present_fraction = size(I1,1)/region_size;
      [I2 J2] = find(edge_mask(color_mask == U(i)));
      edge_present_fraction = size(I2,1)/region_size;    
         
      region_texture = 0;
      region_texture = region_is_common_texture(color_mask, U, i, rgb_image, median_motion_val, edge_present_fraction, raw_h, raw_w);  % 0 = no common texture, 1 = sky, 2 = tree leaves/grass
           
      %if ( ((edge_present_fraction >= e_threshold1)&(edge_present_fraction <= e_threshold2 )) | (3 == region_texture) | ( (median_motion_val >= m_threshold1) & (median_motion_val <= m_threshold2)) | (0 == region_texture) )
      if ( (3 == region_texture) | ( (median_motion_val >= m_threshold1) & (median_motion_val <= m_threshold2) ) | (0 == region_texture) | (2 == region_texture) | (median_motion_val >= m_threshold3) )
          if ( ((0 == region_texture) & ( median_motion_val >= m_threshold1) ) | (3 == region_texture) | (((2 == region_texture) |(1 == region_texture )) & (median_motion_val >= m_threshold3) & (edge_present_fraction < e_threshold3)) | ((median_motion_val >= m_threshold3) & (edge_present_fraction < e_threshold4) ) )              
               %if (3 == region_texture)
                level = level+1; 
                %%% offset needs to be added to nutralize quantization effect 
                obj_map(color_mask == U(i)) = level; %max_motion_val*level;  %median_motion_val
                
                problist_pp(color_mask == U(i)) = median_motion_val; %max_motion_val*level;  %
               %end
          end          
      end
           
      obj_map = reshape(obj_map,size(color_mask));        
      problist_pp = reshape(problist_pp,size(color_mask));   
      
%       if ( 3 ~= region_texture ) 
%           obj_map_temp = obj_map_temp1;    
%           obj_map_temp(color_mask == U(i))= 1;
%           obj_map_temp = reshape(obj_map_temp,size(color_mask));   
%           figure(15), imshow(obj_map_temp,[]);
%           disp('press to see the stat')
%           pause
%           disp(sprintf('region_texture = %d, motion_val = %f, edge_fraction=%f,press to see region %d',region_texture, median_motion_val, edge_present_fraction,i));
%       end
end


%U = unique(obj_map);
%obj_map = U;

end



%%%%%% dteects most common texture based on color analysis
%%%%% discards those regions from possible objects identified
% 0 = no common texture, 1 = sky, 2 = tree leaves/grass, 3 = discard since
% small
function region = region_is_common_texture(color_mask, U, i, rgb_image, median_motion_val, edge_present_fraction, raw_h, raw_w)

    debug = 1;

    region = 0;    
    min_size_thresh = 4;    
    
    %%%%% discrading at this early stage, needs to be verified
    if ( (raw_h < min_size_thresh) | (raw_w < min_size_thresh) )
        %region = 1;
    end
    
    %%%%% discrading at this early stage, needs to be verified
    if ( (raw_h >= (size(color_mask, 1)-5)) | (raw_w >= ( size(color_mask,2) - 5)) )
        region = 1;
    end
    
    copy = rgb_image;
    Chrom = rgb2ycbcr(copy);
    hsvimage = rgb2hsv(copy);

      
    object = 0;
    sky = 1; %%%% may be also water
    leaves = 2; %%%% may be also water and grass
    skin = 3;
    edges = 4;
    
    r_image = rgb_image(:,:,1);
    g_image = rgb_image(:,:,2);
    b_image = rgb_image(:,:,3);
    
    r_image = r_image(:);
    g_image = g_image(:);
    b_image = b_image(:);
    
    
    y_image = Chrom(:,:,1);
    cb_image = Chrom(:,:,2);
    cr_image = Chrom(:,:,3);
    
    y_image = y_image(:);
    cb_image = cb_image(:);
    cr_image = cr_image(:);
    
    %%%% check for sky
    [R C] = find(color_mask == U(i));
      
    median_r = median(r_image(color_mask == U(i))); %median
    median_g = median(g_image(color_mask == U(i))); %median
    median_b = median(b_image(color_mask == U(i))); %median
    
    median_y = median(y_image(color_mask == U(i))); %median
    median_cb = median(cb_image(color_mask == U(i))); %median
    median_cr = median(cr_image(color_mask == U(i))); %median
    
    median_cg = 255 - median_cb/2 - median_cr/2;
    
    avg_y = double(median_r/3) + double(median_g/3) + double(median_b/3);
    
    %avg_y = ((median_r + median_g + median_b))/3;
      
        
%     disp(sprintf('obj_layer =%d, median_r =%d, median_g = %d, median_b = %d', i, median_r,median_b,median_b));
%     disp(sprintf('median_y =%d, median_cb = %d, median_cr = %d',median_y,median_cb,median_cr));
       
    %(median_cb<113) %(median_y<117)  %(110<= median_y) %(median_cr>= 128)
    if ( (120<= median_y) & (median_y<= 130) & (110<= median_cb) & (median_cb<= 120) & (median_cr>= 135) & (median_cr <= 150) )   
        if (1 == debug)
            %disp(sprintf('median_y =%d, median_cb = %d, median_cr = %d',median_y,median_cb,median_cr));
        end
        region = skin;
        
    elseif( (130<= median_r) & (median_r<= 140) & (130<= median_g) & (median_g<= 140) & (90 <= median_g) & (median_b <= 95) )   
        if (1 == debug)
           % disp(sprintf('median_r =%d, median_g = %d, median_b = %d',median_r,median_g,median_b));
        end
        region = skin;
       
    elseif ( ((median_b >= median_g + 30) & (median_g >= median_r)) | ((median_b >= median_g) & (median_g >= median_r)) | ((median_b >= median_r) & (median_r >= median_g)) | ((avg_y >= 200 ) & (abs(median_b - median_g)<=30) & (abs(median_r - median_g)<=30)) )          
        region = sky;
                     
    elseif ( (median_g >= (median_r-5)) & (median_g >= (median_b-5)) )%& (edge_present_fraction > 0.1))  % add offset 20
        region = leaves;
    end
    
    
    if( (edge_present_fraction > 0.4) & (leaves ~= region) )
        region = 0;
    end
    
    if( (median_motion_val > 0.9) & (edge_present_fraction > 0.2) )
        region = edges;
    elseif( (median_motion_val > 0.70)  & (edge_present_fraction > 0.2) & (leaves ~= region) )
        region = 0;
    end
    
%     if (region == 0)
%         disp(sprintf('median_y = %f, median_cb = %f, median_cr = %f, median_cg = %f', median_y, median_cb, median_cr, median_cg));
%         disp(sprintf('avg_y = %f, median_r = %f, median_g = %f, median_b = %f', avg_y, median_r, median_g, median_b));
%     end
   
%     if ( 3 == region ) 
%           obj_map_temp = zeros(size(color_mask));
%           obj_map_temp(color_mask == U(i))= 1;
%           obj_map_temp = reshape(obj_map_temp,size(color_mask));   
%           figure(15), imshow(obj_map_temp,[]);
%           disp('press to see the stat')
%           pause
%           disp(sprintf('median_r =%d, median_g =%d, median_b =%d, median_y =%d, median_cb =%d, median_cr =%d',median_r,median_g,median_b,median_y,median_cb,median_cr));
%     end

end