function color_mask_new = merge_clusters(im, color_mask)

tot_num = size(color_mask,1);
for obj_num = 1: tot_num 
    %obj_num

    [I1 J1] = find(color_mask == obj_num); 
    
    if (obj_num < size(color_mask,1) )
        [I2 J2] = find(color_mask == obj_num+1); 
        
        if ( (abs(max(max(I1)) - max(max(I2))) < 10) &
            (abs(max(max(J1)) - max(max(J2))) < 10 ) &
            (abs(min(min(I1)) - min(min(I2))) < 10) &
            (abs(min(min(J1)) - min(min(J2))) < 10) )
            
            if ( abs(mean(im(I1, J1, 1)) - mean(im(I2,J2,1))) < 10 &
                abs(mean(im(I1, J1, 2)) - mean(im(I2,J2,2))) < 10
                abs(mean(im(I1, J1, 3)) - mean(im(I2,J2,3))) < 10 )
              
                %%% merge
                color_mask(I1, J1) = obj_num+1;
                color_mask(I2, J2) = obj_num+1;
                %color_mask = reshape(color_mask,size(color_mask));                
            end     
         end        
    end
end


U = unique(color_mask_new);
for i = 1:size(U,1)
      color_mask_new(color_mask == U(i)) = i;
      color_mask_new = reshape(color_mask_new,size(color_mask));
end