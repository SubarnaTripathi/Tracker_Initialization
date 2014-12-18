%%%% merge possible clusters
function [obj_map_new, num_obj,k, M] = merge_clusters_new(im, obj_map, obj_map_orig, color_mask, problist_orig, problist_pp, blksize)

    %mov = avifile('example.avi')
    debug = 0;
    merging_dec = 0;
    no_more_change = 0;

    obj_map_new = obj_map;
    prob_list = imresize(problist_pp, blksize );
    bw = zeros(size(obj_map));
    bw1 = bw;

    k = 1;
    
    iter = 0;
    while(1)
        obj_map_new = obj_map;
        iter = iter + 1;
        if (no_more_change == 1)
            break;
        end
    
        U = unique(obj_map);
        num_obj = size(U,1);
              
        
        no_more_change = 1;
        level = 1;
        for i = 1:num_obj
            for j = 1:num_obj
                
                obj_map_new = obj_map;

                if ( i == j )
                    continue;
                end

                if (obj_map(obj_map == U(i))== 0 )
                    %disp('i th region already bck');
                    continue;
                end

                if (obj_map(obj_map == U(j))== 0 )
                    %disp('j th region already bck');
                    continue;
                end
                
                %%%%% no merging if any one region if too concave
                bw = bw1;
                bw(obj_map == U(i))= 1;                 %bw = im2bw(bw);
                [solidity1] = check_ind_solidity(bw);
                bw = bw1;
                bw(obj_map == U(j))= 1;                 %bw = im2bw(bw);
                [solidity2] = check_ind_solidity(bw);
                local_solidity_thresh = 0.05;
                               
                if ( (solidity1 < local_solidity_thresh ) | (solidity2 < local_solidity_thresh ) )
                    %disp('too concave to join');
                    continue;
                end
                

                merging_dec = merge(im, problist_orig, problist_pp, obj_map, obj_map_orig, U,i,j);
                
                if ( 1 == debug )
                    disp(sprintf('i = %d, j= %d merging_dec = %d', i,j, merging_dec));
                end
                
                if ( 1 == merging_dec ) 
                    obj_map_new = obj_map;
                    
%                     bw3 = zeros(size(obj_map));                 
%                     bw3(obj_map == U(i))= 0.4;
%                     bw3(obj_map == U(j))= 0.8;
%                     figure(50); imshow(bw3), title('ui <=> uj');
%                     pause
                    
                    %disp(sprintf('%d th iteration',iter));
                    %disp(sprintf('merging happenning between region %d and %d', i, j));
                    %pause
                    obj_map_new(obj_map == U(i)) = U(i); %level;  %median
                    obj_map_new(obj_map == U(j)) = U(i); %level;  %median
                    
                        
                    no_more_change = 0;
                    level = level+1;
                    
                    disp('within function merge_clusters_new');
                    %pause
                    
                end
                          
                obj_map = reshape(obj_map_new,size(obj_map));                 
                 
                %pause                
                figure(11),imshow(obj_map, []), title('in between merging'); 
                
                %plot_command
                F = getframe();
                M(k) = F;
                %mov = addframe(mov,F);
                
                k = k+1;                                
            end
            %obj_map = reshape(obj_map_new,size(obj_map));    
        end  
    end
    disp(sprintf('total number of iterations = %d',iter));    
    
    obj_map_new = obj_map; 
    
    
    movie(M)

end



   
          







function merging_dec = merge(im, problist_orig, problist_pp, obj_map, obj_map_orig, U, i,j) 

    debug = 1;
    
    connected = 0;
    
    %%%% test : debug
    obj_map_temp1 = zeros(size(obj_map));
    obj_map_temp1(obj_map == U(i))= 1;   
    obj_map_temp1(obj_map == U(j))= 2;
    obj_map_temp1 = reshape(obj_map_temp1,size(obj_map));
    %figure(16), imshow(obj_map_temp1,[]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    merging_dec = 0;

    obj_map_temp = obj_map(:); 
            
    
    [connected, overlap] = check_connectedness(obj_map, U, i,j );
          
    if (connected == 1 )                
        if ( 1 == debug )
            disp(sprintf('segments %d and %d are connected; overlap = %d\n', i,j, overlap)); 
%             bw2 = zeros(size(obj_map));
%             bw2(obj_map == U(i))= 0.4;
%             bw2(obj_map == U(j))= 0.8;     
%             figure(50); imshow(bw2), title('ui <=> uj');
%             pause
        end
        
        %%%%%% merging predicates of evidence of boundary %%%%%%%%
        [intc1, bw1] = calc_internal_diff(im, problist_orig, problist_pp, obj_map, U, i );
        [intc2, bw2] = calc_internal_diff(im, problist_orig, problist_pp, obj_map, U, j ); 
        [Mintc1c2] = calc_min_internal_diff( im, problist_orig, problist_orig, obj_map, U, i, j, intc1, intc2);
        
        [diff, bw] = calc_pair_diff(im, problist_orig, problist_pp, obj_map, U, i, j);
        
        if ( diff <= Mintc1c2)
            merging_dec = 1;
            disp(sprintf('pair diff = %d, minIntdif = %d', diff, Mintc1c2));
        end
        
        %%%% test for solidity 
        [solidity] = check_merged_solidity(bw);
        pair_solidity_thresh = 0.1; %0.2
        if (solidity < pair_solidity_thresh)
            merging_dec = 0;
        end

        %figure(16), imshow(obj_map_temp1,[]);
        disp(sprintf('intc1 = %f, intc2 = %f, Mintc1c2 = %f, diff = %f, solidity = %f, merging_dec = %d', intc1, intc2, Mintc1c2, diff, solidity, merging_dec));
        %pause
    end% end of connectivity test
    
    
    %%%%%% one segment is contained within other
    [r1 c1] = find(obj_map == U(i));
    [r2 c2] = find(obj_map == U(j));
    
    margin = 10;

    if ( ((min(r1)-margin) < min(r2)) & ((min(c1)-margin) < min(c2)) & ((max(r1)+margin) > max(r2)) & ((max(c1)+margin) > max(c2)) )  %ith region contains jth region
        if (1 == debug)
            disp(sprintf('containment: j=%d in i=%d', j, i));
        end
        merging_dec = 1;
    elseif ( ((min(r2)-margin) < min(r1)) & ((min(c2)-margin) < min(c1)) & ((max(r2)+margin) > max(r1)) & ((max(c2)+margin) > max(c1)) )  %jth region contains jth region
        if (1 == debug)
            disp(sprintf('containment: i=%d in j=%d', i, j));
        end
        merging_dec = 1;
    end              
end


function [intc, bw] = calc_internal_diff(im, problist_orig, problist_pp, obj_map, U, index )

        problist_orig_temp = problist_orig(:);
        bw = zeros(size(obj_map));
        bw(obj_map == U(index))= 1;
        bw = im2bw(bw);
        
        bw_m = zeros(size(obj_map));
        bw_m(find(problist_orig))= 1;

        im_r = im(:,:,1);
        im_r = im_r(:);
        im_g = im(:,:,2);
        im_g = im_g(:);
        im_b = im(:,:,3);
        im_b = im_b(:);
                
        range_r = double(max(im_r(bw == 1)) - min(im_r(bw == 1)));
        range_g = double(max(im_g(bw == 1)) - min(im_g(bw == 1)));
        range_b = double(max(im_b(bw == 1)) - min(im_b(bw == 1)));
               
       % range_m = double(max(problist_orig_temp(bw == 1)) - problist_orig_temp(im_b(bw == 1)));
               
        cd_r = double(mean(im_r(bw == 1))); % median
        cd_g = double(mean(im_g(bw == 1))); % median
        cd_b = double(mean(im_b(bw == 1))); % median  
        
       % cd_m = double(mean(problist_orig_temp(bw == 1))); % median  
        
        diff_r = double(max(abs((im_r(bw == 1) - cd_r))));
        diff_g = double(max(abs((im_g(bw == 1) - cd_g))));
        diff_b = double(max(abs((im_b(bw == 1) - cd_b)))); 
        
        color_fact = 0.27; % should be between 0.1 to 0.3
        motion_fact = (1 - 3*color_fact);
        
%         color_fact = 1;
%         motion_fact = 2.5;
        
%         if (range_m ~= 0 )
%             motion_cost = (motion_fact*diff_m/range_m); 
%         else
%             motion_cost = motion_fact*diff_m; 
%         end
        
        %intc = (color_fact*diff_r/range_r + color_fact*diff_g/range_g  + color_fact*diff_b/range_b); % + motion_cost; 
              
        intc = (diff_r/(range_r*3)) + (diff_g/(range_g*3)) + (diff_b/(range_b*3)); 
             
end


function [pair_diff, bw] = calc_pair_diff(im, problist_orig, problist_pp, obj_map, U, i, j)

        obj_map_temp = obj_map(:);                
        problist_orig_temp = problist_orig(:); 
        problist_pp_temp = problist_pp(:); 
        
        bw = zeros(size(obj_map));
        bw(obj_map == U(i))= 1;
        bw(obj_map == U(j))= 1;           
        bw = im2bw(bw);
        
        im_r = im(:,:,1);
        im_r = im_r(:);
        im_g = im(:,:,2);
        im_g = im_g(:);
        im_b = im(:,:,3);
        im_b = im_b(:);
                
        range_r = double(max(im_r(bw == 1)) - min(im_r(bw == 1)));
        range_g = double(max(im_g(bw == 1)) - min(im_g(bw == 1)));
        range_b = double(max(im_b(bw == 1)) - min(im_b(bw == 1)));
        
        %range_m = double(min(obj_map(bw == 1)));  % max        
        %range_m = double(max(problist_pp_temp(bw == 1)) - min(problist_pp_temp(bw == 1)));
           
        cd_r = double(mean(im_r(bw == 1))); % median
        cd_g = double(mean(im_g(bw == 1))); % median
        cd_b = double(mean(im_b(bw == 1))); % median        
        
        diff_r = double(max(abs((im_r(bw == 1) - cd_r))));
        diff_g = double(max(abs((im_g(bw == 1) - cd_g))));
        diff_b = double(max(abs((im_b(bw == 1) - cd_b))));
        
        %diff_m = abs(max(obj_map(obj_map_temp == U(i))) - max(obj_map(obj_map_temp == U(j))));
        %diff_m = abs(max(problist_pp(obj_map_temp == U(i))) - max(problist_pp(obj_map_temp == U(j))));
        
        color_fact = 0.27; % should be between 0.1 to 0.3
        motion_fact = (1 - 3*color_fact);
        
%         color_fact = 1;
%         motion_fact = 2.5;
%         
%         if (range_m ~= 0 )
%             motion_cost = (motion_fact*diff_m/range_m); 
%         else
%             motion_cost = motion_fact*diff_m; 
%         end
        
        pair_diff = (color_fact*diff_r/range_r + color_fact*diff_g/range_g  + color_fact*diff_b/range_b); % + motion_cost;
        
end


function [Mintc1c2] = calc_min_internal_diff( im, problist_orig, problist_pp, obj_map, U, i, j, intc1, intc2)

    obj_map_temp = obj_map(:); 
    k =  500;   %constant
    
    [height width] = size(obj_map);
    area = double(height*width);
        
    bw = zeros(size(obj_map));
    bw(obj_map == U(i))= 1;
    [R C] = find(bw);
    region_size1 = size(R,1);    
    area_per1 = double(100*region_size1/area);  
    %k = intc1;
    t1 = double(k/area_per1);
    
    bw(obj_map == U(i))= 0;
    bw(obj_map == U(j))= 1;
    [R C] = find(bw);
    region_size2 = size(R,1);
    area_per2 = double(100*region_size2/area);
    %k = intc2;
    t2 = double(k/area_per2);
    
    var1 = intc1 + t1;
    var2 = intc2 + t2;
    
    
    % Mintc1c2 is the minimum of var1 and var2
    Mintc1c2 = var1;
    if ( var2 < var1)
        Mintc1c2 = var2;
    end
end


function [solidity] = check_merged_solidity(bw)

 BW = double(bw);

 stats = regionprops(BW, 'Solidity');
 solidity = stats(1).Solidity; 

end



function [t_solidity] = check_ind_solidity(bw)

    [r c v] = find(bw);
    actual_area = double(size(v,1));

    if ( actual_area < 1)
        t_solidity = 0.001;
        return;
    end

    t_height = double(max(r) - min(r));
    t_width = double(max(c) - min(c));

    t_area = (t_width)*(t_height);

    t_solidity = (actual_area/t_area);

end



function [connected, overlap] = check_connectedness(obj_map, U, i,j )

    connected = 0;
    overlap = 0;
    min_dist_thresh = 2;  %2

    [I1 J1] = find(obj_map == U(i));
    [rectx1,recty1,area1,perimeter1] = minboundrect(J1,I1);
       
    
    [I2 J2] = find(obj_map == U(j));
    [rectx2,recty2,area2,perimeter2] = minboundrect(J2,I2);
    
    if ( (size(rectx1,1) ~= 5 ) | (size(rectx2,1) ~= 5) )
        return;
    end
    
    [BW,x1,y1] = roipoly(obj_map,rectx1,recty1);
    [BW,x2,y2] = roipoly(obj_map,rectx2,recty2);
      
    [x0,y0] = intersections(x1,y1,x2,y2);
    
    %overlap_size = size(x0,1)
        
    if (size(x0,1) > 5 )
        overlap = 1;
    end
    
    %%%% if overlapping, then consider connected %%%%%%
    %[overlap] = check_overlapping(rectx1, rectx2, recty1, recty2);
    
   
    
    if (overlap == 1)
        connected = 1;
        return;
    end
    
   
    %%%% if not overlapped, then check for minimum distance 
    min_dist = 9999;
    
    for r1_pt = 1:4
        for r2_pt = 1:4
            
            x1 = uint8(rectx1(r1_pt));
            x2 = uint8(rectx2(r2_pt));
            
            y1 = uint8(recty1(r1_pt));
            y2 = uint8(recty2(r2_pt));
            
            x1 = double(x1);
            y1 = double(y1);
            x2 = double(x2);
            y2 = double(y2);
            
            dist = sqrt((x1-x2)^2 + (y1-y2)^2);
            
            if (dist <= min_dist )
               % disp(sprintf('x1 = %d, x2 = %d, y1 = %d, y2 = %d', x1,x2,y1,y2));
                min_dist = dist;
            end
        end
    end
    
    min_dist
    
    if ( min_dist < min_dist_thresh )
        connected = 1;
    end    
end


function [overlap] = check_overlapping(rectx1, rectx2, recty1, recty2)

    overlap = 0;
    
    overlap_threshold = 5;
    
    ax1 = min(rectx1); ax2 = max(rectx1);    
    bx1 = min(rectx2); bx2 = max(rectx2);
    
    ay1 = min(recty1); ay2 = max(recty1);    
    by1 = min(recty2); by2 = max(recty2);
   
%     if ( (bx1>ax2) | (by1>ay2) | (bx2<ax1) | (by2 <ay1) )
%         overlap = 0;
%     else
%         overlap = 1;
%     end  

    all_x = [ax1 ax2 bx1 bx2];
    all_y = [ay1 ay2 by1 by2];
    [X,Y] = meshgrid(all_x, all_y);

    % find which points are in which rectangle
    in_rect1 = inpolygon(X,Y, rectx1, recty1);
    in_rect2 = inpolygon(X,Y, rectx2, recty2);

    % rectangle overlap if there is at least one point that is in both
    % rectangles
    %overlap = any(in_rect1(:) & in_rect2(:))

    isect_pts = intersect( in_rect1(:), in_rect2(:));
    
    if (size(isect_pts,1) > overlap_threshold )
       overlap = 1; 
    end

end


function d = point_to_line(pt, v1, v2)
    a = v1 - v2;
    b = pt - v2;
    d = norm(cross(a,b)) / norm(a);
end
