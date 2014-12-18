%function [t1,t2,t3,t4,t5,t6] = showframe_cluster(obj_map, k, M, fig_num, rgb,width,height, convex_hull)
function [obj_map_final, M] = showframe_cluster(obj_map,k, M, fig_num, rgb,width,height, convex_hull)

% pause
% movie(M)

obj_map_final = obj_map;

U = unique(obj_map);
tot_num_obj = size(U,1);
convexity_thresh = 0.3;

obj_map_new = obj_map(:);

U = unique(obj_map);
for i = 1:size(U,1)
      %%%% label 0
      obj_map_new(obj_map == U(i)) = i;
      obj_map_new = reshape(obj_map_new,size(obj_map));
end
obj_map = obj_map_new;

fig = rgb;
[H W C] = size(rgb);

figure(fig_num),imshow(fig, []);title(sprintf('show clusters')),hold on
figure(fig_num+1),imshow(fig, []);title(sprintf('show clusters')),hold on
figure(fig_num+2),imshow(fig, []);title(sprintf('show clusters')),hold on

for obj_num = 1:tot_num_obj    
    %obj_num

    [I J] = find(obj_map == obj_num);
    actual_area = double(size(I,1));
     
    if (actual_area < 30 )
        %obj_map_final(obj_map == U(obj_num))= 0;
        obj_map_final(obj_map == obj_num)= 0;
        % obj_num
        % disp('in continue');
        continue;
    end
           
    if ( 1 == convex_hull )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% new code : based on minimum bounded rectangle from convex hull
        theta = 0;
        [rectx,recty,area,perimeter] = minboundrect(J,I);
              
        if ( rectx(1)== rectx(4))
            width = abs(rectx(1)-rectx(2));
            height = abs(recty(1)-recty(4));
        else        
            theta = atan((recty(3)-recty(2))/(rectx(2)- rectx(3)));
            width = abs((rectx(3)-rectx(4))/sin(theta));
            height = abs((recty(2)-recty(3))/sin(theta));
        end
        
        area2 = height*width;
        
        
        %disp(sprintf('area =%d height*width =%d height = %d width = %d theta=%f', area, height*width, height, width, theta));
        
       
        p(1:2,1:2) = [rectx(1) rectx(2); recty(1) recty(2)];
        p(1:2,3:4) = [rectx(2) rectx(3); recty(2) recty(3)];
        p(1:2,5:6) = [rectx(3) rectx(4); recty(3) recty(4)];
        p(1:2,7:8) = [rectx(4) rectx(5); recty(4) recty(5)];
        cx = sum(rectx(:))/5;
        cy = sum(recty(:))/5;
        
       
%         t = [cx, cy, width, height, theta];
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% old code : orientation based on PCA    
        c1= min(J(:));
        c2= max(J(:));
        c_x = (c1+c2)/2;
        width = double(c2 - c1);

        r1= min(I(:));
        r2= max(I(:));
        c_y = (r1+r2)/2;
        height = double(r2 - r1);    

        theta = 0;
        theta = orient_cluster(I, J, obj_num, c_x, c_y, height);    

        p = [c1 c2 c2 c2 c2 c1 c1 c1,
             r1 r1 r1 r2 r2 r2 r2 r1,
             1  1  1  1  1  1  1  1];      

        if(abs(theta) >= 0.02)
            rot_mat = [cos(theta) -sin(theta) 0,
                   sin(theta) cos(theta)  0,
                   0          0           1];

            p(1,:) = p(1,:) - c_x;
            p(2,:) = p(2,:) - c_y;

            new_p = rot_mat*p;
            p= new_p;

            p(1,:) = p(1,:) + c_x;
            p(2,:) = p(2,:) + c_y;
        end
        area = width*height;
        perimeter = 2*(width*height);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %%% remove false positive
    area_fraction = actual_area/area;
    if ( area_fraction < convexity_thresh )
        %obj_map_final(obj_map == U(obj_num))= 0;
        obj_map_final(obj_map == obj_num)= 0;
        disp('too much concave so discarding');
        continue;
    end
        
%     PA = perimeter/area; 
%     if ( PA > 2)
%         obj_map_final(obj_map == U(obj_num))= 0;
%         disp('too long and thin object so discarding');
%         continue;
%     end
   
    if ( (width < 16 ) | (height < 16) )
        %obj_map_final(obj_map == U(obj_num))= 0;
        obj_map_final(obj_map == obj_num)= 0;
         %disp(sprintf('height =%d, width=%d', height, width));
         %obj_num
        disp('very small cluster so discarding');
        continue;
    end
    
    if ((width >= W - 5 ) | (height > H - 5) )
        %obj_map_final(obj_map == U(obj_num))= 0;
        obj_map_final(obj_map == obj_num)= 0;
        %disp(sprintf('height =%d, width=%d', height, width));
        disp('very large cluster so discarding');
        continue;
    end
        
    color = 'RED';    
    
   
    %%%% ploting object boundary 
    im_bw_img = zeros(size(obj_map));
    im_bw_img(obj_map == obj_num)= 1;
    im_bw_img = im2bw(im_bw_img);
    %figure(20), imshow(im_bw_img); 
    [B,L] = bwboundaries(im_bw_img, 'noholes');
    boundary = B{1};
    figure(fig_num); plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth',2);
    hold on
    
    figure(fig_num+2); plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth',2);
    hold on
    
    obj_map_final = reshape(obj_map_final,size(obj_map));
       
    %%%%% draw bounding box %%%%%%%%%%%%%%%%%
    figure(fig_num+1);
    x1=[p(1,1), p(1,2)]; x2=[p(2,1), p(2,2)];
    line(x1,x2, 'color', color),hold on
    x3=[p(1,3), p(1,4)]; x4=[p(2,3), p(2,4)];
    line(x3,x4, 'color', color),hold on
    x5=[p(1,5), p(1,6)]; x6=[p(2,5), p(2,6)];
    line(x5,x6, 'color', color),hold on
    x7=[p(1,7), p(1,8)]; x8=[p(2,7), p(2,8)];
    line(x7,x8, 'color', color),hold on 
    
    k = k+1;
    F = getframe();
    M(k) = F;
    %mov = addframe(mov,F);
    %mov = close(mov);
    
    %pause
    %movie(M);
    
    %movie2avi(mov,'grouping_video.avi', 'compression', 'NONE', 'fps',1);
    
    
    figure(fig_num+2);
    x1=[p(1,1), p(1,2)]; x2=[p(2,1), p(2,2)];
    line(x1,x2, 'color', color),hold on
    x3=[p(1,3), p(1,4)]; x4=[p(2,3), p(2,4)];
    line(x3,x4, 'color', color),hold on
    x5=[p(1,5), p(1,6)]; x6=[p(2,5), p(2,6)];
    line(x5,x6, 'color', color),hold on
    x7=[p(1,7), p(1,8)]; x8=[p(2,7), p(2,8)];
    line(x7,x8, 'color', color),hold on    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
end % end obj_num loop

%     new_fig = getframe(gcf);
%     figure(35),imshow(new_fig.cdata)
%     pause
    %OutSequenceL(fr_cnt+1) = new_fig;
    
    
%pause
%movie(M)

    

