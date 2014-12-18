function fig = show_cluster(Y, height, width, k )

box_color1 = [255, 0, 0]; %RED
box_color2 = [0, 255, 0]; %BLUE
box_color3 = [0, 0, 255]; %GREEN
box_color4 = [255, 255, 0];%YELLOW
box_color5 = [0, 255, 255];%CYAN
box_color6 = [255, 0, 255];%MAGENTA

pts_size = size(Y);
tot_pts = pts_size(1);
fig = zeros(height, width, 3);

for num_pts = 1:tot_pts
%    curr_box_color = sprintf('color%d',Y(num_pts,3));      
%     switch(curr_box_color)
%         case 'color1'
%            fig( Y(num_pts,1), Y(num_pts,2),:) = box_color1;
%         case 'color2'
%            fig( Y(num_pts,1), Y(num_pts,2),:) = box_color1; 
%         case 'color3'
%            fig( Y(num_pts,1), Y(num_pts,2),:) = box_color3; 
%         case 'color4'
%            fig( Y(num_pts,1), Y(num_pts,2),:) = box_color4; 
%         case 'color5'
%            fig( Y(num_pts,1), Y(num_pts,2),:) = box_color5; 
%         case 'color6'
%            fig( Y(num_pts,1), Y(num_pts,2),:) = box_color6;  
%         otherwise
%            disp('Unknown color.')
%     end
     
    color_index = Y(num_pts,3);
    if ( color_index == 1 )
        fig( Y(num_pts,1), Y(num_pts,2),:) = box_color1;
    elseif ( color_index == 2 )    
        fig( Y(num_pts,1), Y(num_pts,2),:) = box_color2;
    elseif ( color_index == 3 )
        fig( Y(num_pts,1), Y(num_pts,2),:) = box_color3;
    elseif ( color_index == 4 )     
        fig( Y(num_pts,1), Y(num_pts,2),:) = box_color4;
    elseif  ( color_index == 5 )
        fig( Y(num_pts,1), Y(num_pts,2),:) = box_color5;
    elseif  ( color_index == 6 )
        fig( Y(num_pts,1), Y(num_pts,2),:) = box_color6;
    end
end

%figure(10),imshow(fig),title('showing only clusters')