function [filtered_output] = eps_Filter(input_motion_mask, edge_mask, blksize)

filtered_output = input_motion_mask;
edge_threshold = 0.4;
[h,w] = size(input_motion_mask);

for i = 1:h
    for j = 1:w 
    
        indices = j+(i-1)*w;
        adjBlocks = [ indices;  %center  
                indices-w-1;    % NW corner
                indices-w;      % top
                indices-w+1;    % NE corner
                indices-1;      % left
                indices+1;      % right
                indices+w-1;    % SW corner
                indices+w;      % bottom
                indices+w+1];   % SE corner
                   
        % mark the valid indices    
        if(i == 1) adjBlocks(2)= 0; adjBlocks(3)= 0; adjBlocks(4)= 0; end   % first row means no top, NE, NW
        if(j == 1) adjBlocks(2) = 0;adjBlocks(5) = 0;adjBlocks(7) = 0; end % first column means no left, SW, NW 
        if(i == h) adjBlocks(7)= 0; adjBlocks(8)= 0; adjBlocks(9)= 0; end   % last row means no bottom, SW, SE
        if(j == w) adjBlocks(4) = 0;adjBlocks(6) = 0;adjBlocks(9) = 0; end % last column means no NE,right, SE 
        
        % pick out only the valid indices 
        validOnes = find((adjBlocks<(h*w)) & (adjBlocks>0));        
        adjGroups = adjBlocks(validOnes);
                            
        %%% if non-edge
        if ( edge_mask(i,j) < edge_threshold)   
               filtered_output(i,j) = median(input_motion_mask(adjGroups));     
               %filtered_output(i,j) = mean(input_motion_mask(adjGroups));     
        end    
    end
end