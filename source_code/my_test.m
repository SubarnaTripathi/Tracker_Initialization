%my_test
% function my_test(input_motion_mask, moving_edge_mask)
% 
% blksize = 4;
% [filtered_output] = eps_filter(input_motion_mask, moving_edge_mask, blksize);
% %[filtered_output] = eps_filter(prob_mask1, moving_edge_mask, blksize);
% figure(100); imshow(input_motion_mask, []), title(' input motion mask')
% figure(200); imshow(moving_edge_mask, []), title(' input edge mask')
% figure(300); imshow(filtered_output, []), title(' smoothed motion mask')



function my_test(M)

mov = avifile('crick_init_grouping.avi', 'compression', 'none', 'fps',1);

for i=1:size(M,2)
    i
    F = M(i);
    if (i ~= 6 )
        mov = addframe(mov,F);
    end
end

mov = close(mov);