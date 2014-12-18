%function [obj_map_final, edge_mask, color_mask, problist, im, M] = tracker_init(edge_mask, problist, problist_pp, im, color_mask)
%function [obj_map_final, edge_mask, color_mask, problist, im, problist_pp, M] = tracker_init(edge_mask)
function [obj_map_final, edge_mask, color_mask, problist, im, problist_pp, M] = tracker_init

InFileName = 'toRoad.avi'; %'crick.avi'; %'foreman_cif.avi'; %'SECOND.avi'; %'coastguard_qcif.avi'; %'films_CIF_part2_fr16.avi'; %'skydiving_CIF.avi'; %'renata_qcif.avi'; %'md.avi'; %'akio_cif.avi'; %'films_CIF_part1_fr15.avi'; %'hand_test.avi'; %'films_CIF_part2_fr16.avi';   

start_frame = 20; %20
end_frame = start_frame+5;
blksize = 4;

[blklist, problist, edge_mask] = motion_segmentation_new (InFileName, blksize, start_frame, end_frame, edge_mask); %color_mask
figure(2); imshow(edge_mask, []); title('edge mask')
%pause

%figure(6); imshow(blklist, []), title('blklist');
figure(7); imshow(problist, []),title('problist');
%pause

[im, color_mask] = colored_region_by_graph_cut_4plane(InFileName,start_frame,blksize,problist); %end_frame
figure(3); imshow(color_mask, []); title('segmentated frame')

%pause

[height width] = size(color_mask);
Y = zeros(height*width, 3);

U = unique(color_mask);
for i = 1:size(U,1)
      color_mask_new(color_mask == U(i)) = i;
      color_mask_new = reshape(color_mask_new,size(color_mask));
end


[obj_map, num_levels, problist_pp] = assign_obj_label(color_mask_new, problist, edge_mask, blksize, im);
figure(10), imshow(obj_map, []), title('obj map')

figure(40), imshow(problist_pp, []), title('problist pp')

pause

%%% obj_map after quantization
max_val = max(max(obj_map));
obj_map = (obj_map/max_val)*255;
obj_map = double(uint8(obj_map));
obj_map_merged = obj_map;

%pause

%%%% merge possible clusters
%%% obj_map after merging and post-processing
imwrite(obj_map, 'obj_map1.ppm', 'pgm', 'encoding', 'ASCII', 'MaxValue',255);
[obj_map_merged, num_obj,k, M] = merge_clusters_new(im, obj_map, obj_map, color_mask, problist, problist_pp, blksize); %problist

figure(11), imshow(obj_map_merged, []), title('obj map after merging')
imwrite(obj_map_merged, 'obj_map.pgm', 'pgm', 'encoding', 'ASCII', 'MaxValue',255);


%%% obj_map after discarding spurious regions
fig_num  = 20;
[obj_map_final, M] = showframe_cluster(obj_map_merged, k, M, fig_num, im, width, height,1); % last parameter is for BB through convex-hull or positional PCA
figure(11), imshow(obj_map_final, []), title('final obj map during tracker initialization')
end







