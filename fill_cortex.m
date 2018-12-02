function [ filled_mask ] = fill_cortex(input_mask)


%Fill cortex - meant to fill between the edges

%input_mask = smoothed_seg;
[x,y,z] = size(input_mask);

filled_mask = zeros(x,y,z);
%take axial slices which is the z dimension
for k = 1:z
    %in each back to front column = set x value
    for i = 1:x 
        %find first and last idx
        fst_y = find(input_mask(i,:,k),1);
        lst_y = find(input_mask(i,:,k),1,'last'); 
        %set = 1 in between idx
        filled_mask(i,fst_y:lst_y,k) = 1;
    end
end 
