function [ center_mask ] = center_remove(bin_img,shrink_vox)

% want to isolate the top center of the brain to work to remove sagital
% sinus. This provides a subject specific framework to find the upper
% portion of the center of the brain. 

%take off center brain slice to avoid brain split
% need bin_img
%  bin_img = bin_img_fill;
[x y z] = size(bin_img);

slice = squeeze(bin_img((x/2+10),:,:));

%rows are back to front brain, columns are bottom to top brain.
top_idx = zeros(1,size(slice,1));
for i = 1:size(top_idx,2)
    if max(slice(i,:)) <1 %if there is no tissue in that row
        top_idx(i) = 0;
    else %if there is, find the top idx for it.
        top_idx(i) = find(slice(i,:),1,'last');
    end 
end

% figure;
% plot(top_idx)

smooth_idx = smooth(top_idx,20);
% figure;
% plot(smooth_idx)

% Shrink voxels number
%shrink_vox = 18; 

fst = find(top_idx,1);
lst =find(top_idx,1,'last'); 
smooth_idx(1:fst+shrink_vox) = 0;
smooth_idx(lst-shrink_vox:end) = 0;
smooth_idx = round(smooth_idx);
%take it down 15
smooth_idx = smooth_idx - shrink_vox;
smooth_idx(smooth_idx<0) = 0;

slice_mask = zeros(size(slice));
%fill in with 1's

for i = 1:size(slice_mask,1)
    if smooth_idx(i) > 1 %if there is no tissue in that row
        slice_mask(i,1:smooth_idx(i)) = 1;
    end 
end

%take the mask and make it the size of the image. 
center_mask = ones(size(bin_img));
for i = (x/2-20):(x/2+20)
    center_mask(i,:,:) = slice_mask;
end

% %troubleshooting visualization
% figure;
% imshow3Dfull(center_mask)
% 
% check_result  = (1-center_mask) .* bin_img_fill;
% figure;
% imshow3Dfull(check_result)
% 
% figure;
% imagesc(slice_mask);


