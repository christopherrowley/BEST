function [ blurred_seg ] = split_and_smooth(seg_to_blur, kernel)
%input is a brain label with 1's on one side, and 2's on the other
% kernel = 3D gaussian kernel you want to use

% seg_to_blur = half_brain;
% kernel = [1 20 1];

right_blur = seg_to_blur == 1;
left_blur = seg_to_blur == 2;

right_blur = imgaussfilt3(double(right_blur),kernel);
right_blur = round(right_blur); %rounds it 
left_blur = imgaussfilt3(double(left_blur),kernel);

left_blur2 = 2*round(left_blur);
%left_blur = 2*left_blur.*(left_blur>=0.5);
blurred_seg = left_blur2+right_blur;

%max(max(max(blurred_seg)))

% figure;
% imshow3Dfull(blurred_seg)

end
