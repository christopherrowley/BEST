function norm_img3 = PV_divison(img)
%%%%%%%%%%%%%%%
%Image normalization
%%%%%%%%%%%%%%%

%goal is to remove myelin and b1 driven changes in GM intensity and the
%voxels in which they are PV'd with CSF

%%%%%%%%%%%%%%%
%try in 3D

check_image = img; 
img_size = size(check_image); 
img_max = zeros(img_size);

ref_img = zeros(img_size(1)+2,img_size(2)+2,img_size(3)+2); %pad the reference image for the moving sampling
ref_img(1:img_size(1),1:img_size(2),1:img_size(3)) = check_image; 

for i = 1:img_size(1)
    for j = 1:img_size(2)
        for k= 1:img_size(3)
            val = max(max(max(ref_img(i:i+2,j:j+2,k:k+2)))); %max val of a sliding window of 3
             img_max(i,j,k) = val;  %set that value in the matrix 
        end
    end
end


img_max_cent = zeros(img_size);
img_max_cent(2:img_size(1),2:img_size(2),2:img_size(3)) = img_max(1:(img_size(1)-1),1:(img_size(2)-1),1:(img_size(3)-1));

Iblur3 = imgaussfilt3(img_max_cent,10);

%adjust Iblur to remove minimal values
val = mean(max(max(Iblur3)));
Iblur3(Iblur3 < val) = val;

norm_img3 = check_image./Iblur3;



















