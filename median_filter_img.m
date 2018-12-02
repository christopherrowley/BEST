function unpad = median_filter_img(img,windowsize)
%%%%%%%%%%%%%%%
%Median filter based on ordfilt3 by:
% (c) Toby Collins, Edinburgh University, UK, September 2008#
%     toby.collins@gmail.com
%     2008
%%%%%%%%%%%%%%%

img = init_seg;

if ~(mod(windowsize,2)==1)
    error('median_filter_img: windowsize must be odd.');
end

%values for use later
w = (windowsize-1)/2;
[x,y,z] = size(img); 
xend = x+2*w;
yend = y+2*w;
zend = z+2*w;

%pad image
pad = zeros((xend),(yend),(zend));
pad(1+w:x+w,1+w:y+w,1+w:z+w) = img; 


figure;
imshow3Dfull(pad)


%need to change the 3D window into a 1D string of value then median
med_img = zeros(size(pad));

tic
for i = 1+w:x
    for j = 1+w:y
        for k= 1+w:z
            vals = reshape(pad(i-w:i+w,j-w:j+w,k-w:k+w),[],1); %max val of a sliding window of 3
            med_img(i,j,k) = median(vals);  %set that value in the matrix 
        end
    end
end
toc % this took 404 seconds. 

% time wise, compare to the special result of a binary image, where it is
% equivalent to rounding the mean of the sampling window. 

med_img2 = zeros(size(pad));

tic
for i = 1+w:x
    for j = 1+w:y
        for k= 1+w:z
            vals = reshape(pad(i-w:i+w,j-w:j+w,k-w:k+w),[],1); %max val of a sliding window of 3
            med_img2(i,j,k) = round(mean(vals));  %set that value in the matrix 
        end
    end
end
toc %499 secons

dif_img = med_img - med_img2;

smoothed_seg = ordfilt3(img,'med',5);

figure;
imshow3Dfull(med_img)
figure;
imshow3Dfull(med_img2)

unpad = med_img(1+w:x+w,1+w:y+w,1+w:z+w);



















