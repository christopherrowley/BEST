
img = norm_img .*init_seg;
[h, edges] = histcounts(img(img>0.005));
edges = edges(2:end); %x -axis, remove extra value

figure; %trouble shooting, visualize it
plot(edges',h');


temp = img_sharp;
temp(temp>3) = 3;
[h, edges] = histcounts(temp(temp>0.05));
edges = edges(2:end); %x -axis, remove extra value

figure; %trouble shooting, visualize it
plot(edges',h');