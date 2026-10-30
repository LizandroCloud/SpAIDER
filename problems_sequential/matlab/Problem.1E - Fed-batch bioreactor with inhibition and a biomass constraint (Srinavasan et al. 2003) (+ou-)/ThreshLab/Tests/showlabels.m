
labelmap = [(0:255)' (0:255)' (0:255)']/255;
figure(1)
x1 = zeros(size(input));
x1(:,:) = 256-label(j,1,:,:)*256;
image(x1)
colormap(labelmap)

figure(2)
x2 = zeros(size(input));
x2(:,:) = 256-label(j,2,:,:)*256;
image(x2)
colormap(labelmap)

figure(3)
x3 = zeros(size(input));
x3(:,:) = 256-label(j,3,:,:)*256;
image(x3)
colormap(labelmap)

j = j-1;

