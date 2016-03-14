function [result] = testpj2(images)

im1 = imread(images{1});
%im1 = rgb2gray(im1);
im2 = imread(images{2});
%im2 = rgb2gray(im2);
im3 = imread(images{3});
%im3 = rgb2gray(im3);

%im1 = images{3};
%im2 = images{3};
%im3 = images{3};
%harris compute R function and use non-max suppression to get corner pts 
%[R1, x1, y1] = harris('DSC_0281.JPG',3,5000,8,1);
%[R2, x2, y2] = harris('DSC_0282.JPG',3,5000,8,1);
%[R3, x3, y3] = harris('DSC_0283.JPG',3,5000,8,1);

[x1, y1, ~] = harris1(im1);
[x2, y2, ~] = harris1(im2);
[x3, y3, ~] = harris1(im3);
%im1 = rgb2gray(im1);
%im2 = rgb2gray(im2);
%im3 = rgb2gray(im3);

corn1 = [x1  y1]';
corn2 = [x2  y2]';
corn3 = [x3  y3]';

RADIUS = 10;
ENLARGE_FACTOR = 1.5;
circles = [x1 y1 RADIUS*ones(length(x1), 1)];
sift1 = findSift(im1, circles, ENLARGE_FACTOR);
siftLocations1 = [x1 y1];
 
circles = [x2 y2 RADIUS*ones(length(x2), 1)];
sift2 = findSift(im2, circles, ENLARGE_FACTOR);
siftLocations2 = [x2 y2];

circles = [x3 y3 RADIUS*ones(length(x3), 1)];
sift3 = findSift(im3, circles, ENLARGE_FACTOR);
siftLocations3 = [x3 y3];

%%Match sift1 w sift2
MATCH_THRESHOLD = 0.5;
 matches1 = match(sift1, sift2, MATCH_THRESHOLD);
  corres1 = siftLocations1(matches1(:, 1), :);
  corres2 = siftLocations2(matches1(:, 2), :);
  
 RANSAC_N = 200;
 RANSAC_EPSILON = 3;
 RANSAC_FINISH_THRES = 0.6;
    
  %find homography
%[H1, corrPtIdx] = ransac1(corres1, corres2, RANSAC_N, RANSAC_EPSILON, RANSAC_FINISH_THRES);
  
[H1 corrPtIdx inln] = findHomography(corres1',corres2');
 mappoint(im1, corres1', corres2', corrPtIdx, 'Image 1 - 2');


matches2 = match(sift3, sift2, MATCH_THRESHOLD);
  corres3 = siftLocations3(matches2(:, 1), :);
  corres4 = siftLocations2(matches2(:, 2), :);
  %find homography
%[H2, corrPtIdx2] = ransac1(corres4, corres4, RANSAC_N, RANSAC_EPSILON, RANSAC_FINISH_THRES);
[H2 corrPtIdx2 inlnj] = findHomography(corres3',corres4');
mappoint(im3, corres3', corres4', corrPtIdx2, 'Image 3 - 2');


images = {im1 im2 im3};
bbox = getBoundingBox([0 0 0 0], images, {H1, eye(3), H2});
l = abs(bbox(1));    t = abs(bbox(3));
bbox = bbox + [-l/4 l -t/2 -t/2];
bb_xmin = bbox(1);
bb_xmax = bbox(2);
bb_ymin = bbox(3);
bb_ymax = bbox(4);




%im1 = rgb2gray(im1);

%im2 = rgb2gray(im2);

%im3 = rgb2gray(im3);
%{
%map 1 to 2
h = inv(H1);
[xi, yi] = meshgrid(bb_xmin:bb_xmax, bb_ymin:bb_ymax);
xx = (h(1,1)*xi+h(1,2)*yi+h(1,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
yy = (h(2,1)*xi+h(2,2)*yi+h(2,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
im = double(im1);
foo = uint8(interp2(im,xx,yy));
%figure; imshow(foo)

h = eye(3);
[xk, yk] = meshgrid(bb_xmin:bb_xmax, bb_ymin:bb_ymax);
xx = (h(1,1)*xk+h(1,2)*yk+h(1,3))./(h(3,1)*xk+h(3,2)*yk+h(3,3));
yy = (h(2,1)*xk+h(2,2)*yk+h(2,3))./(h(3,1)*xk+h(3,2)*yk+h(3,3));
imk = double(im2);
foo1 = uint8(interp2(imk,xx,yy));

%map 3 to 2
hj = inv(H2);
[xj, yj] = meshgrid(bb_xmin:bb_xmax, bb_ymin:bb_ymax);
xxj = (hj(1,1)*xj+hj(1,2)*yj+hj(1,3))./(hj(3,1)*xj+hj(3,2)*yj+hj(3,3));
yyj = (hj(2,1)*xj+hj(2,2)*yj+hj(2,3))./(hj(3,1)*xj+hj(3,2)*yj+hj(3,3));
imj = double(im3);
foo2 = uint8(interp2(imj,xxj,yyj));
figure;
imshow([foo foo1 foo2]);
%}
%im1 = imread(images{1});
%im1 = rgb2gray(im1);
%im2 = imread(images{2});
%im2 = rgb2gray(im2);
%im3 = imread(images{3});

Im2w = vgg_warp_H(double(im2), eye(3), 'linear', bbox);
Im1w = vgg_warp_H(double(im1), H1, 'linear', bbox);
Im3w = vgg_warp_H(double(im3), H2, 'linear', bbox);
   figure();
imagesc(double(max(max(Im1w,Im2w), Im3w))/255)
result = (double(max(max(Im1w,Im2w), Im3w))/255);
%{
%pts1 = corn1;
%pts2 = corn2(:,1:size(pts1,2));
%pts3 = corn3(:,1:size(pts1,2));

i=length(x1);
j=length(x2);
Cnox = zeros(i,j);
for m = 1:i
    %3x3 neighborhood
    patch1 = im1((x1(m)-1):(x1(m)+1),(y1(m)-1):(y1(m)+1));
    for n = 1:j
        patch2 = im2((x2(n)-1):(x2(n)+1),(y2(n)-1):(y2(n)+1));
        c = normxcorr2(patch1, patch2);
        Cnox(m,n) = c(3,3); 
        %count=count+1;
    end
end
%x1,y1 are x,y of pt in first img, and x2, y2 are x,y of pt in second img
t=0.96;
%%a and b are corresponding points from 2 
[a, b]=find(Cnox>t);

countx=length(a);
pts1 = corn1(:, a);
pts2 = corn2(:, b);

%%map im1 to im2
[H corrPtIdx inln] = findHomography(pts1,pts2);

H1 = solveHomo(pts1(:,corrPtIdx), pts2(:,corrPtIdx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%map im3 to im2

i=length(x3);
j=length(x2);
Cnoxj = zeros(i,j);
for m = 1:i
    %3x3 neighborhood
    patch3 = im3((x3(m)-1):(x3(m)+1),(y3(m)-1):(y3(m)+1));
    for n = 1:j
        patch2 = im2((x2(n)-1):(x2(n)+1),(y2(n)-1):(y2(n)+1));
        c = normxcorr2(patch3, patch2);
        Cnoxj(m,n) = c(3,3); 
        %count=count+1;
    end
end

t=0.97;
%%c and d are corresponding points from 2 
[c, d]=find(Cnoxj>t);
%countx=length(c);
pts3 = corn3(:, c);
pts4 = corn2(:, d);

[Hj corrPtIdxj inlnj] = findHomography(pts3,pts4);

H2 = solveHomo(pts3(:,corrPtIdxj), pts4(:,corrPtIdxj));

%%map 1 to 2

%figure; imshow([foo foo2]) %%%1 to 2

%figure; imshow(foo2); %%3 to 2
  
 %Im2w = vgg_warp_H(double(im2), eye(3), 'linear', bbox);
 %Im1w = vgg_warp_H(double(im1), H, 'linear', bbox);
 %Im3w = vgg_warp_H(double(im3), Hj, 'linear', bbox);
 images = {im1 im2 im3};
bbox = getBoundingBox([0 0 0 0], images, {H, eye(3), Hj});
l = abs(bbox(1));    t = abs(bbox(3));
bbox = bbox + [-l/4 l -t/2 -t/2];
bb_xmin = bbox(1);
bb_xmax = bbox(2);
bb_ymin = bbox(3);
bb_ymax = bbox(4);


%map 1 to 2
h = inv(H1);
[xi, yi] = meshgrid(1:512,1:340);
xx = (h(1,1)*xi+h(1,2)*yi+h(1,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
yy = (h(2,1)*xi+h(2,2)*yi+h(2,3))./(h(3,1)*xi+h(3,2)*yi+h(3,3));
im = double(im1);
foo = uint8(interp2(im,xx,yy));
%figure; imshow(foo)

h = eye(3);
[xk, yk] = meshgrid(1:512,1:340);
xx = (h(1,1)*xk+h(1,2)*yk+h(1,3))./(h(3,1)*xk+h(3,2)*yk+h(3,3));
yy = (h(2,1)*xk+h(2,2)*yk+h(2,3))./(h(3,1)*xk+h(3,2)*yk+h(3,3));
imk = double(im2);
foo1 = uint8(interp2(imk,xx,yy));

%map 3 to 2
hj = inv(H2);
[xj, yj] = meshgrid(1:512,1:340);
xxj = (hj(1,1)*xj+hj(1,2)*yj+hj(1,3))./(hj(3,1)*xj+hj(3,2)*yj+hj(3,3));
yyj = (hj(2,1)*xj+hj(2,2)*yj+hj(2,3))./(hj(3,1)*xj+hj(3,2)*yj+hj(3,3));
imj = double(im3);
foo2 = uint8(interp2(imj,xxj,yyj));
figure;
imshow([foo foo1 foo2]);
%figure;
%imagesc((max(max(Im1w,Im2w), Im3w)));

   % Im2w = vgg_warp_H(double(imread('DSC_0281.JPG')), eye(3), 'linear', bbox);
   % Im1w = vgg_warp_H(double(imread('DSC_0282.JPG')), H, 'linear', bbox);
   % Im3w = vgg_warp_H(double(imread('DSC_0283.JPG')), Hj, 'linear', bbox);
    
 %   figure();
  %  imagesc(double(max(max(Im1w,Im2w), Im3w))/255);
%}
end