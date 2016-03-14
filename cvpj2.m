im1 = imread('DSC_0281.JPG');
im1 = rgb2gray(im1);
im2 = imread('DSC_0282.JPG');
im2 = rgb2gray(im2);
%harris compute R function and use non-max suppression to get corner pts 
[R2, x2, y2] = harris('DSC_0282.JPG',3,100000,4,1);
[R1, x1, y1] = harris('DSC_0281.JPG',3,100000,4,1);
corn1 = [x1  y1]';
corn2 = [x2  y2]';

i=length(x1);
j=length(x2);
Cnox = zeros(i,j);
%count=0;
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
t=0.95;
%%a and b are corresponding points from 2 
[a, b]=find(Cnox>t);
countx=length(a);
for p = 1:(countx-3)
%i=1;
%find a set of 4 points that match! change countx to increase number of pts
%while countx > 
%    t = t+0.001;
 %   [a, b]=find(Cnox>t);
%    i=i+1;
 %   countx=length(a);
%end
%write the A matrix to compute homography h

A = [x1(a(p))   y1(a(p))  1     0           0         0    -x1(a(p))*x2(b(p))    -y1(a(p))*x2(b(p));
    0           0         0     x1(a(p))    y1(a(p))  1    -x1(a(p))*y2(b(p))    -y1(a(p))*y2(b(p));
    x1(a(p+1))   y1(a(p+1))   1     0           0         0    -x1(a(p+1))*x2(b(p+1))    -y1(a(p+1))*x2(b(p+1));
    0            0        0     x1(a(p+1))    y1(a(p+1))  1    -x1(a(p+1))*y2(b(p+1))    -y1(a(p+1))*y2(b(p+1));
    x1(a(p+2))   y1(a(p+2))   1     0           0         0    -x1(a(p+2))*x2(b(p+2))    -y1(a(p+2))*x2(b(p+2));
    0           0         0     x1(a(p+2))    y1(a(p+2))  1    -x1(a(p+2))*y2(b(p+2))    -y1(a(p+2))*y2(b(p+2));
    x1(a(p+3))   y1(a(p+3))   1     0           0         0    -x1(a(p+3))*x2(b(p+3))    -y1(a(p+3))*x2(b(p+3));
    0           0         0     x1(a(p+3))    y1(a(p+3))  1    -x1(a(p+3))*y2(b(p+3))    -y1(a(p+3))*y2(b(p+3))];
    
X = [x2(b(p)) ; y2(b(p)); x2(b(p+1)); y2(b(p+1)); x2(b(p+2)); y2(b(p+2)); x2(b(p+3)); y2(b(p+3))];
h = A\X;
%AllH(:,p) = h;
H = [h(1:3) h(4:6) [h(7:8); 1]]';
%project using homography and compute the dist
s = size(corn1,2);
temp= [corn1;ones(1,s)];
pts3 = H*temp;
pts3 = pts3(1:2,:)./repmat(pts3(3,:),2,1);

d(:, p) = sum((corn2(:,1:104)-pts3).^2,1);
%inlierid = find(d < 4.0);
end
[r c] = find(d < 4.0);  %threshold to find min distance
j = 1;
%get the unique points!!!
for i = 2:length(r)
   if (r(i) ~= r(i-1)) 
       inlinerid(j,:) = r(i);
       j = j+1;
   end
end
pt = Cnox(inlinerid,:);

%%%get the matching pt again!!!
[a,b] = find(pt >t);
count=length(a);
chosenpt1 = corn1(:,a);
chosenpt2 = corn2(:,b);
%%find homography again
a = chosenpt1;
b = chosenpt2;
x1 = a(1,:); 
x2 = b(1,:);
y1 = a(2,:);
y2 = b(2,:);
for p = 1:(count-3)
    
%i=1;
%find a set of 4 points that match! change countx to increase number of pts
%while countx > 
%    t = t+0.001;
 %   [a, b]=find(Cnox>t);
%    i=i+1;
 %   countx=length(a);
%end
%write the A matrix to compute homography h

A =[x1(p)       y1(p)     1     0           0         0    -x1(p)*x2(p)        -y1(p)*x2(p);
    0           0         0     x1(p)       y1(p)     1    -x1(p)*y2(p)        -y1(p)*y2(p);
    x1(p+1)     y1(p+1)   1     0           0         0    -x1(p+1)*x2(p+1)    -y1(p+1)*x2(p+1);
    0           0         0     x1(p+1)     y1(p+1)   1    -x1(p+1)*y2(p+1)    -y1(p+1)*y2(p+1);
    x1(p+2)     y1(p+2)  1     0           0         0    -x1(p+2)*x2(p+2)    -y1(p+2)*x2(p+2);
    0           0         0     x1(p+2)     y1(p+2)   1    -x1(p+2)*y2(p+2)    -y1(p+2)*y2(p+2);
    x1(p+3)     y1(p+3)   1     0           0         0    -x1(p+3)*x2(p+3)    -y1(p+3)*x2(p+3);
    0           0         0     x1(p+3)    y1(p+3)   1    -x1(p+3)*y2(p+3)    -y1(p+3)*y2(p+3)];
    
X = [x2(p) ; y2(p); x2(p+1); y2(p+1); x2(p+2); y2(p+2); x2(p+3); y2(p+3)];
h = A\X;
AllH(:,p) = h;
H = [h(1:3) h(4:6) [h(7:8); 1]]';
%project using homography and compute the dist
s = size(chosenpt1,2);
temp= [chosenpt1;ones(1,s)];
pts3 = H*temp;
pts3 = pts3(1:2,:)./repmat(pts3(3,:),2,1);
dis(:, p) = sum((chosenpt2(:,1:s)-pts3).^2,1);
%inlierid = find(d < 4.0);
end



%patch1 = im1((row1(1)-1):(row1(1)+1),(col1(1)-1):(col1(1)+1));
%patch2 = im2((row2(6)-1):(row2(6)+1),(col2(6)-1):(col2(6)+1));
%patch3 = im2((row2(6)-2):(row2(6)+2),(col2(6)-2):(col2(6)+2));

%c1 = normxcorr2(patch1, patch2);

% [features1, valid_points1] = extractFeatures(im1, double(corn1));
% [features2, valid_points2] = extractFeatures(im2, double(corn2));
% indexPairs = matchFeatures(features1, features2);
%c = normxcorr2(R1, R2);
%indexPairs = matchFeatures(R1, R2);
%matchedPoints1 = corn1(c(:, 1), :);
%matchedPoints2 = corn2(c(:, 2), :);
% matchedPoints1 = valid_points1(indexPairs(:, 1), :);
% matchedPoints2 = valid_points2(indexPairs(:, 2), :);
%figure; showMatchedFeatures(im1, im2, matchedPoints1, matchedPoints2);
