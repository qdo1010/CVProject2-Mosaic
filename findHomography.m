function [H corrPtIdx inln] = findHomography(pts1,pts2)
% [H corrPtIdx] = findHomography(pts1,pts2)
%	Find the homography between two planes using a set of corresponding
%	points. PTS1 = [x1,x2,...;y1,y2,...]. RANSAC method is used.
%	corrPtIdx is the indices of inliers.

coef.minPtNum = 4;
coef.iterNum = 500;
coef.thDist = 4;
%use Ransac
[H corrPtIdx inln] = ransac(pts1,pts2,coef);

end

