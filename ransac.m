 function [f inlierIdx inlrNum] = ransac(x,y,ransacCoef)
%[f inlierIdx] = ransac1( x,y,ransacCoef,funcFindF,funcDist )
%	Use RANdom SAmple Consensus to find a fit from X to Y.
%	X is M*n matrix including n points with dim M, Y is N*n;
%	The fit, f, and the indices of inliers, are returned.
%
%	RANSACCOEF is a struct with following fields:
%	minPtNum,iterNum,thDist,thInlrRatio
%	MINPTNUM is the minimum number of points with whom can we 
%	find a fit, For homography, it's 4.
%	ITERNUM is the number of iteration, 
%	distance threshold and ROUND(THINLRRATIO*n) is the inlier number threshold.


minPtNum = ransacCoef.minPtNum;
iterNum = ransacCoef.iterNum;
thDist = ransacCoef.thDist;
ptNum = size(x,2);

inlrNum = zeros(1,iterNum);
fLib = cell(1,iterNum);

for p = 1:iterNum
	% 1. fit using  random points
	sampleIdx = randIndex(ptNum,minPtNum); %get 4 random pts
	f1 = solveHomo(x(:,sampleIdx),y(:,sampleIdx));
	
	% 2. count the inliers, if more than thInlr, refit; else iterate
	dist = calcDist(f1,x,y);
	inlier1 = find(dist < thDist);
	inlrNum(p) = length(inlier1);
	%if length(inlier1) < 16, continue; end
	fLib{p} = solveHomo(x(:,inlier1),y(:,inlier1));
end

% 3. choose the coef with the most inliers
[~,idx] = max(inlrNum);
f = fLib{idx}; %get the H at that location
dist = calcDist(f,x,y);
inlierIdx = find(dist < thDist);

end