
function [R, row, col] = harris(im, sigma, thresh, radius, disp)
    im = imread(im);
    im = rgb2gray(im);
    error(nargchk(2,5,nargin));
    
    dx = [-1 0 1; -1 0 1; -1 0 1]; % Derivative masks
    dy = dx';
    
    Ix = imfilter(double(im), double(dx));    % Image derivatives
    Iy = imfilter(double(im), double(dy));    

    % Generate Gaussian filter of size 6*sigma (+/- 3sigma) and of
    % minimum size 1x1.
    g = fspecial('gaussian',max(1,fix(6*sigma)), sigma);
    
    Ix2 = imfilter(double(Ix.^2), double(g)); % Smoothed squared image derivatives
    Iy2 = imfilter(double(Iy.^2), double(g));
    Ixy = imfilter(double(Ix.*Iy), double(g));
    
   % R = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps); % Harris corner measure

    % Alternate Harris corner measure used by some.  Suggested that
    % k=0.04 - I find this a bit arbitrary and unsatisfactory.
   R = (Ix2.*Iy2 - Ixy.^2) - 0.04*(Ix2 + Iy2).^2; 

    if nargin > 2   % We should perform nonmaximal suppression and threshold
	
	% Extract local maxima by performing a grey scale morphological
	% dilation and then finding points in the corner strength image that
	% match the dilated image and are also greater than the threshold.
	sze = 2*radius+1;                   % Size of mask.
	mx = ordfilt2(R,sze^2,ones(sze)); % Grey-scale dilate.
	Rm = (R==mx)&(R>thresh);       % Find maxima.
	
	[row,col] = find(Rm);                  % Find row,col coords.
	
	if nargin==5 & disp      % overlay corners on original image
	    figure, imagesc(im), axis image, colormap(gray), hold on
	    plot(col,row,'rs'), title('corners detected');
	end
    
    else  % leave R as a corner strength image and make r and c empty.
% 	r = []; c = [];
    end
    