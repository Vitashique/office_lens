function [g] = my_lens ( f, x )

%problems: 12, 13, 14, 16
% im_g = imread('data/Office_Lens_Input_1.jpg');
% x = [1000, 1000];

%getting size of the input image
% Ro = size(f,1);
% Co = size(f,2);

% convert to grayscale
im_g = rgb2gray (f);

% histogram equalization. 
%change equalization if possible to get better grade
hist   = imhist(im_g(:,:)); 
maxIntensity = 255;  
cdf_v(1)= hist(1);
for i=2:maxIntensity+1
    cdf_v(i) = hist(i) + cdf_v(i-1);
end

cdf_v = cdf_v/double(numel(im_g))*255.0;

% converts image into uint8
im_e = uint8( cdf_v ( im_g+1));

% do not crop.  convert to color
%g = cat ( 3, im_e, im_e, im_e );

% Threshold gets a binary (black & white) image (only 0's and 1's)
threshold = 100;
binaryIm = im_e > threshold; %white = ">", black = "<"

% imfill "fills image regions and holes" (from matlab documentation) with white
% Basically gets rid of all the text on the paper
% binaryIm = imfill(binaryIm, 'holes');

% Displays the binary image.
% imshow(binaryIm);

% s = regionprops(binaryIm, 'BoundingBox');
% imCrop = imcrop(binaryIm, s.BoundingBox);

%detects the 4 corners in the binary image and returns them in a matrix
% matrix=corner(binaryIm,4);
% corner1x=matrix(1);
% corner1y=matrix(5);
% disp(corner1x);
% disp(matrix);
% I2 = imcrop(binaryIm,[2290 1170 40 40]);
% %imshow(binaryIm), figure, imshow(I2)

% binaryImage = poly2mask(matrix);
% measurements = regionprops(binaryImage, 'BoundingBox');
% cropped = imcrop(f, measurements.BoundingBox);
label = bwlabel(binaryIm, 8); 
%imshow(label, []); 
measure = regionprops(label, im_e, 'all');
%size = size(measure, 1);

% bwboundaries() returns an array which contains the row/column coordinates of the whole object
boundaries = bwboundaries(binaryIm,4); %bwboundaries traces the boundaries/edges of objects
sizeBoundaries = size(boundaries,1); %size() "returns the sizes of each boundary image"
for k = 1 : sizeBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end

pixels = measure.PixelIdxList;

intensities = [measure.MeanIntensity];
areas = [measure.Area];
%intensities have to be between 150-220 and areas over 1000
allowIntensities = (intensities > 150) & (intensities < 220);
allowArea = areas > 1000; % Take the small objects.
keepThese = find(allowIntensities & allowArea);
% Extract whatever matched our ranges and eliminate the noise
keepObject = ismember(label, keepThese);
labeled = bwlabel(keepObject, 8);
%imshow(labeled, []);

%find the paper's bounding box to make a rectangle and then crop it
boundingBox = measure.BoundingBox; % Get list of pixels
cropped = imcrop(im_e, boundingBox);

g = cropped;

%need code for cropping, contrast/clarity, imrotate

% resize
g = imresize ( g, x );
% g = imresize ( g, sizeBoundaries );

% displays final output/cropped image
%imshow();

