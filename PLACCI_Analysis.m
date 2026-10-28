%% Analyzing PLACCI fluorescent image: MATLAB version
% Author: Shaoti (Tim) Sui
% Date: Oct. 27, 2026
% Email: timstsui@gmail.com

% This file contain a method to extract cell cycle phase proportion
% information from a snapshot of PLACCI reporter line. The default image
% would be a single-layered image containing three channels. This contains
% a few more steps than directly using FRETENATOR2 in ImageJ, which is more
% optimal for the PLACCI snapshot.

%% Load data.
% Since we used Zeiss image system, we need to open a .czi file, which
% requires a 3rd-party package: BioFormats-MATLAB (bfmatlab). Acquire the
% package from: https://www.openmicroscopy.org/bio-formats/downloads/.
% Add the folder into the path could be done when the folder of the package 
% is in the current shown directory.
addpath(genpath('./'))

% Read .czi file
data = bfopen('Managed/LD-05.czi');
data = data{1,1};

% The first layer contains the general cell label/label for S-G2 phase, the
% H3.1-mCherry
Cell = data{1,1};

% The second layer contains the G2-M label CYCP1;1-YFP
G2M = data{2,1};

% The third layer contains the G1 label CDT1a-CFP
G1 = data{3,1};

% Sanity check
%imagesc(Cell)

%% Denoise
% First we would extract nucleus position from H3.1 layer. This is not to 
% find the exact location, but rather filter out those non-nucleus noise 
% that will influence our result.

% NLM-denoise
% We first use Anscombe transform to stabilize the variance of the image.
% This transform Poisson noise in a microscopy image to Gaussian noise.
Cell_ans = 2 * sqrt(max(double(Cell) + 3/8, 0));

% NLM remove all Gaussian-like noise
Cell_nlm = imnlmfilt(Cell_ans / max(Cell_ans(:)));

%figure
%imshowpair(Cell,Cell_nlm,'montage')

% Remove background trend
A = ones(32);
A = A./sum(A,"all");
Cell_bg = imfilter(Cell_nlm, A, 'replicate');
Cell_debg = Cell_nlm - Cell_bg;

% Do Gaussian filtration
Cell_gaus = imgaussfilt(Cell_debg,4);
Cell_gaus = max(Cell_gaus, 0);
Cell_gaus = Cell_gaus./max(Cell_gaus(:));

%% Fill the holes
% Since H3.1-mCherry is used to label cells, a bad thing is about their
% distribution in the nucleus: H3.1 cannot be found in nucleolus, but the
% meristem cells are filled with nucleolus. To fix this, we need to perform
% the following codes to fill the holes.
% An alternative is to fill the holes over the Otsu-thresholded image, but
% using the image with intensity data would benefit the watershed algorithm
% a lot, so we do not apply Otsu thresholding, and instead use gray image
% to fill holes. This then require us to slice the gray image to multiple
% layer of BW image, each containing pixel with intensity < some level.
% This part is contributed by Claude largely.

lyr = 64; % Given that the largest intensity is ~64 in original picture
nucolusFrac = 0.8; % The approximate portion of nucleolus per nucleus
r = 8; % The approximate radius of a nucleus, in pixel
nucArea = pi*r^2;

Cell_gaus_lyr  = round(Cell_gaus .* (lyr - 1));
nucolusArea = round(nucolusFrac * nucArea);
Cell_close = zeros(size(Cell_gaus_lyr));
for lev = 1:lyr - 1
    lyrmap = Cell_gaus_lyr >= lev;
    if ~any(lyrmap(:)) % Skip if there is nothing in a layer
        break;
    end
    smallDark = ~lyrmap & ~bwareaopen(~lyrmap, nucolusArea);
    Cell_close = Cell_close + double(lyrmap | smallDark);
end
Cell_close = Cell_close ./ max(Cell_close(:));
Cell_close = imclose(Cell_close, strel('disk', 4));
Cell_close = imgaussfilt(Cell_close,1.5);
%Cell_close = Cell_close .* bwareaopen(~(Cell_close <= 0.3*mean(Cell_close(:))), round(0.2*nucArea));


%% Watershed segmentation to acquire single nucleus

Cell_bw = imbinarize(Cell_close, graythresh(Cell_close)/1.2);
Cell_bw = bwareaopen(Cell_bw, round(0.2*nucArea));
Cell_dist = bwdist(~Cell_bw);
Cell_shed = imhmin(-Cell_close.*Cell_bw, 0.02);
Cell_label = watershed(Cell_shed);
Cell_label(~Cell_bw) = 0;

overlayfilt = labeloverlay(Cell_gaus,Cell_label);

% As a benchmark...
Cell_gaus_label = watershed(-Cell_gaus);
overlayfilt_bm = labeloverlay(Cell_gaus, Cell_gaus_label);

%% Sanity check
figure; tiledlayout(1,2)
nexttile
imagesc(Cell_gaus)
colormap("hot")

nexttile
imshow(overlayfilt)

%nexttile
%imshow(overlayfilt_bm)
FigEle

%% Extract the intensity on each layer

