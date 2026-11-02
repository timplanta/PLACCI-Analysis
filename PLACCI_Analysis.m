%% Analyzing PLACCI fluorescent image: MATLAB version
% Author: Shaoti (Tim) Sui
% Date: Oct. 27, 2026
% Email: timstsui@gmail.com

% This file contain a method to extract cell cycle phase proportion
% information from a snapshot of PLACCI reporter line. The default image
% would be a single-layered image containing three channels. This contains
% a few more steps than directly using FRETENATOR2 in ImageJ, which is more
% optimal for the PLACCI snapshot.

% Parameter structure: see the explanation on the default parameter setup
% part of the script.

function Idx = PLACCI_Analysis(path, param)

if nargin == 1
    param = struct(); % Use default parameter first!!
    param.GausStd = 4; % Standard deviation of Gaussian filtration on H3.1 layer after denoising step. This is quite optimum now.
    param.holeFillLayer = 64; % Amount of layers to separate the original grey-scale image into. For hole filling. In my case, this is 64 since that is around the maximum intensity in the picture.
    param.nucolusFrac = 0.8; % The areal proportion of nucleolus (void in every red nucleus) to the whole nucleus.
    param.cellRadius = 8; % The approximate radius of each cell on the picture, in pixel.
    param.H = 0.02; % The depth level tolerated in each local minima. Used for h-min algorithm to smooth out noisy local minima before watershed. See following MATLAB document for details: https://www.mathworks.com/help/images/ref/imhmin.html#buo13ca-1-h
end

%% Load data.
% Since we used Zeiss image system, we need to open a .czi file, which
% requires a 3rd-party package: BioFormats-MATLAB (bfmatlab). Acquire the
% package from: https://www.openmicroscopy.org/bio-formats/downloads/.
% Add the folder into the path could be done when the folder of the package 
% is in the current shown directory.

%addpath(genpath('./'))
%path = 'Managed/LD-01.czi';

% Read .czi file
echo off
data = bfopen(path);
data = data{1,1};
echo on;

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
Cell_gaus = imgaussfilt(Cell_debg,param.GausStd);
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

lyr = param.holeFillLayer; % Given that the largest intensity is ~64 in original picture
nucolusFrac = param.nucolusFrac; % The approximate portion of nucleolus per nucleus
r = param.cellRadius; % The approximate radius of a nucleus, in pixel
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
Cell_close = imclose(Cell_close, strel('disk', param.cellRadius/2));
Cell_close = imgaussfilt(Cell_close,1.5); % Double filtration to eliminate some additional noise.
%Cell_close = Cell_close .* bwareaopen(~(Cell_close <= 0.3*mean(Cell_close(:))), round(0.2*nucArea));


%% Watershed segmentation to acquire single nucleus

Cell_bw = imbinarize(Cell_close, graythresh(Cell_close)/1.2);
Cell_bw = bwareaopen(Cell_bw, round(0.2*nucArea));
Cell_shed = imhmin(-Cell_close.*Cell_bw, param.H);
Cell_label = watershed(Cell_shed);
Cell_label(~Cell_bw) = 0;

%overlayfilt = labeloverlay(Cell_gaus,Cell_label);

% As a benchmark...
%Cell_gaus_label = watershed(-Cell_gaus);
%overlayfilt_bm = labeloverlay(Cell_gaus, Cell_gaus_label);

%% Sanity check
%figure; tiledlayout(1,2)
%nexttile
%imagesc(Cell_gaus)
%colormap("hot")

%nexttile
%imshow(overlayfilt)

%nexttile
%imshow(overlayfilt_bm)
%FigEle

%% Extract the intensity on each layer

% With a relatively accurate cell labeling, we now proceed to other layers
% and see the cell cycle phase of each cell.

% Measure mean fluorescence for each segmented nucleus across channels.
% Notice that for H3.1-mCherry channel, we have to reduce the intensity to
% see which cell has high, 2x intensity; for other sensors, we would first
% normalize them.

G2M_ans = 2 * sqrt(max(double(G2M) + 3/8, 0));
G2M_debg = imnlmfilt(G2M_ans / max(G2M_ans(:)));

G1_ans = 2 * sqrt(max(double(G1) + 3/8, 0));
G1_debg = imnlmfilt(G1_ans / max(G1_ans(:)));

Cell_bg = imgaussfilt(Cell_nlm, 64);
G2M_bg = imgaussfilt(G2M_debg, 64);
G1_bg = imgaussfilt(G1_debg, 64);

CellStats = findMeanIntens(Cell_label, Cell_nlm);
CellNoise = findMeanIntens(Cell_label, Cell_bg);
G2MStats = findMeanIntens(Cell_label, G2M_debg);
G2MNoise = findMeanIntens(Cell_label, G2M_bg);
G1Stats = findMeanIntens(Cell_label, G1_debg);
G1Noise = findMeanIntens(Cell_label, G1_bg);

% We would use the signal-noise ratio to compare across channel and find
% which channel contains signal.
SNRMat = [CellStats, G2MStats, G1Stats]./[CellNoise, G2MNoise, G1Noise];

SNRz = (log(SNRMat))./std(log(SNRMat),[],2);

[~, Idx] = max(SNRz,[], 2);
for i = 1:length(Idx)
    if Idx(i) == 1
        if SNRz(i,2) >= 1
            Idx(i) = 2;
        elseif SNRz(i,3) >= 1
            Idx(i) = 3;
        end
    end
end

%% False-color plot of PLACCI result

Idx_map = [0; Idx];
Cell_phase = Idx_map(Cell_label+1);
colorbar = [219 33 0
            245 221 10
            104 224 237] ./ 255;

result = labeloverlay(Cell_gaus, Cell_phase, Colormap=colorbar);

%figure;
%imshow(result)
%FigEle

[filepath,name,~] = fileparts(path);
if ~exist([filepath, '/Digitalized_PLACCI/'], 'dir')
    mkdir([filepath, '/Digitalized_PLACCI/']);
end

imwrite(result, [filepath, '/Digitalized_PLACCI/', name, '.tif'], 'tif')



%%

%figure
%tiledlayout(1,3)
%nexttile
%histogram(SNRz(:,1),50)
%nexttile
%histogram(SNRz(:,2),50)
%nexttile
%histogram(SNRz(:,3),50)


end
%%
function stat = findMeanIntens(label, map)
    stat = regionprops(label, map, 'MeanIntensity');
    stat = [stat.MeanIntensity]';
end