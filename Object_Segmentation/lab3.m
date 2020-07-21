%addpath lab_library/;
%Add the directory and its subdirectories 
%addpath(genpath('lab_library/'))
function [] = lab3()
% init_method = 'random';
init_method = 'kmeans++';
K = 10;  
scale_factor = 0.5; 
image_sigma = 3.0; 
kmeans_quetions(K,scale_factor ,image_sigma , init_method)

% spatial_bandwidth = 10.0;  % spatial bandwidth
% colour_bandwidth = 5.0;   % colour bandwidth
% % % % % image = 'tiger3.jpg';
% % % % % for i = 2:2:10
% % % % %     for j = 1:1:10
% % % % %         mean_shift_segmentation(image, i, j) 
% % % % %     end 
% % % % % end
% % % norm_cut_segmentation() 
% % colour_bandwidth = 10.0; 
% % radius = 7;              
% % ncuts_thresh = 0.2;      
% % scale_factor = 0.4;     
% % image_sigma = 2.0;    
% % image = 'tiger1.jpg';
% % min_area = 250;        
% % max_depth = 8; 
% % for radius = 3:1:10
% %     norm_cut_segmentation(image,colour_bandwidth , radius,ncuts_thresh, scale_factor, image_sigma, min_area, max_depth) 
% % end

% % for colour_bandwidth = 10:10:100
% %     for radius = 3:2:10
% %         for ncuts_thresh = 0.1:0.1:0.4
% %             for min_area = 100:50:250
% %                 norm_cut_segmentation(image,colour_bandwidth , radius,ncuts_thresh, scale_factor, image_sigma, min_area, max_depth) 
% %             end 
% %         end
% %     end
% % end

% image = 'orange.jpg';

% % % image = 'tiger3.jpg';
% % % alpha = 15;                 % maximum edge cost
% % % sigma = 20;                % edge cost decay factor
% % % K = 16;
% % % for k = 1:1:3
% % %     graphcut_segmentation(image, k, alpha, sigma);
% % % end

% graphcut_segmentation(image, K, alpha, sigma);
% for alpha = 4:2:10
%     for sigma = 4:2:12
%         graphcut_segmentation(image, alpha, sigma);
%     end
% end

end 

function [ segm, prior ] = graphcut_segmentation(image, K, alpha, sigma)
scale_factor = 0.5;          % image downscale factor
area = [ 80, 110, 570, 300 ] % image region to train foreground with
% K = 16;                      % number of mixture components
% alpha = 8.0;                 % maximum edge cost
% sigma = 10.0;                % edge cost decay factor

I = imread(image);
I = imresize(I, scale_factor);
Iback = I;
area = int16(area*scale_factor);
[ segm, prior ] = graphcut_segm(I, area, K, alpha, sigma);


Inew = mean_segments(Iback, segm);
I = overlay_bounds(Iback, segm);
figure;
imwrite(Inew,'graphcut1.png')
imwrite(I,'graphcut2.png')
imwrite(prior,'graphcut3.png')
subplot(2,2,1); imshow(Inew);
subplot(2,2,2); imshow(I);
subplot(2,2,3); imshow(prior);

end 

function [prob]  = mixture_prob(I, K ,maxIter, mask)


[gmm] = EM_algorithm(I, K, maxIter, mask);
mixing_coefficients = gmm.pi
[im_h, im_w, ~] = size(I);
Ivec = single(reshape(I, im_h*im_w, 3)); % Flattened out pixels
eps = 1e-8; % Small number for numerical stability 

prob = zeros(im_h*im_w, K);
for k = 1:K
%     gmm_fg.cov(:,:,k)
    %prob(:,k) = mvnpdf(Ivec1, gmm_fg.mu(k,:), gmm_fg.cov{k}) + eps;
    temp = gmm.cov(:,:,k);
    s = (temp + temp.') / 2;
    prob(:,k) = (gmm.pi(k))  .* mvnpdf(Ivec, gmm.mu(k,:), s) + eps;
end
% logl_fg = bsxfun(@plus, -log(prob), -log(gmm_fg.pi));
% [~,I] = min(logl_fg, [], 2);
% 
% prob = 0;
prob = sum(prob,2);
probb_examples = prob(1:10,:)
end 


function [gmm] = EM_algorithm(I, K, L, mask)
% Learn GMM components for a set of pixels with the EM algorithm
% Args: I = All pixels (RGB values)
%       K = Number of clusters
%       L = Number of iterations for EM 
%       T = Trimap with pixel indeces
% Returns: gmm = GMM components (means, covariance matrix, weighting coeff.)
%%%
[im_h, im_w] = size(mask);
T = reshape(mask, im_h*im_w, 1);
T = find(T);
% Store all pixels from trimap
[im_h, im_w, ~] = size(I);
Ivec = single(reshape(I, im_h*im_w, 3));
IvecCCC_size = size(Ivec)
T_size = size(T)
Ivec1 = Ivec(T,:);
N = size(Ivec1,1);
D = size(Ivec1,2);


% disp("Ivec : ")
% size(Ivec1)  % Covariance matrix 

responsability = zeros(N,K);
eps = 1e-8;

% Initialize weights, mean and covariance with Kmeans
[segm, centers] = kmeans(Ivec1, K, 'Distance', 'cityblock', 'Replicates', 5);
disp("Mean : ")
mu_k = centers;  % mean 
disp("Sigma : ")
%Sigma = cell(K,1)  % Covariance matrix  zeros(2,2,3)
Sigma = zeros(3,3,K);
disp("Mixing coeff. : ")
w_k = ones(1,K)/K; % mixing coefficients 
for k = 1:K
    %Sigma{k} = cov(Ivec1((segm(:) == k), :));
    Sigma(:,:,k) = single(cov(Ivec1((segm(:) == k), :)))
    w_k(k) = sum((segm(:) == k)) / numel(segm)
end
disp("Sigma 2 : ")
% Sigma  % Covariance matrix 

% EM algorithm
fprintf('\n EM algorithm \n')
for l = 1:L

    % E-step. Compute the responsibilities from Eq. 18.
    for k = 1:K
        size_sigmaaaaaa = size(Sigma(:,:,k))
        sigmund = Sigma(:,:,k) + eps
        kkkkkkkkkkk = k
%         ivec1111111 = Ivec1
%         mu_kkkkkkkkk = mu_k(k,:)
        temp = Sigma(:,:,k);
        s = (temp + temp.') / 2
        p = mvnpdf(Ivec1,mu_k(k,:), s);
        p = w_k(k)*p ;
        responsability(:,k) = p ;
    end
    
    %Normalize 
    for i = 1:N
        responsability(i,:) = responsability(i,:) / (sum(responsability(i,:)) + eps );
    end
    % M-step. Compute the means, covariances and mixture weights from Eq. 19, 20, and 21.

    % Update mixing coefficients
    w_k = mean(responsability,1);
    % Update means
    for k = 1:K
        mu_k(k,:) = sum(Ivec1.*repmat(responsability(:,k),1,D),1)./ (sum(responsability(:,k)) + eps );
    end
    % update covariances
    for k = 1:K
%         size(repmat(mu_k(k,:),N,1))
        k_updateeeeeeeeeeeeeeeeeeee = k
        Xm = Ivec1 - repmat(mu_k(k,:),N,1);
        temp = (Xm.*repmat(responsability(:,k),1,D))'*Xm
        responsabilityyyyyyyyyyyyyyyyyyyy = sum(responsability(:,k))
        Sigma(:,:,k) = temp ./ ( sum(responsability(:,k)) + eps )  ;
        Sigma(:,:,k) = Sigma(:,:,k) + eps;
    end
    
%     fprintf('Loop %i \n', l)
end
fprintf('EM done \n')

%%% Store GMMs in struct object
gmm.mu = mu_k;
gmm.cov = Sigma;
gmm.pi = w_k;

end



function [ segm, prior ] = graphcut_segm(I, area, K, alpha, sigma)
% area = [ minx, miny, maxx, maxy ]

[h,w,c] = size(I);
dw = area(3) - area(1) + 1;
dh = area(4) - area(2) + 1;
mask = uint8([zeros(area(2)-1,w); zeros(dh,area(1)-1), ones(dh,dw), ...
	     zeros(dh,w-area(3)); zeros(h-area(4),w)]);
Ivec = single(reshape(I, size(I,1)*size(I,2), 3));

grey = single(rgb2gray(I));
h = fspecial('gaussian', [7,7], 0.5);
grey = imfilter(grey, h);
h = fspecial('sobel');
dx = imfilter(grey, h/4);
dy = imfilter(grey, h/4');
grad = sqrt(dx.^2 + dy.^2);
edge = (alpha*sigma)*ones(size(grey)) ./ (grad + sigma);

Ivec_size = size(Ivec)
I_size = size(I)
mask_size = size(mask)
mask_inv = 1 - mask;
mask_inverse_size = size(mask_inv)
L = 3;
prob = mixture_prob(I, K ,L, mask);
% [gmm] = EM_algorithm(I, K, L, mask);

tic
for l=1:3
    iterationss = l
    fprintf('Find Gaussian mixture models...\n');
    fprob = mixture_prob(I, K, 10, mask);
    bprob = mixture_prob(I, K, 10, 1-mask);
    prior = reshape(fprob ./ (fprob + bprob), size(I,1), size(I,2), 1);
    toc

    fprintf('Find minimum cut...\n');
    [u, erriter, i] = cmf_cut(prior, edge);
    mask = uint8(u>0.5);
    toc

end

segm = int16(u>0.5) + 1;
end


function [u, erriter, i] = cmf_cut(ur, alpha)

[rows, cols] = size(ur);
imgSize = rows*cols;

% define the required parameters:
%
%   - alpha: the penalty parameter to the total-variation term.
%       For the case without incorporating image-edge weights, alpha is given
%       by the constant everywhere. For the case with image-edge weights,
%       alpha is given by the pixelwise weight function:
%
%       For example, alpha(x) = b/(1 + a*| nabla f(x)|) where b and a are positive
%       constants and |nabla f(x)| gives the strength of the local gradient.
%
%   - cc: gives the step-size of the augmented Lagrangian method.
%       The optimal range of cc is [0.3, 3].
%
%   - errbound: the error bound for convergence.
%
%   - numIter: the maximum iteration number.
%
%   - steps: the step-size for the graident-projection step to the
%       total-variation function. The optimal range of steps is [0.1,
%       0.17].
%

% alpha = alphafac*ones(rows,cols); 
cc = 0.3;
errbound = 1e-4;
numIter = 100;
steps = 0.16;

% build up the data terms
ulab(1) = 0.0; %0.15;
ulab(2) = 1.0; %0.6;
Cs = abs(ur - ulab(1));
Ct = abs(ur - ulab(2));

urb = 0.8*ur + 0.1;
Cs = -log(1.0 - urb);
Ct = -log(urb);


% set the initial values
%   - the initial value of u is set to be an initial cut, see below.
%   - the initial values of two terminal flows ps and pt are set to be the
%     specified legal flows.
%   - the initial value of the spatial flow fiels p = (pp1, pp2) is set to
%   be zero.

u = double((Cs-Ct) >= 0);
ps = min(Cs, Ct);
pt = ps;

pp1 = zeros(rows, cols+1);
pp2 = zeros(rows+1, cols);
divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);

erriter = zeros(numIter,1);

for i = 1:numIter

	% update the spatial flow field p = (pp1, pp2):
    %   the following steps are the gradient descent step with steps as the
    %   step-size.
    
    pts = divp - (ps - pt  + u/cc);
    pp1(:,2:cols) = pp1(:,2:cols) + steps*(pts(:,2:cols) - pts(:,1:cols-1)); 
    pp2(2:rows,:) = pp2(2:rows,:) + steps*(pts(2:rows,:) - pts(1:rows-1,:));
    
    % the following steps give the projection to make |p(x)| <= alpha(x)
    
    gk = sqrt((pp1(:,1:cols).^2 + pp1(:,2:cols+1).^2 + pp2(1:rows,:).^2 + pp2(2:rows+1,:).^2)*0.5);
    gk = double(gk <= alpha) + double(~(gk <= alpha)).*(gk ./ alpha);
    gk = 1 ./ gk;
    
    pp1(:,2:cols) = (0.5*(gk(:,2:cols) + gk(:,1:cols-1))).*pp1(:,2:cols); 
    pp2(2:rows,:) = (0.5*(gk(2:rows,:) + gk(1:rows-1,:))).*pp2(2:rows,:);
    
    divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);
    
    % updata the source flow ps
    
    pts = divp + pt - u/cc + 1/cc;
    ps = min(pts, Cs);
    
    % update the sink flow pt
    
    pts = - divp + ps + u/cc;
    pt = min(pts, Ct);

	% update the multiplier u
    
	erru = cc*(divp + pt  - ps);
	u = u - erru;
    
    % evaluate the avarage error
    
    erriter(i) = sum(sum(abs(erru)))/imgSize; 
   
    if (erriter(i) < errbound)
        break;
    end
end

msg = sprintf('number of iterations = %u. \n', i);
disp(msg);

end 

function [] = norm_cut_segmentation(image,colour_bandwidth , radius,ncuts_thresh, scale_factor, image_sigma, min_area, max_depth) 
% colour_bandwidth = 20.0; % color bandwidth
% radius = 3;              % maximum neighbourhood distance
% ncuts_thresh = 0.2;      % c
% utting threshold
% min_area = 200;          % minimum area of segment
% max_depth = 8;           % maximum splitting depth
% scale_factor = 0.4;      % image downscale factor
% image_sigma = 2.0;       % image preblurring scale
% I = imread('tiger2.jpg');
I = imread(image);
I = imresize(I, scale_factor);
Iback = I;
d = 2*ceil(image_sigma*2) + 1;
h = fspecial('gaussian', [d d], image_sigma);
I = imfilter(I, h);

segm = norm_cuts_segm(I, colour_bandwidth, radius, ncuts_thresh, min_area, max_depth);
Inew = mean_segments(Iback, segm);
I = overlay_bounds(Iback, segm);
imwrite(Inew,'normcuts1.png')
imwrite(I,'normcuts2.png')
figure;
% subplot(1,2,1); imshow(Inew);
% subplot(1,2,2); imshow(I);
imshow(I);
title(sprintf('r = %d, col_b = %d, ncut = %d, area = %d,depth = %d ',radius, colour_bandwidth, ncuts_thresh, min_area, max_depth ))

end


function [EV, EVal] = ncuts(A, n_ev)
% Computes the n_ev smallest (non-zero) eigenvectors and eigenvalues of the 
% of the Laplacian of A

D = sparse(1:size(A,1), 1:size(A,1), full(sum(A, 1)), size(A,1), size(A,2));

opts.issym = 0;
opts.isreal = 1;
opts.disp = 0;
nvec = n_ev+1;

[EV, EVal] = eigs((D - A) + (10^-10) * speye(size(D)), D, nvec, 'sm',opts);

[junk, sortidx] = sort(diag(EVal), 'descend');
EV = EV(:,sortidx(end-1:-1:1));
v = diag(EVal);
EVal = v(sortidx(end-1:-1:1));

EV = bsxfun(@rdivide, EV, sqrt(sum(EV.^2,1))); % makes the eigenvectors unit norm

end 

function [A] = ncuts_affinity(im, XY_RADIUS, RGB_SIGMA)

sz = [size(im,1), size(im,2)];

% Find all pairs of pixels within a distance of XY_RADIUS
rad = ceil(XY_RADIUS);
[di,dj] = ndgrid(-rad:rad, -rad:rad);
dv = (dj.^2 + di.^2) <= XY_RADIUS.^2;
di = di(dv);
dj = dj(dv);

[i,j] = ndgrid(1:size(im,1), 1:size(im,2));
i = repmat(i(:), 1, length(di));
j = repmat(j(:), 1, length(di));
i_ = bsxfun(@plus, i, di');
j_ = bsxfun(@plus, j, dj');
v = (i_ >= 1) & (i_ <= size(im,1)) & (j_ >= 1) & (j_ <= size(im,2));
pair_i = sub2ind(sz, i(v), j(v));
pair_j = sub2ind(sz, i_(v), j_(v));

% Weight each pair by the difference in RGB values, divided by RGB_SIGMA
RGB = double(reshape(im, [], size(im,3)))/RGB_SIGMA;
W = exp(-sum((RGB(pair_i,:) - RGB(pair_j,:)).^2,2));

% Construct an affinity matrix
A = sparse(pair_i, pair_j, W, prod(sz), prod(sz));
end

function segm = norm_cuts_segm(I, colour_bandwidth, radius, ncuts_thresh, max_area, max_depth)
tic
[nRow, nCol, c] = size(I);
N = nRow * nCol;
V = reshape(I, N, c); % connect up-to-down way. Vertices of Graph

fprintf('Compute affinity matrix...\n');
W = ncuts_affinity(I, radius, colour_bandwidth);
toc 

fprintf('Solve eigenvalue problems to find partitions...\n');
Seg = (1:N)'; 
[Seg, Id, Ncut] = ncuts_partition(Seg, W, ncuts_thresh, max_area, 'ROOT', max_depth, 1);

segm = zeros(nRow*nCol, 1);
for i=1:length(Seg)
    segm(Seg{i}) = i;
    fprintf('%s. Ncut = %f\n', Id{i}, Ncut{i});
end

segm = uint32(reshape(segm, nRow, nCol));
toc

end 


function [Seg, Id, Ncut] = ncuts_partition(I, W, sNcut, sArea, id, maxDepth, depth)
% NcutPartition - Partitioning
%
% Synopsis
%  [sub ids ncuts] = ncuts_partition(I, W, sNcut, sArea, [id])
%
% Description
%  Partitioning. This function is called recursively.
%
% Inputs ([]s are optional)
%  (vector) I        N x 1 vector representing a segment to be partitioned.
%                    Each element has a node index of V (global segment).
%  (matrux) W        N x N matrix representing the computed similarity
%                    (weight) matrix.
%                    W(i,j) is similarity between node i and j.
%  (scalar) sNcut    The smallest Ncut value (threshold) to keep partitioning.
%  (scalar) sArea    The smallest size of area (threshold) to be accepted
%                    as a segment.
%  (string) [id]     A label of the segment (for debugg)
%
% Outputs ([]s are optional)
%  (cell)   Seg      A cell array of segments partitioned.
%                    Each cell is the each segment.
%  (cell)   Id       A cell array of strings representing labels of each segment.
%                    IDs are generated as children based on a parent id.
%  (cell)   Ncut     A cell array of scalars representing Ncut values
%                    of each segment.
%
% Requirements
%  NcutValue
%
% Authors
%  Naotoshi Seo <sonots(at)sonots.com>
%
% License
%  The program is free to use for non-commercial academic purposes,
%  but for course works, you must understand what is going inside to use.
%  The program can be used, modified, or re-distributed for any purposes
%  if you or one of your group understand codes (the one must come to
%  court if court cases occur.) Please contact the authors if you are
%  interested in using the program without meeting the above conditions.

% Changes
%  10/01/2006  First Edition
% Compute D
N = length(W);
d = sum(W, 2);
D = spdiags(d, 0, N, N); % diagonal matrix

% Step 2 and 3. Solve generalized eigensystem (D -W)*S = S*D*U (12).
% (13) is not necessary thanks to smart matlab.
% Get the 2 smallests ('sm')
warning off; % let me stop warning
[EV, EVal] = ncuts(W, 2); %%%
% 2nd smallest (1st smallest has all same value elements, and useless)
U2 = EV(:, 2);

% Step 3. Refer 3.1 Example 3.
% Bipartition the graph at point that Ncut is minimized.
t = mean(U2);
% ncut = ncuts_value(t, U2, W, D)
% t = fminsearch('ncuts_value', t, [], U2, W, D);
t = fminsearch(@ncuts_value, t, [], U2, W, D);
A = find(U2 > t);
B = find(U2 <= t);

% Step 4. Decide if the current partition should be divided
% if either of partition is too small, stop recursion.
% if Ncut is larger than threshold, stop recursion.
ncut = ncuts_value(t, U2, W, D);
fprintf('Cutting ncut=%.3f sizes=(%d,%d)\n', ncut, length(A), length(B));
if (length(A) < sArea || length(B) < sArea) || ncut > sNcut || depth > maxDepth
    Seg{1}   = I;
    Id{1}   = id; % for debugging
    Ncut{1} = ncut; % for duebuggin
    return;
end

% Seg segments of A
[SegA IdA NcutA] = ncuts_partition(I(A), W(A, A), sNcut, sArea, [id '-A'], maxDepth, depth+1);
% I(A): node index at V. A is index at the segment, I
% W(A, A); % weight matrix in segment A

% Seg segments of B
[SegB IdB NcutB] = ncuts_partition(I(B), W(B, B), sNcut, sArea, [id '-B'], maxDepth, depth+1);

% concatenate cell arrays
Seg   = [SegA SegB];
Id   = [IdA IdB];
Ncut = [NcutA NcutB];
end

% ncuts_value
% ncuts_value
function [ncut] = ncuts_value(t, U2, W, D)
% NcutValue - 2.1 Computing the Optimal Partition Ncut. eq (5)
%
% Synopsis
%  ncut = ncuts_value(T, U2, D, W);
%
% Inputs ([]s are optional)
%  (scalar) t        splitting point (threshold)
%  (vector) U2       N x 1 vector representing the 2nd smallest
%                     eigenvector computed at step 2.
%  (matrix) W        N x N weight matrix
%  (matrix) D        N x N diagonal matrix
%
% Outputs ([]s are optional)
%  (scalar) ncut     The value calculated at the right term of eq (5).
%                    This is used to find minimum Ncut.
%
% Authors
%  Naotoshi Seo <sonots(at)sonots.com>
%
% License
%  The program is free to use for non-commercial academic purposes,
%  but for course works, you must understand what is going inside to use.
%  The program can be used, modified, or re-distributed for any purposes
%  if you or one of your group understand codes (the one must come to
%  court if court cases occur.) Please contact the authors if you are
%  interested in using the program without meeting the above conditions.

% Changes
%  10/01/2006  First Edition
x = (U2 > t);
x = (2 * x) - 1; % convert [1 0 0 1 0]' to [1 -1 -1 1 -1]' to follow paper's way
d = diag(D);
k = sum(d(x > 0)) / sum(d);
b = k / (1 - k);
y = (1 + x) - b * (1 - x);
ncut = (y' * (D - W) * y) / ( y' * D * y );
end


function [] = mean_shift_segmentation(image, spatial_bandwidth, colour_bandwidth) 
scale_factor = 0.5;       % image downscale factor
% spatial_bandwidth = 10.0;  % spatial bandwidth
% colour_bandwidth = 5.0;   % colour bandwidth
num_iterations = 40;      % number of mean-shift iterations
image_sigma = 2.0;        % image preblurring scale
% I = imread('orange.jpg');
I = imread(image);
% I = imread('tiger1.jpg');
I = imresize(I, scale_factor);
Iback = I;
d = 2*ceil(image_sigma*2) + 1;
h = fspecial('gaussian', [d d], image_sigma);
I = imfilter(I, h);

segm = mean_shift_segm(I, spatial_bandwidth, colour_bandwidth, num_iterations);
Inew = mean_segments(Iback, segm);
I = overlay_bounds(Iback, segm);
% imwrite(Inew,'meanshift1.png')
% imwrite(I,'meanshift2.png')
% subplot(1,2,1); imshow(Inew);
% subplot(1,2,2); imshow(I);
figure;
imshow(I);
title(sprintf('spatial_b = %d, col_b = %d ',spatial_bandwidth, colour_bandwidth ))

end 

%%%

function segm = mean_shift_segm(I, spatial_bandwidth, colour_bandwidth, num_iterations)
tic
fprintf('Find colour channels with K-means...\n');
K = 16;
% kmeans_segm(image, K, maxIteration, seed, init_method)
[ segm, centers ] = kmeans_segm(I, K, 10, 4321, 'kmeans');
toc

centers(isnan(centers)) = 0.0;
%imshow(overlay_bounds(I, segm))
%pause

[ height, width, depth ] = size(I);
idx = reshape(segm, [height, width]);
maps = zeros(height, width, K, 'single');
mapx = zeros(height, width, K, 'single');
mapy = zeros(height, width, K, 'single');
[X, Y] = meshgrid(1:width, 1: height);
for k = 1:K
    maps(:,:,k) = (idx == k);
    mapx(:,:,k) = maps(:,:,k).*X;
    mapy(:,:,k) = maps(:,:,k).*Y; 
end
s = 2*ceil(1.5*spatial_bandwidth) + 1;
h = fspecial('gaussian', [s s], spatial_bandwidth);
mapsw = reshape(imfilter(maps, h), [width*height,K]) + 1e-6;
mapsx = reshape(imfilter(mapx, h), [width*height,K]);
mapsy = reshape(imfilter(mapy, h), [width*height,K]);
toc

fprintf('Search for high density points...\n');
constC = -0.5/(colour_bandwidth^2);
x = reshape(X, [width*height, 1]);
y = reshape(Y, [width*height, 1]);
Ic = single(reshape(I, [width*height, 3]));
wei = exp(constC*pdist2(Ic, centers));
for l = 1:num_iterations
    p = (round(x)-1)*height + round(y);
    ww = mapsw(p,:) .* wei;
    w = sum(ww, 2);
    u = bsxfun(@rdivide, ww*centers, w);
    x = bsxfun(@rdivide, sum(mapsx(p,:).*wei, 2), w);
    y = bsxfun(@rdivide, sum(mapsy(p,:).*wei, 2), w);
    wei = bsxfun(@rdivide, ww, w);
    x = max(min(x, width), 1);
    y = max(min(y, height), 1);
    toc
end

fprintf('Assign high density points to pixels...\n');
XY = [ x, y ];
thr = 10.0;
val = 0;
mask = zeros(height*width, 1, 'int16');
for y=1:height
    for x=1:width
        p = round(x-1)*height + round(y);
        if mask(p)==0
    	    stack = [ p ];
            val = val + 1;
            mask(p) = val;
            while length(stack)>0
	        p0 = stack(end);
                xy = XY(p0,:);
                x0 = floor((p0-1)/height) + 1;
                y0 = mod(p0-1, height) + 1;
	        stack = stack(1:end-1);
                if x0<width & mask(p0+height)==0 & sum((xy-XY(p0+height,:)).^2)<thr
                    stack = [ stack, p0+height ]; 
                    mask(p0+height) = val;
                end
                if x0>1 & mask(p0-height)==0 & sum((xy-XY(p0-height,:)).^2)<thr
                    stack = [ stack, p0-height ]; 
                    mask(p0-height) = val;
                end
                if y0<height & mask(p0+1)==0 & sum((xy-XY(p0+1,:)).^2)<thr
                    stack = [ stack, p0+1 ]; 
                    mask(p0+1) = val;
                end
                if y0>1 & mask(p0-1)==0 & sum((xy-XY(p0-1,:)).^2)<thr
                    stack = [ stack, p0-1 ]; 
                    mask(p0-1) = val;
                end
            end
        end
    end
end
segm = reshape(mask, [height,width]);
toc
end 



%%%%%%%%%%%%%%%%%%%%%% TODO  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [dist] = distance(p1, p2)
% Function to compute euclidean distance 
    dist = sum((p1 - p2).^2);
end


function [centroids] = kmeans_pp_init(I, K)

    I = im2double(I);
    F = reshape(I,size(I,1)*size(I,2),3);  % Color Features
    centroids = zeros(K,3);
    centroids(1,:)  = F(ceil(rand()*size(F,1)), :);
    f = realmax;
    for k = 1:K-1
        dist = ones(size(F,1),1).* (-f);
        for i = 1:size(F,1)
            point = F(i,:);
            d = realmax; 
            for j = 1:k
                temp_dist = distance(point, centroids(j,:));
                d = min(d, temp_dist);
            end
            dist(i,:) = d;      
        end
%         dist_size = size(dist)
        [maximum, index] = max(dist);
%         next_centroid = F(index,:) 
%         centroidddddds = centroids
        centroids(k+1,:) = F(index,:) ;
    end
    
end


function [] = kmeans_quetions(K,scale_factor ,image_sigma , init_method) 

% clear command windows
% clc
% clear all
% close all

% K = 8;               % number of clusters used
L = 50;              % number of iterations
seed = 14;           % seed used for random initialization
% scale_factor = 1.0;  % image downscale factor
% image_sigma = 1.0;   % image preblurring scale
% I = imread('tiger3.jpg');
% I = imread('tiger1.jpg');
I = imread('orange.jpg');
I = imresize(I, scale_factor);
Iback = I;
d = 2*ceil(image_sigma*2) + 1;
h = fspecial('gaussian', [d d], image_sigma);
I = imfilter(I, h);

tic
[ segm, centers ] = kmeans_segm(I, K, L, seed, init_method);
toc

% size_Img_begin = size(Iback)
% size_Img_segm = size(segm)
Inew = mean_segments(Iback, segm);
I = overlay_bounds(Iback, segm);
imwrite(Inew,'kmeans1.png')
imwrite(I,'kmeans2.png')
figure;
subplot(1,2,1); imshow(Inew);
subplot(1,2,2); imshow(I);

end



function [segm, centroids] = kmeans_segm(image, K, maxIteration, seed, init_method)
rng(seed)
% K-means Segmentation (option: K (Number of Clusters))
% I = imread(image); 
I = im2double(image);
F = reshape(I,size(I,1)*size(I,2),3);                 % Color Features
% K-means
% disp("Random centroids : ")
% CENTS = F( ceil(rand(K,1)*size(F,1)) ,:)            % Cluster Centers
% disp("Kmeans++ centroids : ")
% centroids = kmeans_pp_init(I, K)
disp("K Centroids : ")
tic
if init_method(1) == 'k'
    centroids = kmeans_pp_init(I, K)
else
    centroids = F( ceil(rand(K,1)*size(F,1)) ,:)
end
toc 
eps = 1e-5;
distanceLabels   = zeros(size(F,1),K+2);                         % Distances and Labels
for n = 1:maxIteration
   for i = 1:size(F,1)
      for j = 1:K  
        distanceLabels(i,j) = norm(F(i,:) - centroids(j,:));      
      end
      [Distance, CN] = min(distanceLabels(i,1:K));                
      distanceLabels(i,K+1) = CN;                               
      distanceLabels(i,K+2) = Distance;                         
   end
   centroids_before = centroids;
   for i = 1:K
      A = (distanceLabels(:,K+1) == i);                       
      centroids(i,:) = mean(F(A,:));                    
      if sum(isnan(centroids(:))) ~= 0                   
         NC = find(isnan(centroids(:,1)) == 1);          
         for Ind = 1:size(NC,1)
            centroids(NC(Ind),:) = F(randi(size(F,1)),:);
         end
      end
   end
   centroids_after = centroids;
   S = centroids_after - centroids_before;
   S = max(S(:));
   Tobreak = S < eps;
    if Tobreak
        disp("Break Converged !! ")
        disp(n)
        break; 
    end
end

X = zeros(size(F));
for i = 1:K
    idx = find(distanceLabels(:,K+1) == i);
    X(idx,:) = repmat(centroids(i,:),size(idx,1),1); 
end
output_image = reshape(X,size(I,1),size(I,2),3);

figure;
imshow(output_image);  
title(['Kmeans',' : ',num2str(K)]);


X = zeros(size(F,1),1);
for i = 1:K
    idx = find(distanceLabels(:,K+1) == i);
    X(idx,:) = repmat(i,size(idx,1),1); 
end
segm = reshape(X,size(I,1),size(I,2),1);

end 


function [imout] = overlay_bounds(im, segm)
[h,w] = size(segm);
imcut = segm(2:h-1,2:w-1);
diff = [ ones(1,w) ; ones(h-2,1), ((imcut==segm(3:h,2:w-1)) & (imcut==segm(2:h-1,3:w)) & (imcut==segm(1:h-2,2:w-1)) & (imcut==segm(2:h-1,1:w-2))), ones(h-2,1); ones(1,w) ];
mask = repmat(uint8(diff), [1,1,3]);
imout = (im .* mask);
imout(:,:,1) = imout(:,:,1) + (1 - mask(:,:,1))*255;

end 

function [Inew] = mean_segments(I, segm)

N = max(max(segm));
[h, w, c] = size(I);
Ic = single(reshape(I, [h*w, c]));
sc = reshape(segm, [h*w,1]);
cols = zeros(N, c);
nums = zeros(N, 1);
for i=1:h*w
    s = sc(i);
    cols(s,:) = cols(s,:) + Ic(i,:);
    nums(s) = nums(s) + 1;
end
cols = bsxfun(@rdivide, cols, nums);
Inew = uint8(reshape(cols(sc,:), [h,w,c]));

end 
