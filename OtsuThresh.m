function [level, em] = OtsuThresh(datahist, datarange, noisemode)
%OTSUTHRESH data thresholding using Otsu's method. Modified from MATLAB's
%built-in GRAYTHRESH to take a numerical histogram as input rather than an
%image.
%
%   LEVEL = OTSUTHRESH(I) computes a global threshold (LEVEL) that can be
%   used to convert an intensity image to a binary image with IM2BW. LEVEL
%   is a normalized intensity value that lies in the range [0, 1].
%   GRAYTHRESH uses Otsu's method, which chooses the threshold to minimize
%   the intraclass variance of the thresholded black and white pixels.
%
%   [LEVEL, EM] = OTSUTHRESH(I) returns effectiveness metric, EM, as the
%   second output argument. It indicates the effectiveness of thresholding
%   of the input image and it is in the range [0, 1]. The lower bound is
%   attainable only by images having a single gray level, and the upper
%   bound is attainable only by two-valued images.

% Reference:
% N. Otsu, "A Threshold Selection Method from Gray-Level Histograms,"
% IEEE Transactions on Systems, Man, and Cybernetics, vol. 9, no. 1,
% pp. 62-66, 1979.

% Variables names are chosen to be similar to the formulas in
% the Otsu paper.

if ~isempty(datahist)
    p = datahist./sum(datahist);
    
    % filter noise in histogram
    if ~exist('noisemode','var')
        noisemode = 0.005;
    end
    
    if isnumeric(noisemode)
        % use threshold-based noise filtering
        p(p < noisemode) = 0;
    elseif ischar(noisemode) && strcmp(noisemode,'square')
        p = p.^2;
    end
    
    % renormalize PDF after noise filtering
    p = p ./ sum(p);
    
    % compute otsu threshold
    num_bins = length(datahist);
    omega = cumsum(p);
    mu = cumsum(p .* (1:num_bins));
    mu_t = mu(end);
    
    sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));
    
    % Find the location of the maximum value of sigma_b_squared.
    % The maximum may extend over several bins, so average together the
    % locations.  If maxval is NaN, meaning that sigma_b_squared is all NaN,
    % then return 0.
    maxval = max(sigma_b_squared);
    isfinite_maxval = isfinite(maxval);
    if isfinite_maxval
        idx = mean(find(sigma_b_squared == maxval));
        window = floor(idx):(floor(idx)+1);
        level = interp1(window, datarange(window), idx);
    else
        level = nan;
    end
else
    level = nan;
    isfinite_maxval = false;
end

% compute the effectiveness metric
if nargout > 1
    if isfinite_maxval
        em = maxval/(sum(p.*((1:num_bins).^2)) - mu_t^2);
    else
        em = 0;
    end
end


