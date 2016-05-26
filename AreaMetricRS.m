function [a, onMean, onMed, onGeoMean] = ...
    AreaMetricRS(qhist, bhist, datarange)
% AreaMetricRS computes the fraction of a distribution that is above and
% outside a reference "basal" distribution, after rescaling left half of basal
% distribution to optimally overlay with query distribution.
%
% Updated 20160328
if any(isnan(bhist))
    a = nan;
    return
end

bhist = reshape(bhist,[],1);
qhist = reshape(qhist,[],1);

idxmean = floor(sum(bhist.*(1:length(bhist))')./sum(bhist));

% scaled off distribution
yref = bhist(1:idxmean);
yque = qhist(1:idxmean);
objfunc = @(a) sum((yque-a.*yref).^2);
a0 = fminbnd(objfunc,0,1);
% a0 = 1; % standard area metric

% inferred on distribution
histdiff = qhist - a0*bhist;
histdiff(histdiff<0) = 0;
histdiff(1:idxmean) = 0;

a = sum(histdiff);

if a==0 || ~exist('datarange','var')
    onMean = nan;
    onMed = nan;
    onGeoMean = nan;
    return
end

% Compute on mean
datarange = reshape(datarange,[],1);
if numel(datarange) == numel(qhist)+1   % histcounts
    binCenters = datarange(1:end-1) + mean(diff(datarange));
else    % histc
    binCenters = datarange;
end

onPdf = histdiff./a;
onMean = log10(sum(onPdf.*10.^binCenters));
onGeoMean = sum(onPdf.*binCenters);

onCdf = [0; cumsum(onPdf)];
idx = find(onCdf>0.5,1);
onMed = interp1(onCdf((idx-1):idx), datarange((idx-1):idx), 0.5);