function [a, onMean, onMed, onGeoMean] = AreaMetricHist(qhist, bhist, datarange)
% AreaMetricHist computes the fraction of a distribution that is above and
% outside a reference "basal" distribution.
%
% Updated 20160328
if any(isnan(bhist))
    a = nan;
    return
end

histdiff = qhist - bhist;
histdiff(histdiff<0) = 0;

bhist = reshape(bhist,[],1);

idxmean = floor(sum(bhist.*(1:length(bhist))')./sum(bhist));
histdiff(1:idxmean) = 0;

a = sum(histdiff);

if a==0
    onMean = nan;
    onMed = nan;
    onGeoMean = nan;
    return
end

onPdf = histdiff./a;
binCenters = datarange(1:end-1) + mean(diff(datarange));
onMean = log10(sum(onPdf.*10.^binCenters));
onGeoMean = sum(onPdf.*binCenters);

onCdf = [0 cumsum(onPdf)];
idx = find(onCdf>0.5,1);
onMed = interp1(onCdf((idx-1):idx), datarange((idx-1):idx), 0.5);