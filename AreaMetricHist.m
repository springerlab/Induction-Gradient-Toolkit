function [a, onMean, onMed, onGeoMean] = AreaMetricHist(qhist, bhist, datarange)
% AreaMetricHist used to do something different from AreaMetric, but now
% does the same thing, so it points to AreaMetric
%
% Updated 20160328

if exist('datarange','var')
    [a,onMean,onMed,onGeoMean] = AreaMetric(qhist,bhist,datarange);
else
    [a,onMean,onMed,onGeoMean] = AreaMetric(qhist,bhist);
end
