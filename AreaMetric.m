function [a,onPop] = AreaMetric(query, ref)
% AreaMetric computes the fraction of a distribution that is above and
% outside another reference distribution.
%
% Updated 20160119
if isempty(ref) || fcsisempty(ref)
    a = nan;
    return
end

qhist = 
histdiff = qhist - rhist;
histdiff(histdiff<0) = 0;

rhist = reshape(rhist,1,[]);

idxmean = floor(sum(rhist.*(1:length(rhist)))./sum(rhist));
histdiff(1:idxmean) = 0;

a = sum(histdiff);