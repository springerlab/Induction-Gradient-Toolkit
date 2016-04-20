function summaryMetrics = SummarizeReps(repMetrics)
% SummarizeReps computes the mean, standard deviation, and number of
% replicates for each data variable in the input table repMetrics. Missing
% data are assumed to be NaNs.
%
% The output table summaryMetrics will contain variables named with a
% metric and suffixed with '_mean', '_sd', or '_n'.
%
% 20160121

% only summarize numeric data
idxNum = varfun(@isnumeric, repMetrics);
idxNum = idxNum{1,:};
numericVars = repMetrics.Properties.VariableNames(idxNum);

for ivar = 1:length(numericVars)
    varName = numericVars{ivar};
    
    repMetrics.([varName '_mean']) = nanmean(repMetrics.(varName),2);
    repMetrics.([varName '_sd']) = nanstd(repMetrics.(varName),[],2);
    repMetrics.([varName '_n']) = sum(~isnan(repMetrics.(varName)),2);
    
    repMetrics.(varName) = [];
end

summaryMetrics = repMetrics;