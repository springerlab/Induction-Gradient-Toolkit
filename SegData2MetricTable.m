function metricTable = SegData2MetricTable(segData, expID)
% SegData2MetricTable loops through segmented FC data for each strain of an
% experiment and calculates metrics for that strain, saving the results in
% a table.
%
% Updated 20160119
metricTable = table;

for istr = 1:length(segData)
    qdata = segData(istr).query;
    
    m = table;
    m.name = {segData(istr).name};
    m.info = {segData(istr).info};
    m.expID = {expID};
    m = [m, CalcYfpMetrics(qdata)]; % metrics calculated by this function
    
    metricTable = [metricTable; m];
end
