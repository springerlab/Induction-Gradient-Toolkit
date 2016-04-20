function repMetrics = GatherReps(metricTable, idVar)
% GatherReps identifies and gathers together all the replicate measurements
% in the input table metricTable and saves them as a single table variable
% with multiple columns, 1 for each replicate. Missing data (i.e. strains
% for which fewer than the max number of replicates exist) are filled in as
% NaNs.
%
% Numeric variables in metricTable are removed. Text variables are kept,
% but only for the first replicate of a given strain.
%
% Last updated 20160121

if ~exist('idVar','var')
    idVar = 'name';
end

% remove non-scalar data
idxNonscalar = varfun(@(x) isnumeric(x) && size(x,2)>1, ...
    metricTable);
idxNonscalar = idxNonscalar{1,:};
metricTable(:,idxNonscalar) = [];

% get ready to combine scalar data
idxScalar = varfun(@(x) isnumeric(x) && size(x,2)==1, ...
    metricTable);
idxScalar = idxScalar{1,:};
scalarVars = metricTable.Properties.VariableNames(idxScalar);

% find replicates
n = 1;  % max number of replicates
istr = 1;
while istr <= height(metricTable)
    strNamesRemaining = metricTable.name((istr+1):end);
    idxRep = istr + find(strcmp(metricTable.name(istr), strNamesRemaining));
    
    if numel(idxRep) >= 1
        % increase number of replicates by padding with nans
        ndiff = numel(idxRep)+1 - n;
        if ndiff > 0
            for ivar = 1:length(scalarVars)
                m = metricTable.(scalarVars{ivar});
                metricTable.(scalarVars{ivar}) = [m, nan(size(m,1),ndiff)];
            end
            n = n + ndiff;
        end
        
        % populate additional replicate measurements
        for ivar = 1:length(scalarVars)
            for irep = 1:length(idxRep)
                metricTable.(scalarVars{ivar})(istr,irep+1) = ...
                    metricTable.(scalarVars{ivar})(idxRep(irep),1);
            end
        end
        
        metricTable(idxRep,:) = [];
    end
    istr = istr + 1;
end

repMetrics = metricTable;