function strainData = LoadSegData(rootDir, expID)
% 20160421
if ischar(expID)
    expID = {expID};
end

strainData = [];
for iexp = 1:length(expID)
    fns = dir([rootDir '*' expID{iexp} '*']);
    if isempty(fns)
        continue;
    end
    
    datadir = [rootDir fns(1).name '/'];
    fns = dir([datadir 'segmented data ' expID{iexp} '.mat']);
    if isempty(fns)
        continue;
    end
    
    results = load([datadir fns(1).name]);
    
    strainData = [strainData; results.strainData];
end