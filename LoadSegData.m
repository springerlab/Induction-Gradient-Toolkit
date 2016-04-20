function results = LoadSegData(rootDir, expID)
% 20160412
fns = dir([rootDir '*' expID '*']);
if isempty(fns)
    results = [];
    return
end

rootDir = [rootDir fns(1).name '/'];
fns = dir([rootDir 'segmented data ' expID '.mat']);
if isempty(fns)
    results = [];
    return
end

results = load([rootDir fns(1).name]);