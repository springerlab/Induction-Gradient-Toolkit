function [metrics] = PlotYfpProfile(strainData)
% 20160114

ncond = length(strainData.yfp);

ygrid = linspace(-3.5,.5,100);
xgrid = (1:ncond) - 0.5;

yfphist = nan(ncond,99);

for icond = 1:ncond
    yfpnorm = log10(strainData.yfp{icond}./strainData.ssc{icond});
    
    % MATLAB 2015b-style histogram
    ntotal = sum(yfpnorm >= min(ygrid) & yfpnorm <= max(ygrid));
    yfphist(icond,:) = histcounts(yfpnorm, ygrid) ./ ntotal;
end

% pad data so that pcolor doesn't cut off last row+column of data
zpadded = padarray(yfphist',[1,1],'post');
xpadded = [xgrid, 2*xgrid(end) - xgrid(end-1)];

pcolor(xpadded, ygrid, zpadded);

colormap(flipud(bone));
shading flat

hold all
box off

xlim([min(xpadded)-0.01*range(xpadded), max(xpadded)]);
ylim([min(ygrid)-0.01*range(ygrid), max(ygrid)]);