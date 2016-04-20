function [metrics] = CalcYfpMetrics(strainData, makePlots)
% 20160114
if ~exist('makePlots','var')
    makePlots = [];
end

ncond = length(strainData.yfp);

ygrid = linspace(-3.5,.5,100);
xgrid = (1:ncond) - 0.5;

yfphist = nan(ncond,99);
onFrac = nan(ncond,1);
areaMetric = nan(ncond,1);
stdDev = nan(ncond,1);
thresh = nan(ncond,1);
eff = nan(ncond,1);

clear yfpnorm_ref;
for icond = 1:ncond
    if isempty(strainData.yfp{icond})
        onFrac(icond) = nan;
        areaMetric(icond) = nan;
        stdDev(icond) = nan;
        [thresh(icond), eff(icond)] = nan;
    else
        yfpnorm = log10(strainData.yfp{icond}./strainData.ssc{icond});
        
        % MATLAB 2015b-style histogram
        ntotal = sum(yfpnorm >= min(ygrid) & yfpnorm <= max(ygrid));
        yfphist(icond,:) = histcounts(yfpnorm, ygrid) ./ ntotal;
        
        onFrac(icond) = sum(yfpnorm > -2)./numel(yfpnorm);
        areaMetric(icond) = AreaMetric(yfphist(icond,:),yfphist(1,:));
        stdDev(icond) = std(yfpnorm);
        [thresh(icond), eff(icond)] = OtsuThresh(yfphist(icond,:), ygrid);
    end
end


metrics = table;

% scalar metrics
onFrac([1 end]) = nan;
idx = find(onFrac<0.5,1,'last');
if ~isempty(idx) && idx < ncond && ~any(isnan(onFrac(idx:(idx+1))))
    metrics.glu50_onFrac = 2-interp1(onFrac(idx:(idx+1)), xgrid(idx:(idx+1)), 0.5);
else
    metrics.glu50_onFrac = nan;
end

areaMetric([1 end]) = nan;
idx = find(areaMetric<0.5,1,'last');
if ~isempty(idx) && idx < ncond && ~any(isnan(areaMetric(idx:(idx+1))))
    metrics.glu50_areaMetric = 2-interp1(areaMetric(idx:(idx+1)), xgrid(idx:(idx+1)), 0.5);
else
    metrics.glu50_areaMetric = nan;
end

stdDev([1 end]) = nan;
metrics.sdMax = max(stdDev);
metrics.effMax = max(eff);

% array metrics
metrics.onFrac = onFrac';
metrics.stdDev = stdDev';
metrics.areaMetric = areaMetric';
metrics.thresh = thresh';
metrics.eff = eff';

%% plots
if ~isempty(makePlots)
    
    % histogram series heatmap
    
    % prevent pcolor from omitting last row+column
    zpadded = padarray(yfphist',[1,1],'post');
    xpadded = [xgrid, 2*xgrid(end) - xgrid(end-1)];
    
    pcolor(xpadded, ygrid, zpadded);
    
    colormap(flipud(bone));
    shading flat
    
    hold all
    box off
    ax = gca;
    ax.XTick = 2:3:11;
    ax.XTickLabel = {'0','-3','-6','-9'};
    
    xlim([min(xpadded)-0.01*range(xpadded), max(xpadded)]);
    ylim([min(ygrid)-0.01*range(ygrid), max(ygrid)]);
    
    % plot metrics
    xgrid = xgrid + 0.5;
    colors = lines;
    
    if any(strcmp(makePlots,'onfrac'))
        plot(xgrid, onFrac*range(ygrid) + min(ygrid),'o-','color',colors(1,:));
    end
    
    if any(strcmp(makePlots,'glu50_onfrac'))
        plot(2.5-metrics.glu50_onFrac([1 1]), min(ygrid)+[0 0.5].*range(ygrid),'.--','color',colors(1,:));
    end
    
    if any(strcmp(makePlots,'areametric'))
        plot(xgrid, areaMetric*range(ygrid) + min(ygrid),'o-','color',colors(2,:));
    end
    if any(strcmp(makePlots,'glu50_area'))
        plot(2.5-metrics.glu50_areaMetric([1 1]), min(ygrid)+[0 0.5].*range(ygrid),'.--','color',colors(2,:));
    end
    
    if any(strcmp(makePlots,'stddev'))
        plot(xgrid, stdDev+min(ygrid),'o-','color',colors(3,:));
    end
    
    if any(strcmp(makePlots,'otsuthresh'))
        plot(xgrid, thresh,'o-','color',colors(4,:));
    end
    
    if any(strcmp(makePlots,'otsueff'))
        plot(xgrid, eff*range(ygrid) + min(ygrid),'o-','color',colors(5,:));
    end
end
