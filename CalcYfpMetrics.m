function [metrics, h] = CalcYfpMetrics(strainData, makePlots)
% 20160114
if ~exist('makePlots','var')
    makePlots = [];
end

ncond = length(strainData);

ygrid = linspace(-3.5,.5,100);
xgrid = (1:ncond) - 0.5;

yfphist = nan(ncond,99);
onFrac = nan(ncond,1);
areaMetric = nan(ncond,1);
onMean = nan(ncond,1);
onMed = nan(ncond,1);
onGeoMean = nan(ncond,1);
meanMetric = nan(ncond,1);
onMed = nan(ncond,1);
stdDev = nan(ncond,1);
thresh = nan(ncond,1);
eff = nan(ncond,1);

for icond = 1:ncond
    if isempty(strainData{icond})
        onFrac(icond) = nan;
        areaMetric(icond) = nan;
        stdDev(icond) = nan;
        thresh(icond) = nan;
        eff(icond) = nan;
    else
        yfpnorm = log10(strainData{icond}.yfp./strainData{icond}.ssc);
        
        % MATLAB 2015b-style histogram
        ntotal = sum(yfpnorm >= min(ygrid) & yfpnorm <= max(ygrid));
        yfphist(icond,:) = histcounts(yfpnorm, ygrid) ./ ntotal;
        
        onFrac(icond) = sum(yfpnorm > -2.5)./numel(yfpnorm);
        [areaMetric(icond), onMean(icond), onMed(icond), onGeoMean(icond)] = ...
            AreaMetricHist(yfphist(icond,:),yfphist(1,:),ygrid);
        meanMetric(icond) = log10(mean(10.^yfpnorm));
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

thresh_lo = 0.2;
thresh_hi = 0.8;
idx = find(areaMetric<thresh_lo,1,'last');
if ~isempty(idx) && idx < ncond && ~any(isnan(areaMetric(idx:(idx+1))))
    glu_lo = 2-interp1(areaMetric(idx:(idx+1)), xgrid(idx:(idx+1)), thresh_lo);
else
    glu_lo = nan;
end
idx = find(areaMetric<thresh_hi,1,'last');
if ~isempty(idx) && idx < ncond && ~any(isnan(areaMetric(idx:(idx+1))))
    glu_hi = 2-interp1(areaMetric(idx:(idx+1)), xgrid(idx:(idx+1)), thresh_hi);
else
    glu_hi = nan;
end
metrics.areaMetricRange = glu_lo - glu_hi;

onMean_thresh = onMean(end) - 1;
onMean([1 end]) = nan;
idx = find(onMean<onMean_thresh,1,'last');
if ~isempty(idx) && idx < ncond && ~any(isnan(onMean(idx:(idx+1))))
    metrics.onMean10 = ...
        2-interp1(onMean(idx:(idx+1)), xgrid(idx:(idx+1)), onMean_thresh);
else
    metrics.onMean10 = nan;
end

mean_thresh = meanMetric(end) - 1;
meanMetric([1 end]) = nan;
idx = find(meanMetric<mean_thresh,1,'last');
if ~isempty(idx) && idx < ncond && ~any(isnan(meanMetric(idx:(idx+1))))
    metrics.meanMetric10 = ...
        2-interp1(meanMetric(idx:(idx+1)), xgrid(idx:(idx+1)), mean_thresh);
else
    metrics.meanMetric10 = nan;
end

stdDev([1 end]) = nan;
metrics.sdMax = max(stdDev);
metrics.effMax = max(eff);
metrics.sdInt = nansum(stdDev);

% array metrics
metrics.onFrac = onFrac';
metrics.stdDev = stdDev';
metrics.areaMetric = areaMetric';
metrics.thresh = thresh';
metrics.eff = eff';

%% plots
h = [];
iplot = 1;
if ~isempty(makePlots)
    % prevent pcolor from omitting last row+column
    zpadded = padarray(yfphist',[1,1],'post');
    xpadded = [xgrid, 2*xgrid(end) - xgrid(end-1)];
    
    hold all
    
    if any(strcmp(makePlots,'yfpprofile'))
        % histogram series heatmap
        pcolor(xpadded, ygrid, zpadded);
        
        colormap(flipud(bone));
        shading flat
        caxis([0 0.15]);
    end
    box off
    ax = gca;
    ax.XTick = 2:3:11;
    ax.XTickLabel = {'0','-3','-6','-9'};
    ax.YTick = -3:0;
    ax.YTickLabelMode = 'auto';
    
    xlim([min(xpadded)-0.01*range(xpadded), max(xpadded)]);
    ylim([min(ygrid)-0.01*range(ygrid), max(ygrid)]);
    
    % plot metrics
    xgrid = xgrid + 0.5;
    colors = lines(9);
    
    if any(strcmp(makePlots,'onfrac'))
        h(iplot) = plot(xgrid, onFrac*range(ygrid) + min(ygrid),'-','color',colors(1,:));
        iplot = iplot + 1;
    end
    
    if any(strcmp(makePlots,'glu50_onfrac'))
        h(iplot) = plot(2.5-metrics.glu50_onFrac([1 1]), min(ygrid)+[0 0.5].*range(ygrid),'.--','color',colors(1,:));
        iplot = iplot + 1;
    end
    
    if any(strcmp(makePlots,'areametric'))
        h(iplot) = plot(xgrid, areaMetric*range(ygrid) + min(ygrid),'-','color',colors(2,:));
        iplot = iplot + 1;
    end
    if any(strcmp(makePlots,'glu50_area'))
        h(iplot) = plot(2.5-metrics.glu50_areaMetric,...
            min(ygrid)+0.5.*range(ygrid),'o','color',colors(2,:));
        h(iplot) = plot(2.5-metrics.glu50_areaMetric([1 1]),...
            min(ygrid)+[0 0.5].*range(ygrid),'--','color',colors(2,:));
        iplot = iplot + 1;
    end
    if any(strcmp(makePlots,'areametric_range'))
        h(iplot) = plot(2.5-[glu_lo glu_hi], ...
            min(ygrid)+[thresh_lo thresh_hi].*range(ygrid),'o','color',colors(2,:));
        iplot = iplot + 1;
    end
    if any(strcmp(makePlots,'stddev'))
        h(iplot) = plot(xgrid, stdDev+min(ygrid),'-','color',colors(3,:));
        iplot = iplot + 1;
    end
    
    if any(strcmp(makePlots,'otsuthresh'))
        h(iplot) = plot(xgrid, thresh,'-','color',colors(4,:));
        iplot = iplot + 1;
    end
    
    if any(strcmp(makePlots,'otsueff'))
        h(iplot) = plot(xgrid, eff*range(ygrid) + min(ygrid),'-','color',colors(5,:));
        iplot = iplot + 1;
    end
    if any(strcmp(makePlots,'meanmetric'))
        h(iplot) = plot(xgrid, meanMetric,'-','color',colors(6,:));
        iplot = iplot + 1;
    end
    if any(strcmp(makePlots,'onmean'))
        y = onMean;
        idx = find(areaMetric>=0.1,1);
        h(iplot) = plot(xgrid(1:idx), y(1:idx),'--','color',colors(1,:));
        h(iplot) = plot(xgrid(idx:end), y(idx:end),'-','color',colors(1,:));
        iplot = iplot + 1;
    end
    
    if any(strcmp(makePlots,'onmean10'))
        h(iplot) = plot(2.5-metrics.onMean10, onMean_thresh,...
            'o','color',colors(1,:));
        h(iplot) = plot(2.5-metrics.onMean10([1 1]),...
            [min(ygrid) onMean_thresh],'--','color',colors(1,:));
        iplot = iplot + 1;
    end
    
    if any(strcmp(makePlots,'onmed'))
        y = onMed;
        y(areaMetric<0.1) = nan;
        h(iplot) = plot(xgrid, y,'-','color',colors(8,:));
        iplot = iplot + 1;
    end
    if any(strcmp(makePlots,'ongeomean'))
        y = onGeoMean;
        y(areaMetric<0.1) = nan;
        h(iplot) = plot(xgrid, y,'-','color',colors(9,:));
        iplot = iplot + 1;
    end
end
