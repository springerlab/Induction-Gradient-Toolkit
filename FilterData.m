function datafilt = FilterData(data,fthresh,nContours)

data = fcsselect(data, data.fsc<9999 & data.ssc < 9999);

%             x = data.fsc;
%             y = data.ssc;
%             xi = linspace(0,1e4,50);
%             yi = linspace(0,1e4,50);
x = log10(data.fsc);
y = log10(data.ssc);
xi = linspace(2.5,4,50);
yi = linspace(2.5,4,50);
%             d=fcsdensity(x,y,'xi',xi,'yi',yi);
%             hold all
%             xlim([2.5 4]);
%             ylim([2.5 4]);
%             caxis([0 0.01]);

d = hist3([x y],'edges',{xi' yi'})./numel(x);

d = imgaussfilt(d,1.5);
d(d<=5e-5) = nan;
%             d = log10(d);

%             fthresh = 0.5;
%             nContours = 20;

allThresh = logspace(log10(max(d(:))),log10(min(d(:)))+1,nContours);
for ith = 1:length(allThresh)
    bw = d>allThresh(ith);
    f = sum(d(d>allThresh(ith)));
    if f>fthresh
        [B,L] = bwboundaries(bw);
        idx = zeros(length(x),1);
        
        for i=1:length(B)
            pixels = B{i};
            xb = xi(pixels(:,1));
            yb = yi(pixels(:,2));
            
            idx = idx | inpolygon(x,y,xb,yb);
            datafilt = fcsselect(data,idx);
            
            %                         plot(xb,yb,'color','r','LineWidth',1);
        end
        break;
    end
end