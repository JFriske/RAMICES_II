set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',20);


% T=tiledlayout(4,2);
files = "../Output/Pollution/" + ["Active"] + "/StellarCatalogue.dat";

plotter(files,1)
function plotter(files,i)
    q = figure(i);
    clf;
    T = tiledlayout('flow');
    for f = files
        g = readtable(f,"ReadVariableNames",true);
        disp("Loaded")
        cut = g.FeH < -10;
        g(cut,:) = [];
        disp("Cut")

%         nexttile;
%         histogram2(g.Radius,g.BirthRadius,15,'FaceColor','flat',"ShowEmptyBins",true);
%         set(gca,'zscale','log');
          set(gca,'colorscale','log')
          xlabel("Final Annulus");
          ylabel("Birth Annulus");
          colorbar;
        view(2) 
    %     
%         nexttile;
%         histogram2(g.MeasuredAge, g.Radius - g.BirthRadius,20,'FaceColor','flat',"ShowEmptyBins",true)
%         set(gca,'zscale','log');
        set(gca,'colorscale','log')
        xlabel("Age (Gyr)");
        ylabel("Outward Migration (kpc)");
        view(2) 
        
        
%         nexttile;
        delta = g.MgH - g.FeH;
        cutter = ~((delta > 0.6) | (delta < -0.5) | (g.FeH < -5));
    %     histogram2(g.FeH(cutter),delta(cutter),[20,30])

%         spatialDelta = g.BirthRadius;
        n = height(g);
        
%         scaling = 
        
        xDelta = normrnd(0,1,n,1) .*  0.0;
        yDelta = normrnd(0,1,n,1) .* 0.0;
%         scatter(g.FeH,delta,3,g.BirthRadius,'filled')
        
        disp("Plotted 1")
        Ncols = 100000;
        zeroed = hot(Ncols);
        zeroed(1,:) = [0,0,0];
        nexttile;

        
        colorbar
        x = g.FeH + xDelta;
        y = delta + yDelta;
        [N,X,Y] = histcounts2(x,y,300);
        N = N';


        colormap(zeroed)
        image([min(X),max(X)],[min(Y),max(Y)],N,'CDataMapping','scaled')
      
        set(gca,'YDir','normal')
        set(gca,'ColorScale','log')
        disp("Plotted 2")
        xlim([-2.4,0.6]);
%         ylim([-0.2,0.4]);
        colorbar
        grid on;

%         nexttile;
%         histogram(y,30)
        
%         set(gca,'yscale','log')
          clear x y N xDelta yDelta; 
    %     scatter(g.FeH,delta,1)

    %    
    %    
        nexttile;
        color = g.BMag - g.VMag;
        scatter(color,g.VMag,3,g.MeasuredAge,'filled');
        histogram2(color(cutter),g.VMag(cutter),'FaceColor','flat');
        set(gca,'ydir','reverse');
        set(gca,'ColorScale','log')
        colorbar
        view(2);


%         nexttile;
%         histogram(g.Mass)
%         disp("Plotted 3")
%         set(gca,'yscale','log');

        thickSampler = (g.MeasuredAge < 3);
        z0 = 0.02;
        kappa = 0.3;
        pow = 0.66;
        
        thickScale = mean( z0 + kappa * g.MeasuredAge(thickSampler).^pow)
        thinScale = mean(z0 + kappa * g.MeasuredAge(~thickSampler).^pow)
        
    

        
    end
end
    