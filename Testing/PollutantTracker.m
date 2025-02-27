set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',23);

files = "../Output/Pollutant/" + ["Active","Quiescent","Active_NoInflow","Quiescent_NoInflow","Active_NoDispersion","Quiescent_NoDispersion"] + "/Enrichment_Absolute_ColdGas.dat";
% files = files(1)

names = [ "Active","Quiescent","Active, No Inflow","Quiescent, No Inflow","Active, No Diffusion","Quiescent, No Diffusion"];
times = [0.1,3,10];

figure(1);
clf;
T = tiledlayout(1,3,'Padding','None','TileSpacing','None');
for file = files
    file
    track(file,times)
end

fs = 24;
nexttile(1);
legend(names,"FontSize",22);
xlabel(T,"Radius (kpc)","Interpreter","latex","FontSize",fs)
ylabel(T,"Pollutant Density (Arbitrary units)","Interpreter","latex","FontSize",fs)
function track(fileName,times)

    cs = colororder;
    col = cs(1,:);
    if (contains(fileName,"Inflow"))
        col = cs(2,:);
    end
    if (contains(fileName,"Dispersion"));
        col = cs(3,:);
    end
    style = "-";
    if (contains(fileName,"Quiescent"))
        style = "--";
    end
    


    floor = 1e-10;
    opts = detectImportOptions(fileName);

    opts.VariableTypes(:) = {'double'};

    f = readtable(fileName,opts);

    simTimes = unique(f.Time);
    tSim = [];
	maxVal = 0;
    for i = 1:length(times)
        t = times(i);
        [~,I] = min(abs(simTimes - t));
        tSim(end+1) = I;
%         simTimes(I)
%         [t,I,simTimes(I)]
        nexttile(i);
        cut = (f.Time == simTimes(I));    
        focus = f(cut,:);
%         focus(1:3,:)
        
        hold on;
        v = focus.Total_Eu ./ focus.Total_H;
%         m = mean(v);
%         v = sqrt( (v - m).^2);
        
%         totalEu = sum(v)
		if i == 1
			maxVal = max(v);
            maxH = max(focus.Total_H);
        end
        v = v/maxVal;
%         v(v<floor) = floor * 0.8;
        plot(focus.RingRadius,v,style,'Color',col,'LineWidth',2);
%         plot(focus.RingRadius,focus.Total_H/maxH,':','Color','k','HandleVisibility','Off');
        hold off;
        set(gca,'yscale','log')
        xlim([0,20]);
        ylim([floor,1]);

        grid on;
        title("$t = " + num2str(t) + "$ Gyr");
        
        
    end
    

end