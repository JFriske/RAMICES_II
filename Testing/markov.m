
Nvec = 400;
basis = linspace(0,200,Nvec)';

m = 10*basis.*exp(-0.02 * basis) .* (sin(basis/15).^2+0.01);
% m = basis;
M = matrix(m);

conservationTest = M * m;



blue = finder(basis,m,10,1);
red = finder(basis,m,141,10);
green = finder(basis,m,60,2);
N = 100;
ts = 10.^linspace(-1,3,N);

clf;
 plot(basis,blue,'b');
 hold on;
 plot(basis,red,'r');
 plot(basis,green,'g');
 plot(basis,m,'k')
 hold off;
 ylim([1e-5,max(m)*1.2]);
 set(gca,'yscale','log');
 drawnow;

 delay = 0.05;
 h = figure(1);
frame = getframe(h); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
filename = "dissipation.gif";
 imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay); 
for i = 1:N
    mixer = expm((M - eye(Nvec))*ts(i));
   yR= mixer * red;
   yB = mixer * blue;
   yG = mixer * green;
%    hold on;
    plot(basis,m,'k');
    hold on;
   plot(basis,yR,'r');
   plot(basis,yB,'b');
   plot(basis,yG,'g');
%    [trapz(basis,m), trapz(basis,y), trapz(basis,dissipator)]
 set(gca,'yscale','log');
    ylim([1e-5,max(m)*1.2]);
    title("t = " + num2str(round(ts(i),3,'significant')) + "Gyr");
    drawnow
   hold off;
         frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay); 
end
% xlim([6,8])
function M = matrix(mvec)
    n = length(mvec);
    M = zeros(n);
    Mbar = max(mvec);
    k = 100;
    for i = 1:n
        
        for j = 1:n
            if i == j
               fac =  k/Mbar;
               summer = 0;
               if i < n
                   summer = summer+mvec(i+1);
               end
               if i > 1
                   summer = summer+mvec(i-1);
               end
               M(i,j) = 1.0-fac * summer;
               
               
            elseif abs(j-i) == 1
                fac =  k/Mbar;
                M(i,j) = fac * mvec(i);
            end
        end
    end

end

function diss = finder(basis,m,mu,sigma)
    [~,mID] = min(abs(basis -mu));
    diss = m(mID) * exp( -0.5*(basis - mu).^2/(sigma)^2);
end