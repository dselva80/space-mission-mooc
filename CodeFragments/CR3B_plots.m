%%
mu = 0.3;
x = linspace(-1.5,1.5,1000);
[x,y] = meshgrid(x,x);
r1 = sqrt((x+mu).^2+y.^2);
r2 = sqrt((x - (1-mu)).^2+y.^2);
U = -(x.^2+y.^2)/2 - ((1-mu)./r1 + mu./r2);
figure(1)
clf
contour(x,y,log(abs(U)),75)

tmp = -log(abs(U));
tmp(tmp < -0.7) = NaN;
figure(2)
clf
surface(x,y,tmp)
shading flat; view(72,66)

%%
f = 3;
%mov = avifile('hill_curve_anim.avi');
figure(f)
clf
feh = logspace(log10(abs(max(U(:)))),log10(abs(min(U(:)))/100),200);
for j=length(feh):-1:2
    figure(3)
    buh = zeros(size(U)); 
    buh(U > -feh(j)) = max(feh) - feh(j)+min(feh);
    imagesc(buh,[min(feh),max(feh)]);
    set(gca,'XTick',[],'YTick',[])
    pause(0.1); 
    %mov = addframe(mov,f);
end
%mov = close(mov); 
%% Earth-moon system
mu = 1/81;
x = linspace(-1.2,1.3,500);
[x,y] = meshgrid(x,x);
r1 = sqrt((x+mu).^2+y.^2);
r2 = sqrt((x - (1-mu)).^2+y.^2);
U = -(x.^2+y.^2)/2 - ((1-mu)./r1 + mu./r2);
figure(3)
clf
contour(x,y,log(abs(U)),[logspace(log10(0.4),log10(0.42),5),0.43:0.01:0.49,0.5:0.1:1,2,3,6,log(1/2*(3+mu-mu^2))])
text(0.79,0,'$L_1$','FontSize',16,'Interpreter','LaTex')
text(0.5-mu,sqrt(3)/2,'$L_4$','FontSize',16,'Interpreter','LaTex','HorizontalAlignment','center')
text(0.5-mu,-sqrt(3)/2,'$L_5$','FontSize',16,'Interpreter','LaTex','HorizontalAlignment','center')
text(1.1697,0,'$L_2$','FontSize',16,'Interpreter','LaTex','HorizontalAlignment','center')
text(-1,0,'$L_3$','FontSize',16,'Interpreter','LaTex','HorizontalAlignment','center')

%%
mu = 0.5;
x = linspace(-1.5,1.5,100);
f = x - (1 -mu)*(x+mu)./abs(x+mu).^3 - mu*(x - 1+mu)./abs(x - 1 + mu).^3;
figure(4)
clf
plot(x,f,'Linewidth',2)
hold on
plot([-1.5,1.5],[0,0],'k--')
ylim([-20,20])
%set(gca,'Xtick',[],'Ytick',[])

%%
D= linspace(0,100,1000);
r1 = D/2 - 1 + sqrt((D/2 - 1).^2 - (1+2*D).*(1 - D));
r2 = D/2 - 1 - sqrt((D/2 - 1).^2 - (1+2*D).*(1 - D));
r1(imag(r1) ~= 0) = NaN;
r2(imag(r2) ~= 0) = NaN;
figure(5)
clf
plot(D,sqrt(r1),D,-sqrt(r1),D,sqrt(r2),D,-sqrt(r2))

%%
msun = 132712440018;
mearth = 398600.4418;
mu = mearth/(msun+mearth);
tmp = @(x) x - (1 -mu)*(x+mu)./abs(x+mu).^3 - mu*(x - 1+mu)./abs(x - 1 + mu).^3;
l2 = fsolve(tmp,0.5);
%%
%     1 2 3  4  5  6
%Z = [x,y,z,xd,yd,zd]
dZ = @(t,Z) [Z(4:6);...
    2*Z(5)+Z(1)-(1-mu)*(Z(1)+mu)/((Z(1)+mu)^2 + Z(2)^2+Z(3)^2)^(3/2)-...
                    mu*(Z(1)-1+mu)/((Z(1)-(1-mu))^2 + Z(2)^2+Z(3)^2)^(3/2);...
    -2*Z(4)+Z(2)-(1-mu)*Z(2)/((Z(1)+mu)^2+Z(2)^2+Z(3)^2)^(3/2)-...
                    mu*Z(2)/((Z(1)-(1-mu))^2 + Z(2)^2+Z(3)^2)^(3/2);...
    -(1-mu)*Z(3)/((Z(1)+mu)^2+Z(2)^2+Z(3)^2)^(3/2)-mu*Z(3)/((Z(1)-(1-mu))^2+Z(2)^2+Z(3)^2)^(3/2)];

C0 = l2^2/2 - ((1-mu)/sqrt((l2+mu)^2) + mu/sqrt((l2-1+mu)^2));

figure(1)
clf
hold on
for j = 1:10000
    j
    x0 = l2+randn(1)*0.1;
    y0 = randn(1)*0.1;
    K = 2*(C0 + (x0^2+y0^2)/2 + ((1-mu)/sqrt((x0+mu)^2+y0^2) + mu/sqrt((x0-1+mu)^2+y0^2)));
    if K <= 0; 
        disp('bah')
           continue; 
    end
    a = rand(1);
    xd0 = sqrt(a*K);
    yd0 = sqrt((1-a)*K);
    [t,Z] = ode113(dZ,[0,10],[x0;y0;0;xd0;yd0;0]);
    plot(Z(:,1),Z(:,2),'.')
end

%%

[t,Z] = ode113(dZ,[0,100],[l2+eps;ones(5,1)*eps*10],odeset('AbsTol',1e-12));
    plot(Z(:,1),Z(:,2),'.')
%%
Gamma=0.5;
deq=@(t,x) [x(2);x(1)-0.3*x(2)-(x(1))^3+Gamma*cos(1.25*t)];
options=odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,xx]=ode45(deq,0:(2/1.25)*pi:(4000/1.25)*pi,[1,0]);
plot(xx(:,1),xx(:,2),'.','MarkerSize',1)
fsize=15;
axis([-2 2 -2 2])
xlabel('x','FontSize',fsize)
ylabel('y','FontSize',fsize)
title('Poincare Section of the Duffing System')