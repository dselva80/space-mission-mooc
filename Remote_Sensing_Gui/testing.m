% 

% plot world horizon on legend axis
x_world = -1:.1:3;
y_world = -.03*(x_world-1).^2+.3;
plot(x_world,y_world,'LineWidth',2)

hold on;
% plot FOR lines on legend axis
x_for = 0:.1:2.5;
y_for = -.9*abs(x_for-1)+1;
plot(x_for,y_for,'--r')
x_forlabel = .81:.01:1.19;
y_forlabel = (x_forlabel-1).^2+.8;
plot(x_forlabel,y_forlabel,'r')
text(.5,.83,'FOR','Color','r','FontSize',13)
plot([0 2.5], [.1 .1],'r')

% plot FOV lines on legend axis
x_fov = .953:.01:1.053;
y_fov = -19*abs(x_fov-1)+1;
absy = plot(x_fov,y_fov,'k');
x_fovlabel = .92:.01:1.08;
y_fovlabel = (x_fovlabel-1).^2+.48;
plot(x_fovlabel,y_fovlabel,'k');
text(1.18,.48,'FOV','Color','k','FontSize',13)
x1 = .55;
twist = .16;
width =.8;
y1 = 0.01;
height = .18;
x2 = x1+width;
y2 = y1+height;
plot([x1 x1+twist x2 x2-twist x1],[y1 y2 y2 y1 y1],'k');

% plot pixels IFOV
po =  .07;
plot([x1+po x1+twist-po x2-po x2+po-twist x1+po],[y1+po y2-po y2-po y1+po y1+po],'k');
rotate(absy,[0 1 0],40)
pwidth = .08;
px = [x1+po+pwidth x1+pwidth+po+twist/4];
py = [y1+po y2-po];
plot(px, py,'k')
px2 = [x1+po+pwidth*2 x1+pwidth*2+po+twist/4];
plot(px2, py,'k');
px3 = [x1+po+pwidth*3 x1+pwidth*3+po+twist/4];
plot(px3, py,'k');
px4 = [x1+po+pwidth*4 x1+pwidth*4+po+twist/4];
plot(px4, py,'k');
px5 = [x1+po+pwidth*5 x1+pwidth*5+po+twist/4];
plot(px5, py,'k');



axis('off')
plot(1,1,'sk','MarkerFaceColor','k','MarkerSize',20)
