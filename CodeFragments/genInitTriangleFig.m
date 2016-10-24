function genInitTriangleFig(r1v,r2v,f)

figure(f)
hold on
plot(0,0,'.')
set(gca,'FontName','Times','FontSize',14)
set(gca,'DefaulttextFontName','Times','DefaulttextFontSize',...
    12,'DefaulttextInterpreter','LaTeX','Box','on')
hold on
plot(r1v(1),r1v(2),'.')
plot([0,r1v(1)],[0,r1v(2)],'k')
plot(r2v(1),r2v(2),'.')
plot([0,r2v(1)],[0,r2v(2)],'k')
plot([r1v(1),r2v(1)],[r2v(2),r2v(2)],'k')
axis equal
hold off
set(gca,'FontSize',18,'FontName','Times')

%annotate
shim = max(diff(axis))*0.02/3;
text(0, -shim,'F','HorizontalAlignment','center','VerticalAlignment','top','FontSize',18,'FontName','Times')
text(r1v(1)-shim,r1v(2),'$$P_1$$','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',18,'FontName','Times')
text(r2v(1)+shim,r2v(2),'$$P_2$$','VerticalAlignment','bottom','FontSize',18,'FontName','Times')