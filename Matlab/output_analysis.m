%%
clear all; close all; clc;
load('../output.dat');

n_elem=64;
h=1/n_elem;
output=output(1:(n_elem+1)^2,:);
output=sortrows(output,[1 2 3]);

y=output(output(:,1)==0,2);
x=output(output(:,2)==0,1);
[X,Y]=meshgrid(x,y);

sol=output(:,3);
sol=reshape(sol,[n_elem+1 n_elem+1]);
contourf(X, Y, sol);
hold on;
colorbar;
colormap jet;
titletxt = 'Contour plot of VEM solution ($k=1$)';
title(titletxt,'fontsize',12,'interpreter','latex');
xlabel({'$$x$$'},'fontsize',12,'interpreter','latex')
ylabel({'$$y$$'},'fontsize',12,'interpreter','latex')

%%
clear all; close all; clc;
load('../output.dat');

x = output(:,1);
y = output(:,2);
v = output(:,3);

n_elem=64; k=5;
n_pts=n_elem*k;
xx=linspace(0,1,n_pts); yy=linspace(0,1,n_pts);
[X,Y] = meshgrid(xx, yy);
sol = griddata(x,y,v,X,Y);

contourf(X, Y, sol);
hold on;
colorbar;
colormap jet;
titletxt = 'Contour plot of VEM solution ($k=2$)';
title(titletxt,'fontsize',12,'interpreter','latex');
xlabel({'$$x$$'},'fontsize',12,'interpreter','latex')
ylabel({'$$y$$'},'fontsize',12,'interpreter','latex')
%% OLD VERSION
%esempio
%X=[1.7 2.3 3.5 4.0 -5 6 7 8 9]';
%Y=[0 2.3 3 4 5 6 7 10 9]';
%U=[0.5 0.3 0.2 0.1 0.9 0.5 0.7 0.6 0.1]';
%plot3(X,Y,U,'bo','markerfacecolor','b','markersize',4);
%grid on;
%[xq,yq] = meshgrid(1.0:0.1:9.0, 1.0:0.1:9.0);
%vq = griddata(X,Y,U,xq,yq);
%mesh(xq,yq,vq);
%%
%aux=load('sinsin.dat');
%X=aux(:,1); Y=aux(:,2); U=aux(:,3);
%figure
%plot3(X,Y,U,'bo','markerfacecolor','b','markersize',4);
%grid on;
%hold on;
%fff=@(x,y) sin(x).*sin(y);
%fff=@(x,y) x.^2+y.^2;
%F=fff(X,Y);
%plot3(X,Y,F,'rs','markerfacecolor','r','markersize',4);
%xlabel('x'); ylabel('y'); zlabel('u');
%legend('approximated','exact');
%[value,index]=max(abs(F-U));
%point=[X(index,1); Y(index,1)]
%disp(['Infinity norm: ', num2str(max(abs(F-U)))]);