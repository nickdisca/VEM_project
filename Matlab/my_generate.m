%%Generatore di mesh uniforme e salva su file

%%
%Genera la mesh
clear all; close all;
nelx=8; nely=8;
dx = 1/nelx; dy = 1/nely;
[X,Y] = meshgrid(dx/2:dx:1,dy/2:dy:1);
P = [X(:) Y(:)];
[Node,Element,Supp,Load,P] = PolyMesher(@MbbDomain,1500,0,P);
%se cambio la mesh devo cambiare anche MbbDomain

%%
hold on
axis([0 1 0 1])
axis tight
titletxt = 'Uniform grid';
title(titletxt,'fontsize',12,'interpreter','latex');
xlabel({'$$x$$'},'fontsize',12,'interpreter','latex')
ylabel({'$$y$$'},'fontsize',12,'interpreter','latex')

%%
%Nodes
figure
plot(Node(:,1),Node(:,2),'bo','markersize',2,'markerfacecolor','b');
%Elements
%Convert cell into matrix
aux=cell2mat(Element);
for i=1:length(aux)
    for j=1:4
        A(i,j)=aux(i,j);
    end
end
clear('aux');
%Boundary
boundary=[]; j=1;
for i=1:size(Node,1)
    if (1-Node(i,1)<=1e-10 || 1-Node(i,2)<=1e-10 || Node(i,1)<=1e-10 || Node(i,2)<=1e-10)
        %boundary(j,1)=Node(i,1); boundary(j,2)=Node(i,2);
        boundary(j,1)=i;
        j=j+1;
    end
end
figure
plot(Node(boundary,1),Node(boundary,2),'bo','markersize',2,'markerfacecolor','b');
%%
%Save in output file
file=fopen('prova.dat','w');
for i=1:size(Node,1)
    fprintf(file,'%f %f \n',Node(i,:));
end
for i=1:size(A,1)
    fprintf(file,'%d %d %d %d \n',A(i,:));
end
for i=1:size(boundary,1)
    fprintf(file,'%d \n',boundary(i));
end
fclose(file);