%%Generatore di mesh NON uniforme (Voronoi) e salva su file
clear all
close all
%%
%Genera la mesh
clear all
nelx=8; nely=8;
dx = 1/nelx; dy = 1/nely;
[X,Y] = meshgrid(dx/2:dx:1,dy/2:dy:1);
p = [X(:) Y(:)];
d = dRectangle(p,0,1,0,1); %quadrato unitario
[Node,Element,Supp,Load,P] = PolyMesher(@MbbDomain,128,100);

%%
hold on
axis([0 1 0 1])
axis tight
titletxt = 'Non uniform grid';
title(titletxt,'fontsize',12,'interpreter','latex');
xlabel({'$$x$$'},'fontsize',12,'interpreter','latex')
ylabel({'$$y$$'},'fontsize',12,'interpreter','latex')
%%
%Nodes
figure
file=fopen('128.dat','w');
plot(Node(:,1),Node(:,2),'bo','markersize',2,'markerfacecolor','b');
for i=1:size(Node,1)
    fprintf(file,'%f %f \n',Node(i,:));
end
%%
%Elements
%Convert cell into matrix
for i=1:length(Element)
    aux=Element(i,1);
    myElem=cell2mat(aux);
    myElem=mat2str(myElem);
    myElem=myElem(2:end-1);
    myElem=[myElem, ' ','\n'];
    fprintf(file,myElem);
end

%%
%Boundary
boundary=[]; j=1;
for i=1:size(Node,1)
    if (1-Node(i,1)<=1e-6 || 1-Node(i,2)<=1e-6 || Node(i,1)<=1e-6 || Node(i,2)<=1e-6)
        %boundary(j,1)=Node(i,1); boundary(j,2)=Node(i,2);
        boundary(j,1)=i;
        j=j+1;
    end
end
figure
plot(Node(boundary,1),Node(boundary,2),'bo','markersize',2,'markerfacecolor','b');
for i=1:size(boundary,1)
    fprintf(file,'%d \n',boundary(i));
end
fclose(file);

%% compute diameter
distance=@(p1,p2) sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
H=0;
for i=1:length(Element)
    aux=Element(i,1);
    myElem=cell2mat(aux);
    diam=0.0;
    for j=1:length(myElem)
        for k=1:length(myElem)
            p1=Node(myElem(j),:); p2=Node(myElem(k),:);
        diam=max(diam,distance(p1,p2));
        end
    end
    disp(diam);
    H=max(H,diam);
end
stri=['Global diameter is h=', num2str(H)];
disp(stri);