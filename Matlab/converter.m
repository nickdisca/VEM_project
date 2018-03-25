%%
%Nodes are ok
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

plot(Node(boundary,1),Node(boundary,2),'bo','markersize',2,'markerfacecolor','b');
%%
%Save in output file
file=fopen('100.dat','w');
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

%Save in .csv format
%??