function [xyz,mesh,nn,tri,bnd,map] = importmesh2(fileToRead1)
%READ FILE
fprintf('Read file %s\n',fileToRead1);
DELIMITER = ' ';
HEADERLINES = 4;
nodes = importdata(fileToRead1, DELIMITER, HEADERLINES);
nn=nodes.data(1);
HEADERLINES = nn+7;
elements = importdata(fileToRead1, DELIMITER, HEADERLINES);
%END READ FILE

%PROCESS MESH INFORMATION
%nn  - number of nodes in mesh
%xyz - (nn,2) array of point data
%ien - elements from *.msh file
%tri - connectivity mesh for visualization
xyz=nodes.data(2:length(nodes.data));
xyz=reshape(xyz,4,[])';
xyz=xyz(:,2:3);
ien=elements.data(2:length(elements.data));
tri=zeros(elements.data(1),3);

i=2;
bnd=[];
e=0;
%This "while" loop determines the number of triangular elements and where in
%the "gmsh" mesh the triangular elements begin.
while i<length(ien)
    if ien(i)==15 %15 corresponds to point elements
        i=i+6;
        tri(length(tri),:)=[];
    elseif ien(i)==1 %1 corresponds to line elements
        bnd=[bnd ien(i+4)];
        bnd=[bnd ien(i+5)];
        i=i+7;
        tri(length(tri),:)=[];
    elseif ien(i)==2 %corresponds to the tri's in our mesh
        e=e+1;
        tri(e,:)=ien(i+4:i+6);
        i=i+8;
    else disp('error')
    end
end
bnd=sort(unique(bnd));

xmin=min(xyz(:,1));
xmax=max(xyz(:,1));
ymin=min(xyz(:,2));
ymax=max(xyz(:,1));
fprintf('Mesh bounds: x=[%g, %g], y=[%g, %g] \n',xmin, xmax,ymin,ymax)
%END PROCESS MESH INFORMATION

%CREATE NEIGHBOR INFORMATION
mesh = struct('coeff',{},'bnd',{},'nb1',{},'nbr',{},'dx',{},'dy',{});
for n=1:nn
    mesh(n).nbr = [];%neighbors
    mesh(n).dx  = [];
    mesh(n).dy  = [];
end

for i=1:length(tri)
    for j=1:3
        for k=1:2
            node1=tri(i,j);
            node2=tri(i,mod(j+k-1,3)+1);
            if isempty(find(mesh(node1).nb1==node2))
                mesh(node1).nb1=[mesh(node1).nb1 node2];
            end
        end
    end
end

count=0;
skip=[]; %some *.msh files have nodes with no connectivity.
for n=1:nn
    if isempty(mesh(n).nb1)
        skip=[skip n];
        count=count+1;
    end
end
fprintf('Nodes with no connectivity: %g \n\n',count)

minn=1e5;
maxn=0;
lt9=0;
for n=1:nn
    nbo=mesh(n).nb1;
    for t=1:length(nbo)
        nb2=mesh(nbo(t)).nb1;
        for i=1:length(nb2)
            if isempty(find(mesh(n).nbr==nb2(i),1)) && nb2(i)~=n
                mesh(n).nbr=[mesh(n).nbr nb2(i)];
            end
        end
    end
    if length(mesh(n).nbr)<minn && ismember(n,skip)==0
        minn=length(mesh(n).nbr);
    end
    if length(mesh(n).nbr)<9
        lt9=lt9+1;
    end
    if length(mesh(n).nbr)>maxn
        maxn=length(mesh(n).nbr);
    end
end
mesh = rmfield(mesh,'nb1'); 
fprintf('Max neighbors: %g \n',maxn)
fprintf('Min neighbors: %g \n',minn)
fprintf('Nodes with 1st order spacial accuracy: %g \n\n',lt9)
for n=1:nn
    for i=1:length(mesh(n).nbr)
        n2=mesh(n).nbr(i);
        mesh(n).dx=[mesh(n).dx xyz(n2,1)-xyz(n,1)];
        mesh(n).dy=[mesh(n).dy xyz(n2,2)-xyz(n,2)];
    end
end
%END CREATE NEIGHBOR INFORMATION

%CREATE BOUNDARY INFORMATION
for n=1:nn
    mesh(n).bnd=-1;
end

tol=1e-12;
%bnd numbers 1-4 correspond to outer walls, 5 corresponds to inner circle
for i=1:length(bnd)
    e=bnd(i);
    if abs(xyz(e,1)-xmin)<tol
        mesh(e).bnd=1; continue;
    elseif abs(xyz(e,2)-ymin)<tol
        mesh(e).bnd=2; continue;
    elseif abs(xyz(e,1)-xmax)<tol
        mesh(e).bnd=3; continue;
    elseif abs(xyz(e,2)-ymax)<tol
        mesh(e).bnd=4; continue;
    end
    mesh(e).bnd=5;
end
for i=1:length(skip)
    mesh(skip(i)).bnd=99; %bnd=99 will be skipped in "heat.m"
    bnd=[bnd skip(i)];
end

bnd=sort(unique(bnd));
map=zeros(1,(nn-length(bnd)));
inc=1;
for i=1:nn
    if ismember(i,bnd)
        continue
    end
    map(inc)= i;
    inc=inc+1;
end
fprintf('Number of points on the bnd: %g \n',length(bnd))
%END CREATE bnd INFORMATION

%CREATE COEFFICIENT MATRICES
poorc=0;
for n=1:nn
    if mesh(n).bnd>0
        continue
    end
    na=length(mesh(n).nbr);
    A=zeros(na,9);
    dx=mesh(n).dx; dy=mesh(n).dy;
    for i=1:na
        A(i,1:9)=[dx(i) dy(i) ...
            dx(i)*dy(i) dx(i)^2/2 dy(i)^2/2 ...
            dx(i)^2*dy(i)/6 dy(i)^2*dx(i)/6 dx(i)^3/6 dy(i)^3/6];
    end
    if na>9
        temp=(A'*A)\A'; %least squares 
        if rcond(inv((A'*A)))<1e-15
            poorc=poorc+1;
            disp(n)
            disp([mesh(n).dx' mesh(n).dy'])
        end
    elseif na>=5 && na<9
        A=A(:,1:5);
        if rcond((A'*A))<1e-15
            poorc=poorc+1;
            disp(n)
        end
        temp=(A'*A)\A';
    elseif na<5
        disp('neighbor connectivity less than 5 not supported\n')
        return
    else
        temp=inv(A);
    end
    mesh(n).coeff=temp(4,:)+temp(5,:); %dxx+dyy
end
for n=1:nn
    mesh(n).coeff=[mesh(n).coeff sum(mesh(n).coeff)];
end
fprintf('Number of poorly conditioned coefficient matrices: %g \n',poorc)
%END CREATE COEFFICIENT MATRICES

end %importmesh2