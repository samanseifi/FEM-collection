% femcode2.m
 clear all
% [p,t,b] from distmesh tool
% make sure your matlab path includes the directory where distmesh is installed.
 
fd=@(p) sqrt(sum(p.^2,2))-1;
[p,t]=distmesh2d(fd,@huniform,0.05,[-1,-1;1,1],[]);
b1=unique(boundedges(p,t));

%fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.2));
%fh=@(p) 0.05+0.3*dcircle(p,0,0,0.2);
%[p,t]=distmesh2d(fd,fh,0.025,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
       
%[b,b2,b3]=unique(boundedges(p,t));
       
%[xyz,mesh,nn,tri,bnd,map]=importmesh2('plate_hole.msh');
%p = xyz;
%t = tri;
%b = bnd;



% [K,F] = assemble(p,t) % K and F for any mesh of triangles: linear phi's
N = size(p, 1);
T = size(t, 1); % number of nodes, number of triangles
% p lists x,y coordinates of N nodes, t lists triangles by 3 node numbers
K = sparse(2*N, 2*N); % zero matrix in sparse format: zeros(N) would be "dense"
F = zeros(2*N, 1); % load vector F to hold integrals of phi's times load f(x,y)
 
for e = 1:T  % integration over one triangular element at a time
    nodes = t(e, :); % row of t = node numbers of the 3 corners of triangle e
    pe = p(nodes,:); % positions
    x = pe(:, 1);
    y = pe(:, 2);
    Pe = [ones(3,1), pe]; % 3 by 3 matrix with rows=[1 xcorner ycorner]
    Area = abs(det(Pe))/2; % area of triangle e = half of parallelogram area
    C = inv(Pe); % columns of C are coeffs in a+bx+cy to give phi=1,0,0 at nodes
	% now compute 3 by 3 Ke and 3 by 1 Fe for element e
    B = [y(2)-y(3), 0, y(3)-y(1), 0, y(1)-y(2),0;
        0, x(3)-x(2), 0, x(1)-x(3), 0, x(2)-x(1); 
        x(2)-x(3), y(3)-y(2), x(3)-x(1), y(1)-y(3), x(1)-x(2), y(2)-y(1)];    
    B = B/2/Area;
    Ke = Area*(B'*B); % element matrix from slopes b,c in grad
    Fe = Area/6*4; % intral of phi over triangle is volume of pyramid: f(x,y)=4
    % multiply Fe by f at centroid for load f(x,y): one-point quadrature!
    % centroid would be mean(p(nodes,:)) = average of 3 node coordinates
    n = [2*nodes(1)-1  2*nodes(1)  2*nodes(2)-1  2*nodes(2)  2*nodes(3)-1  2*nodes(3)];
    
    K(n,n) = K(n,n) + Ke; % add Ke to 9 entries of global K
    F(n) = F(n) + Fe; % add Fe to 3 components of load vector F
end   % all T element matrices and vectors now assembled into K and F
 
% [Kb,Fb] = dirichlet(K,F,b) % assembled K was singular! K*ones(N,1)=0
% Implement Dirichlet boundary conditions U(b)=0 at nodes in list b
K(2*b,:) = 0; K(2*b-1,:) = 0; K(:,2*b-1) = 0; K(:,2*b) = 0; F(2*b) = 0; F(2*b-1)=0;% put zeros in boundary rows/columns of K and F

K(2*b,2*b) = speye(length(b), length(b)); % put I into boundary submatrix of K
K(2*b-1,2*b-1) = speye(length(b), length(b)); % put I into boundary submatrix of K


Kb = K; Fb = F; % Stiffness matrix Kb (sparse format) and load vector Fb
 
% Solving for the vector U will produce U(b)=0 at boundary nodes
U=Kb\Fb;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N
%U = gmres(Kb, Fb, 20, 1e-9, 30);
% Plot the FEM approximation U(x,y) with values U_1 to U_N at the nodes
trisurf(t, p(:, 1), p(:, 2), 0*p(:,1), sqrt(U(2:2:2*N).^2 + U(1:2:2*N).^2 ),'edgecolor','k','facecolor','interp');
%axis([-1 1 -1 1]), colorbar
view(0, 90)