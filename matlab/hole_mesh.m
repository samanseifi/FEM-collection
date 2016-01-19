function [mesh] = hole_mesh(nh, side_len, hole_radius)
% MAKE_MESH Returns a mesh for a square plate with a center hole.
%   Input arguments: nh - number of elements along half of a side.
%                    side_length - length of a side of the plate.
%                    hole_radius - radius of hole in center of plate.	    
    n_shell = 2*nh;
    nns = 8*nh;  % number of nodes in a shell.            
    circle = hole_radius*make_circle(nns)';
    square = side_len*make_square(nh)';
    mesh.x = [];
    for s = 1:n_shell
        % f goes from 0 to 1 as a linear function.
        f = (s-1)/(n_shell-1);
        mesh.x = [mesh.x, square*(1-f) + circle*f];
    end            
    mesh.conn = [];
    for j=1:n_shell-1
        for i=1:nns-1               
          n0 = i + (j-1)*nns;          
          mesh.conn(:,end+1) = [n0; n0+1; n0+1+nns; n0+nns];           
       end               
       mesh.conn(:,end+1) = [1+(j-1)*nns; 1+(j)*nns; 
                             mesh.conn(3,end); mesh.conn(2,end)];
    end    
    fprintf('Created mesh with %d elements, and %d nodes.\n', ...
        length(mesh.conn), length(mesh.x));
end

function [x] = make_square(nh)    
    x(1:nh+1,1) = 0.5;
    x(1:nh+1,2) = linspace(0, 0.5, nh+1);
    x(nh+1:3*nh+1,1) = linspace(0.5, -0.5, 2*nh+1);
    x(nh+1:3*nh+1,2) = 0.5;    
    x(3*nh+1:5*nh+1,1) = -0.5;
    x(3*nh+1:5*nh+1,2) = linspace(0.5, -0.5, 2*nh+1);    
    x(5*nh+1:7*nh+1,1) = linspace(-0.5, 0.5, 2*nh+1);
    x(5*nh+1:7*nh+1,2) = -0.5;
    x(7*nh+1:8*nh,1) = 0.5;    
    x(7*nh+1:8*nh,2) = linspace(-0.5,-0.5/nh,nh);
end

function [x] = make_circle(nns)
    q = linspace(0, 2*pi, nns+1)';
    x = [cos(q(1:end-1)), sin(q(1:end-1)) + 0.00*sin(7*q(1:end-1)) ];
end