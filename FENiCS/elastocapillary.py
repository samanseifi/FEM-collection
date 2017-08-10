from dolfin import *
import numpy as np

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

user_par = Parameters("user")
user_par.add("bounds_xmin",-0.5)
user_par.add("bounds_xmax",0.5)
user_par.add("bounds_ymin",-0.5)
user_par.add("bounds_ymax",0.5)
user_par.add("bounds_zmin",0.)
user_par.add("bounds_zmax", 2.5)
user_par.add("fe_order_u",1)
user_par.add("fe_order_p",1)
user_par.add("gamma_min",0.)
user_par.add("gamma_max",0.4)
user_par.add("gamma_nsteps",20)
user_par.add("mesh_ref",10)
user_par.add("save_dir","results")
user_par.add("output_type","pvd")
user_par.add("plot",True)

#xmin,xmax = user_par.bounds_xmin,user_par.bounds_xmax
#ymin,ymax = user_par.bounds_ymin,user_par.bounds_ymax
#zmin,zmax = user_par.bounds_zmin,user_par.bounds_zmax
#geom = Box(xmin,ymin,zmin,xmax,ymax,zmax)
#mesh = Mesh(geom,user_par.mesh_ref)

# Create mesh and define function space
mesh = UnitCubeMesh(16, 16, 16)
v1 = VectorFunctionSpace(mesh, "CG", 1)
p1 = FunctionSpace(mesh, "CG", 1)
V = MixedFunctionSpace([p1, v1])

V_u = V.sub(1)
V_p = V.sub(0)

ndim = 3

up = Function(V)
(p, u) = split(up)

dup = TrialFunction(V)
vq = TestFunction(V)
(q, v) = TestFunctions(V)

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# Define Dirichlet boundary (x = 0 or x = 1)
c = Expression(("0.0", "0.0", "0.0"))
r = Expression(("scale*0.0",
                "scale*0.0",
                "scale*0.0"),
                scale = 0.5, y0 = 0.5, z0 = 0.5, theta = pi/3)

bcl = DirichletBC(V_u, c, left)
#bcr = DirichletBC(V, r, right)
#bcs = [bcl, bcr]
bcs = []

# Define functions
#du = TrialFunction(V)            # Incremental displacement
#v  = TestFunction(V)             # Test function
#u  = Function(V)                 # Displacement from previous iteration
B  = Constant((0.0,  0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0,  0.0, 0.0))  # Traction force on the boundary

# Kinematics
d = u.geometric_dimension()
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

N = FacetNormal(mesh)

NansonOp = transpose(cofac(F))

deformed_N = dot(NansonOp,N)

current_element_of_area = sqrt(dot(deformed_N,deformed_N))

gamma=Expression("t",t=1.00)
surface_energy_density = gamma*current_element_of_area
surface_energy = surface_energy_density*ds

# Elasticity parameters
#E, nu = 5.0, 0.3
#mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))
mu, lmbda = 1.0, 1000.0

# Stored strain energy density (compressible neo-Hookean model)
#psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2
psi = mu*(Ic - 3) - (mu + p)*ln(J) - 1/(2*lmbda)*p**2

# Total potential energy
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds + surface_energy
Pi = psi*dx + surface_energy

# Compute first variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, up, vq)

# Compute Jacobian of F
J = derivative(F, up, dup)

# Solve variational problem


gamma_list = np.linspace(user_par.gamma_min,user_par.gamma_max,user_par.gamma_nsteps) # list of values of surface tension for the simulations
# directory and files to save the results
#save_dir = parameters.user.save_dir
#file_u = File(save_dir+"/displacement."+parameters.user.output_type)
#file_p = File(save_dir+"/pressure."+parameters.user.output_type)
# Solve with Newton solver for each value of the surface tension, using the previous solution as a starting point.
file_u = File("displacement.pvd")
file_p = File("pressure.pvd")

for t in gamma_list:
    # update the value of the surface tension
    gamma.t = t
    # solve the nonlinear problem (using Newton solver)
    solve(F == 0, up, bcs, J=J,
          form_compiler_parameters=ffc_options)
    # Save solution to file (readable by Paraview or Visit)
    (p,u) = up.split()
    file_u << (u,t)
    file_p << (p,t)
    # Plot and save png image
    plot_u = plot(u, mode = "displacement",title="Displacement field gamma=%.4f"%t,elevate=25.0)
    #    plot_u.write_png(save_dir+"/displacement_%.4f"%t

# Save solution in VTK format
#file = File("displacement.pvd");
#file << u;

# Plot and hold solution
plot(u, mode = "displacement", interactive = True)
