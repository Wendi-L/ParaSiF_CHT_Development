import os
from pathlib import Path

import basix.ufl
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem, io
from dolfinx.fem import petsc


comm = MPI.COMM_WORLD
rank = comm.rank

if rank == 0:
    print("********** THERMO-ELASTICITY SIMULATION BEGIN **********")


# -----------------------------------------------------------------------------
# Input/output
# -----------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent
MESH_FILE = ROOT / "sphere_validation_mesh" / "Structure_FEniCS.msh"
OUTPUT_DIR = ROOT / "Benchmark_results"
OUTPUT_DIR.mkdir(exist_ok=True)


# -----------------------------------------------------------------------------
# Material and thermal parameters
# -----------------------------------------------------------------------------

E_s = 190e9  # Young's modulus [Pa]
rho_s = 7750.0  # Density [kg/m^3]
nu_s = 0.305  # Poisson ratio [-]

T_ref_value = 0.0  # Reference temperature [K]
T_inner = 200.0  # Temperature on inner spherical surface [K]
T_outer = 0.0  # Temperature on outer spherical surface [K]

alpha = 9.7e-6  # Thermal expansion coefficient [1/K]
c_m = 486.0  # Specific heat capacity [m^2/(K*s^2)]
k_s = 20.0  # Thermal conductivity [kg*m/(K*s^3)]

Time = 2.0
dt_value = 1.0


# Physical tags from sphere_validation_mesh/Structure_FEniCS.geo
INNER_TAG = 1
OUTER_TAG = 2
Z_SYMMETRY_TAG = 3
Y_SYMMETRY_TAG = 4
X_SYMMETRY_TAG = 5


# -----------------------------------------------------------------------------
# Mesh
# -----------------------------------------------------------------------------

def read_gmsh_mesh(filename):
    try:
        # DOLFINx <= 0.7 style used by comet-fenicsx-main.
        gmshio = io.gmshio
    except AttributeError:
        # Newer DOLFINx style used by recent dolfinx-tutorial-main.
        from dolfinx.io import gmsh as gmshio

    try:
        mesh_data = gmshio.read_from_msh(str(filename), comm, rank=0, gdim=3)
    except TypeError:
        # Some releases do not expose the rank keyword on read_from_msh.
        mesh_data = gmshio.read_from_msh(str(filename), comm, gdim=3)

    if isinstance(mesh_data, tuple):
        return mesh_data
    return mesh_data.mesh, mesh_data.cell_tags, mesh_data.facet_tags


domain, cell_tags, facet_tags = read_gmsh_mesh(MESH_FILE)

if facet_tags is None:
    raise RuntimeError(f"No facet tags found in {MESH_FILE}")

tdim = domain.topology.dim
fdim = tdim - 1
gdim = domain.geometry.dim
domain.topology.create_connectivity(fdim, tdim)

dx = ufl.Measure("dx", domain=domain, subdomain_data=cell_tags)
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

if rank == 0:
    print(f"Loaded mesh: {domain.topology.index_map(tdim).size_global} cells")


# -----------------------------------------------------------------------------
# Function spaces
# -----------------------------------------------------------------------------

degree = 2
cell = domain.basix_cell()
u_element = basix.ufl.element("P", cell, degree, shape=(gdim,))
T_element = basix.ufl.element("P", cell, degree)
V = fem.functionspace(domain, basix.ufl.mixed_element([u_element, T_element]))

V_T, _ = V.sub(1).collapse()
V_u_components = [V.sub(0).sub(i).collapse()[0] for i in range(gdim)]

Q = fem.functionspace(domain, ("DG", 1))


# -----------------------------------------------------------------------------
# Boundary conditions
# -----------------------------------------------------------------------------

def constant_function(space, value):
    f = fem.Function(space)
    f.x.array[:] = default_scalar_type(value)
    return f


def component_bc(component, value, tag):
    subspace = V.sub(0).sub(component)
    collapsed = V_u_components[component]
    dofs = fem.locate_dofs_topological(
        (subspace, collapsed), fdim, facet_tags.find(tag)
    )
    return fem.dirichletbc(constant_function(collapsed, value), dofs, subspace)


def temperature_bc(value, tag):
    dofs = fem.locate_dofs_topological((V.sub(1), V_T), fdim, facet_tags.find(tag))
    return fem.dirichletbc(constant_function(V_T, value), dofs, V.sub(1))


bcs = [
    component_bc(0, 0.0, INNER_TAG),
    component_bc(1, 0.0, INNER_TAG),
    temperature_bc(T_inner, INNER_TAG),
    component_bc(0, 0.0, OUTER_TAG),
    component_bc(1, 0.0, OUTER_TAG),
    component_bc(2, 0.0, INNER_TAG),
    component_bc(2, 0.0, OUTER_TAG),
]


# -----------------------------------------------------------------------------
# Variational formulation
# -----------------------------------------------------------------------------

U = fem.Function(V, name="Thermoelastic solution")
U_old = fem.Function(V, name="Previous thermoelastic solution")
dU = ufl.TrialFunction(V)
U_test = ufl.TestFunction(V)

du, dT = ufl.split(dU)
u_old, T_old = ufl.split(U_old)
u_test, T_test = ufl.split(U_test)

I = ufl.Identity(gdim)

E = fem.Constant(domain, default_scalar_type(E_s))
nu = fem.Constant(domain, default_scalar_type(nu_s))
rho = fem.Constant(domain, default_scalar_type(rho_s))
cV = fem.Constant(domain, default_scalar_type(c_m))
k = fem.Constant(domain, default_scalar_type(k_s))
T_ref = fem.Constant(domain, default_scalar_type(T_ref_value))
dt = fem.Constant(domain, default_scalar_type(dt_value))

mu = E / (2.0 * (1.0 + nu))
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
mu_value = E_s / (2.0 * (1.0 + nu_s))
lmbda_value = E_s * nu_s / ((1.0 + nu_s) * (1.0 - 2.0 * nu_s))
kappa_value = alpha * (2.0 * mu_value + 3.0 * lmbda_value)
kappa = fem.Constant(domain, default_scalar_type(kappa_value))


def epsilon(v):
    return ufl.sym(ufl.grad(v))


def sigma(v, theta):
    return lmbda * ufl.tr(epsilon(v)) * I + 2.0 * mu * epsilon(v) - kappa * theta * I


def von_mises(v, theta):
    s = sigma(v, theta)
    dev = s - ufl.tr(s) / 3.0 * I
    return ufl.sqrt(3.0 / 2.0 * ufl.inner(dev, dev))


a = ufl.inner(sigma(du, dT), epsilon(u_test)) * dx
a += (rho * cV / dt) * dT * T_test * dx
a += (kappa * T_ref / dt) * ufl.tr(epsilon(du)) * T_test * dx
a += ufl.dot(k * ufl.grad(dT), ufl.grad(T_test)) * dx

L = (rho * cV / dt) * T_old * T_test * dx
L += (kappa * T_ref / dt) * ufl.tr(epsilon(u_old)) * T_test * dx

PETSC_OPTIONS = {
    "ksp_type": "preonly",
    "pc_type": "lu",
}
if os.environ.get("THERMOELASTIC_USE_MUMPS", "0") == "1":
    PETSC_OPTIONS["pc_factor_mat_solver_type"] = "mumps"

problem = petsc.LinearProblem(a, L, bcs=bcs, u=U, petsc_options=PETSC_OPTIONS)

u_h, T_h = ufl.split(U)
vm_expr = fem.Expression(von_mises(u_h, T_h), Q.element.interpolation_points())
vm_stress = fem.Function(Q, name="von_Mises stress")


# -----------------------------------------------------------------------------
# Time loop and output
# -----------------------------------------------------------------------------

vtk_u = io.VTKFile(comm, str(OUTPUT_DIR / "displacement.pvd"), "w")
vtk_T = io.VTKFile(comm, str(OUTPUT_DIR / "temperature.pvd"), "w")
vtk_vm = io.VTKFile(comm, str(OUTPUT_DIR / "stress.pvd"), "w")

t = 0.0
while t <= Time + 1e-12:
    if rank == 0:
        print(f"Time: {t:g}")

    problem.solve()
    U.x.scatter_forward()

    displacement = U.sub(0).collapse()
    displacement.name = "Displacement"
    temperature = U.sub(1).collapse()
    temperature.name = "Temperature"
    vm_stress.interpolate(vm_expr)
    vm_stress.x.scatter_forward()

    vtk_u.write_function(displacement, t)
    vtk_T.write_function(temperature, t)
    vtk_vm.write_function(vm_stress, t)

    U_old.x.array[:] = U.x.array
    U_old.x.scatter_forward()
    t += dt_value

vtk_u.close()
vtk_T.close()
vtk_vm.close()

if rank == 0:
    print()
    print("********** THERMO-ELASTICITY SIMULATION COMPLETED **********")
