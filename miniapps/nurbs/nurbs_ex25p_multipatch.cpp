//                     MFEM Example nurbs_ex25p_multipatch
//
// Compile with: make nurbs_ex25p_multipatch
//
// Sample runs:  mpirun -np 4 nurbs_ex25p_multipatch -m meshes/two-squares-nurbs.mesh
//               mpirun -np 4 nurbs_ex25p_multipatch -m meshes/two-cubes-nurbs.mesh
//
// Description:  This example tests multi-patch H(curl) NURBS support.
//               It solves the Maxwell equation on a multi-patch NURBS mesh:
//
//                  curl curl E + E = f
//
//               with homogeneous Dirichlet boundary conditions E x n = 0.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Exact solution (for verification)
void E_exact(const Vector &x, Vector &E);
void f_exact(const Vector &x, Vector &f);

int main(int argc, char *argv[])
{
   // Initialize MPI
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   
   // Parse command-line options
   const char *mesh_file = "meshes/two-squares-nurbs.mesh";
   int order = 2;
   int ref_levels = 1;
   bool visualization = true;
   
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element order.");
   args.AddOption(&ref_levels, "-r", "--refine", "Number of refinements.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization", "Enable visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0) { args.PrintUsage(cout); }
      return 1;
   }
   if (myid == 0) { args.PrintOptions(cout); }

   // Read the mesh
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   
   if (myid == 0)
   {
      cout << "Mesh dimension: " << dim << endl;
      cout << "Number of NURBS patches: " << mesh->NURBSext->GetNP() << endl;
   }

   // Refine the mesh
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   // Create parallel mesh
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   if (myid == 0)
   {
      cout << "Parallel mesh created with " << num_procs << " processors." << endl;
      cout << "Testing multi-patch H(curl) NURBS..." << endl;
   }

   // Define finite element collection and space
   // This is the key test: H(curl) on multi-patch NURBS
   FiniteElementCollection *fec = new NURBS_HCurlFECollection(order, dim);
   NURBSExtension *NURBSext = new NURBSExtension(pmesh->NURBSext, order);
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, NURBSext, fec);
   
   HYPRE_BigInt global_size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << global_size << endl;
   }

   // Define boundary conditions (essential BC: E x n = 0 on all boundaries)
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // Set up the linear form (right-hand side)
   VectorFunctionCoefficient f_coef(dim, f_exact);
   ParLinearForm *b = new ParLinearForm(fespace);
   b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f_coef));
   b->Assemble();

   // Define the solution vector
   ParGridFunction x(fespace);
   x = 0.0;

   // Set up the bilinear form
   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new CurlCurlIntegrator);
   a->AddDomainIntegrator(new VectorFEMassIntegrator);
   a->Assemble();

   // Form the linear system
   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   if (myid == 0)
   {
      cout << "Linear system assembled." << endl;
   }

   // Solve the system using PCG with AMS preconditioner
   HypreAMS *ams = new HypreAMS(*A.As<HypreParMatrix>(), fespace);
   
   CGSolver cg(MPI_COMM_WORLD);
   cg.SetRelTol(1e-12);
   cg.SetAbsTol(1e-12);
   cg.SetMaxIter(500);
   cg.SetPrintLevel(1);
   cg.SetPreconditioner(*ams);
   cg.SetOperator(*A);
   cg.Mult(B, X);

   // Recover the solution
   a->RecoverFEMSolution(X, *b, x);

   // Compute and print the L2 error
   VectorFunctionCoefficient E_coef(dim, E_exact);
   real_t err = x.ComputeL2Error(E_coef);
   if (myid == 0)
   {
      cout << "\n|| E_h - E ||_{L^2} = " << err << endl;
   }

   // Save the solution for visualization
   if (visualization)
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << myid;
      sol_name << "sol." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }

   // Clean up
   delete ams;
   delete a;
   delete b;
   delete fespace;
   delete fec;
   delete pmesh;

   if (myid == 0)
   {
      cout << "\nMulti-patch H(curl) test completed successfully!" << endl;
   }

   return 0;
}

// Exact solution: E = (sin(pi*y), sin(pi*x), 0) for 2D
//                 E = (sin(pi*y)*sin(pi*z), sin(pi*x)*sin(pi*z), sin(pi*x)*sin(pi*y)) for 3D
void E_exact(const Vector &x, Vector &E)
{
   int dim = x.Size();
   E.SetSize(dim);
   
   if (dim == 2)
   {
      E(0) = sin(M_PI * x(1));
      E(1) = sin(M_PI * x(0));
   }
   else // dim == 3
   {
      E(0) = sin(M_PI * x(1)) * sin(M_PI * x(2));
      E(1) = sin(M_PI * x(0)) * sin(M_PI * x(2));
      E(2) = sin(M_PI * x(0)) * sin(M_PI * x(1));
   }
}

// Right-hand side: f = curl curl E + E
void f_exact(const Vector &x, Vector &f)
{
   int dim = x.Size();
   f.SetSize(dim);
   
   real_t pi2 = M_PI * M_PI;
   
   if (dim == 2)
   {
      // For 2D: curl curl E + E
      // curl E = dE_y/dx - dE_x/dy = pi*cos(pi*x) - pi*cos(pi*y)
      // curl curl E = (d(curl E)/dy, -d(curl E)/dx)
      //             = (pi^2*sin(pi*y), pi^2*sin(pi*x))
      // f = curl curl E + E = ((1+pi^2)*sin(pi*y), (1+pi^2)*sin(pi*x))
      f(0) = (1.0 + pi2) * sin(M_PI * x(1));
      f(1) = (1.0 + pi2) * sin(M_PI * x(0));
   }
   else // dim == 3
   {
      // For 3D, the computation is more complex
      // Using the fact that for E = (E_x, E_y, E_z):
      // curl curl E = grad(div E) - Laplacian(E)
      // Since div E = 0 for our choice, curl curl E = -Laplacian(E)
      f(0) = (1.0 + 2.0*pi2) * sin(M_PI * x(1)) * sin(M_PI * x(2));
      f(1) = (1.0 + 2.0*pi2) * sin(M_PI * x(0)) * sin(M_PI * x(2));
      f(2) = (1.0 + 2.0*pi2) * sin(M_PI * x(0)) * sin(M_PI * x(1));
   }
}
