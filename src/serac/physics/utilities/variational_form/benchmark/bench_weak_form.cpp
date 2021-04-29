#include <fstream>
#include <iostream>

#include "mfem.hpp"

#include "axom/slic/core/SimpleLogger.hpp"

#include "serac/serac_config.hpp"
#include "serac/physics/operators/stdfunction_operator.hpp"
#include "serac/numerics/expr_template_ops.hpp"
#include "serac/physics/utilities/variational_form/weak_form.hpp"
#include "serac/physics/utilities/variational_form/tensor.hpp"
#include "serac/infrastructure/profiling.hpp"
#include "benchmark.hpp"

#include <gtest/gtest.h>

using namespace std;
using namespace mfem;

int         num_procs, myid;
int         refinements = 0;
const char* mesh_file   = SERAC_REPO_DIR "/data/meshes/star.mesh";
int         nsamples    = 1;

auto setup(int argc, char* argv[])
{
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&refinements, "-r", "--ref", "");
  args.AddOption(&nsamples, "-n", "--num-samples", "Samples per test");

  args.Parse();
  if (!args.Good()) {
    if (myid == 0) {
      args.PrintUsage(cout);
    }
    MPI_Finalize();
    exit(1);
  }
  if (myid == 0) {
    args.PrintOptions(cout);
  }

  mfem::Mesh mesh(mesh_file, 1, 1);
  for (int l = 0; l < refinements; l++) {
    mesh.UniformRefinement();
  }

  return mfem::ParMesh(MPI_COMM_WORLD, mesh);
}

template <int p, int dim>
void weak_form_test(mfem::ParMesh& mesh, H1<p> test, H1<p> trial, Dimension<dim>)
{
  using namespace serac::profiling;
  std::string postfix = serac::profiling::concat("_H1<", p, ">");

  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = H1_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec);

  ParGridFunction x(&fespace);
  x.Randomize();

  Vector X(fespace.TrueVSize());
  x.GetTrueDofs(X);

  ParBilinearForm A(&fespace);

  ConstantCoefficient a_coef(a);
  A.AddDomainIntegrator(new MassIntegrator(a_coef));

  ConstantCoefficient b_coef(b);
  A.AddDomainIntegrator(new DiffusionIntegrator(b_coef));
  {
    SERAC_PROFILE_REGION(concat("mfem_localAssemble", postfix));
    A.Assemble(0);
  }

  A.Finalize();

  mfem::Vector r1, r2;
  mfem::Vector g1, g2;

  SERAC_PROFILE_REGION2(concat("mfem_parallelAssemble", postfix))
  {
    std::unique_ptr<mfem::HypreParMatrix> J(A.ParallelAssemble());

    SERAC_PROFILE_REGION2(concat("mfem_ApplyGradient", postfix)) { g1 = (*J) * x; }
  }

  ParLinearForm       f(&fespace);
  FunctionCoefficient load_func([&](const Vector& coords) { return 100 * coords(0) * coords(1); });

  f.AddDomainIntegrator(new DomainLFIntegrator(load_func));
  {
    SERAC_PROFILE_REGION(concat("mfem_fAssemble", postfix));
    f.Assemble();
  }

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  auto constitutive_model = [&](auto x, auto temperature) {
    auto [u, du_dx] = temperature;
    auto f0         = a * u - (100 * x[0] * x[1]);
    auto f1         = b * du_dx;
    return std::tuple{f0, f1};
  };

  if constexpr (dim == 2) {
    residual.AddAreaIntegral(constitutive_model, mesh);
  }
  if constexpr (dim == 3) {
    residual.AddVolumeIntegral(constitutive_model, mesh);
  }

  {
    SERAC_PROFILE_REGION(concat("mfem_Apply", postfix));
    r1 = A * x - f;
  }

  {
    SERAC_PROFILE_REGION(concat("weakform_AssembleApply", postfix));
    r2 = residual * x;
  }

  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), 1.e-14);

  {
    SERAC_PROFILE_REGION(concat("weakform_GetGradient", postfix));
    mfem::Operator& grad2 = residual.GetGradient(x);
    {
      SERAC_PROFILE_REGION(concat("weakform_ApplyGradient", postfix));
      g2 = grad2 * x;
    }
  }

  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), 1.e-14);
}

template <int p, int dim>
void weak_form_test(mfem::ParMesh& mesh, H1<p, dim> test, H1<p, dim> trial, Dimension<dim>)
{
  using namespace serac::profiling;
  std::string             postfix = serac::profiling::concat("_H1<", p, ",", dim, ">");
  static constexpr double a       = 1.7;
  static constexpr double b       = 2.1;

  auto                  fec = H1_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec, dim);

  ParBilinearForm A(&fespace);

  ConstantCoefficient a_coef(a);
  A.AddDomainIntegrator(new VectorMassIntegrator(a_coef));

  ConstantCoefficient lambda_coef(b);
  ConstantCoefficient mu_coef(b);
  A.AddDomainIntegrator(new ElasticityIntegrator(lambda_coef, mu_coef));
  SERAC_PROFILE_EXPR(concat("mfem_localAssemble", postfix), A.Assemble(0));
  A.Finalize();

  SERAC_MARK_START(concat("mfem_parallelAssemble", postfix));
  std::unique_ptr<mfem::HypreParMatrix> J(A.ParallelAssemble());
  SERAC_MARK_END(concat("mfem_parallelAssemble", postfix));

  LinearForm                f(&fespace);
  VectorFunctionCoefficient load_func(dim, [&](const Vector& /*coords*/, Vector& force) {
    force    = 0.0;
    force(0) = -1.0;
  });

  f.AddDomainIntegrator(new VectorDomainLFIntegrator(load_func));

  SERAC_PROFILE_EXPR(concat("mfem_fAssemble", postfix), f.Assemble());

  ParGridFunction x(&fespace);
  x.Randomize();

  Vector X(fespace.TrueVSize());
  x.GetTrueDofs(X);

  static constexpr auto I = Identity<dim>();

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  auto constitutive_model = [&](auto /*x*/, auto displacement) {
    auto [u, du_dx] = displacement;
    auto f0         = a * u + I[0];
    auto strain     = 0.5 * (du_dx + transpose(du_dx));
    auto f1         = b * tr(strain) * I + 2.0 * b * strain;
    return std::tuple{f0, f1};
  };

  if constexpr (dim == 2) {
    residual.AddAreaIntegral(constitutive_model, mesh);
  }

  if constexpr (dim == 3) {
    residual.AddVolumeIntegral(constitutive_model, mesh);
  }

  mfem::Vector r1 = SERAC_PROFILE_EXPR(concat("mfem_Apply", postfix), A * x - f);

  mfem::Vector r2 = SERAC_PROFILE_EXPR(concat("weakform_AssembleApply", postfix), residual * x);

  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), 1.e-14);

  mfem::Vector g1 = SERAC_PROFILE_EXPR(concat("mfem_ApplyGradient", postfix), (*J) * x);

  mfem::Operator& grad = SERAC_PROFILE_EXPR(concat("weakform_GetGradient", postfix), residual.GetGradient(x));

  mfem::Vector g2 = SERAC_PROFILE_EXPR(concat("weakform_ApplyGradient", postfix), grad * x);

  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), 1.e-14);
}

template <int p, int dim>
void weak_form_test(mfem::ParMesh& mesh, Hcurl<p> test, Hcurl<p> trial, Dimension<dim>)
{
  using namespace serac::profiling;
  std::string postfix = serac::profiling::concat("_Hcurl<", p, ">");

  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = ND_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec);

  ParBilinearForm A(&fespace);

  ConstantCoefficient a_coef(a);
  A.AddDomainIntegrator(new VectorFEMassIntegrator(a_coef));

  ConstantCoefficient b_coef(b);
  A.AddDomainIntegrator(new CurlCurlIntegrator(b_coef));
  SERAC_PROFILE_REGION2(concat("mfem_localAssemble", postfix))
  A.Assemble(0);
  A.Finalize();

  std::unique_ptr<mfem::HypreParMatrix> J(
      SERAC_PROFILE_EXPR(concat("mfem_parallelAssemble", postfix), A.ParallelAssemble()));

  LinearForm                f(&fespace);
  VectorFunctionCoefficient load_func(dim, [&](const Vector& coords, Vector& output) {
    double x  = coords(0);
    double y  = coords(1);
    output    = 0.0;
    output(0) = 10 * x * y;
    output(1) = -5 * (x - y) * y;
  });

  f.AddDomainIntegrator(new VectorFEDomainLFIntegrator(load_func));
  SERAC_PROFILE_REGION2(concat("mfem_fAssemble", postfix))
  f.Assemble();

  ParGridFunction x(&fespace);
  x.Randomize();

  Vector X(fespace.TrueVSize());
  x.GetTrueDofs(X);

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  residual.AddDomainIntegral(
      Dimension<dim>{},
      [&](auto x, auto vector_potential) {
        auto [A, curl_A] = vector_potential;
        auto f0          = a * A - tensor<double, dim>{10 * x[0] * x[1], -5 * (x[0] - x[1]) * x[1]};
        auto f1          = b * curl_A;
        return std::tuple{f0, f1};
      },
      mesh);

  mfem::Vector r1, r2;
  SERAC_PROFILE_REGION2(concat("mfem_Apply", postfix))
  r1 = A * x - f;
  SERAC_PROFILE_REGION2(concat("weakform_AssembleApply", postfix))
  r2 = residual * x;

  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), 1.e-13);

  mfem::Vector g2;
  mfem::Vector g1 = SERAC_PROFILE_EXPR(concat("mfem_ApplyGradient", postfix), (*J) * x);

  SERAC_PROFILE_REGION2(concat("weakform_GetGradient", postfix))
  {
    mfem::Operator& grad = residual.GetGradient(x);

    SERAC_PROFILE_REGION2(concat("weakform_ApplyGradient", postfix)) { g2 = grad * x; }
  }

  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), 1.e-13);
}

template <int dim>
void run_tests(mfem::ParMesh& mesh)
{
  Dimension<dim> d;

  for (int i = 0; i < nsamples; i++) {
    weak_form_test(mesh, H1<1>{}, H1<1>{}, d);
    weak_form_test(mesh, H1<2>{}, H1<2>{}, d);
    weak_form_test(mesh, H1<3>{}, H1<3>{}, d);

    weak_form_test(mesh, H1<1, dim>{}, H1<1, dim>{}, d);
    weak_form_test(mesh, H1<2, dim>{}, H1<2, dim>{}, d);
    weak_form_test(mesh, H1<3, dim>{}, H1<3, dim>{}, d);

    weak_form_test(mesh, Hcurl<1>{}, Hcurl<1>{}, d);
    weak_form_test(mesh, Hcurl<2>{}, Hcurl<2>{}, d);
    weak_form_test(mesh, Hcurl<3>{}, Hcurl<3>{}, d);
  }
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  axom::slic::SimpleLogger logger;

  serac::profiling::initializeCaliper();

  auto mesh = setup(argc, argv);

  if (mesh.Dimension() == 2) {
    run_tests<2>(mesh);
  }
  if (mesh.Dimension() == 3) {
    run_tests<3>(mesh);
  }

  serac::profiling::terminateCaliper();

  MPI_Finalize();

  return 0;
}
