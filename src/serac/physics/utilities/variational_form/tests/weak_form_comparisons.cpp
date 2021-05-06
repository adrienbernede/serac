#include <fstream>
#include <iostream>

#include "mfem.hpp"

#include "axom/slic/core/SimpleLogger.hpp"

#include "serac/serac_config.hpp"
#include "serac/numerics/mesh_utils.hpp"
#include "serac/numerics/expr_template_ops.hpp"
#include "serac/physics/operators/stdfunction_operator.hpp"
#include "serac/physics/utilities/variational_form/weak_form.hpp"
#include "serac/physics/utilities/variational_form/tensor.hpp"
#include "serac/physics/utilities/variational_form/integral.hpp"

#include <gtest/gtest.h>

using namespace std;
using namespace mfem;
using namespace serac;

double tol = 1.e-13;
int    num_procs, myid;

constexpr bool                 verbose = false;
std::unique_ptr<mfem::ParMesh> mesh2D;
std::unique_ptr<mfem::ParMesh> mesh3D;

template <auto offsetV, auto indicesV, auto gatherV>
struct forbidden_restriction {
  //  friend result_fes _get_fes(From & from) { return & from.*fesV; }
  friend mfem::Array<int>& __get_offsets(mfem::ElementRestriction& from) { return from.*offsetV; }
  friend mfem::Array<int>& __get_indices(mfem::ElementRestriction& from) { return from.*indicesV; }
  friend mfem::Array<int>& __get_gatherMap(mfem::ElementRestriction& from) { return from.*gatherV; }
};

// const FiniteElementSpace * __get_fes(mfem::ElementRestriction &);
mfem::Array<int>& __get_offsets(mfem::ElementRestriction&);
mfem::Array<int>& __get_indices(mfem::ElementRestriction&);
mfem::Array<int>& __get_gatherMap(mfem::ElementRestriction&);

template struct forbidden_restriction<&mfem::ElementRestriction::offsets, &mfem::ElementRestriction::indices,
                                      &mfem::ElementRestriction::gatherMap>;

namespace serac {
namespace mfem_ext {
class AssembledSparseMatrix : public mfem::SparseMatrix {
public:
  AssembledSparseMatrix(const mfem::FiniteElementSpace& test,   // test_elem_dofs * ne * vdim x vdim * test_ndofs
                        const mfem::FiniteElementSpace& trial,  // trial_elem_dofs * ne * vdim x vdim * trial_ndofs
                        mfem::ElementDofOrdering        elem_order)
      : mfem::SparseMatrix(test.GetNDofs() * test.GetVDim(), trial.GetNDofs() * trial.GetVDim()),
        test_fes(test),
        trial_fes(trial),
        test_restriction(test, elem_order),
        trial_restriction(trial, elem_order),
        elem_ordering(elem_order)
  {
    GetMemoryI().New(Height() + 1, GetMemoryI().GetMemoryType());

    const int nnz = FillI();
    GetMemoryJ().New(nnz, GetMemoryJ().GetMemoryType());
    GetMemoryData().New(nnz, GetMemoryData().GetMemoryType());
    FillJ();

    // zero initialize the data
    for (int i = 0; i < nnz; i++) {
      A[i] = 0.;
    }
  }

  int  FillI();
  void FillJ();
  void FillData(const mfem::Vector& ea_data);

protected:
  const mfem::FiniteElementSpace& test_fes;
  const mfem::FiniteElementSpace& trial_fes;
  mfem::ElementRestriction        test_restriction;
  mfem::ElementRestriction        trial_restriction;
  mfem::ElementDofOrdering        elem_ordering;
  mfem::Array<int>                ea_map;
};

int AssembledSparseMatrix::FillI()
{
  [[maybe_unused]] auto& test_offsets = __get_offsets(test_restriction);  // offsets for rows.. each row is a test_vdof
  [[maybe_unused]] auto& test_indices =
      __get_indices(test_restriction);  // returns (test_elem_dof , ne) id corresponding to a test_vdof_offset
  [[maybe_unused]] auto& test_gatherMap  = __get_gatherMap(test_restriction);  // returns test_vdof
  [[maybe_unused]] auto& trial_offsets   = __get_offsets(trial_restriction);
  [[maybe_unused]] auto& trial_indices   = __get_indices(trial_restriction);
  [[maybe_unused]] auto& trial_gatherMap = __get_gatherMap(trial_restriction);

  /**
     We expect mat_ea to be of size (test_elem_dof * test_vdim, trial_elem_dof * trial_vdim, ne)
     We assume a consistent striding from (elem_dof, vd) within each element

   */
  const int test_elem_dof  = test_fes.GetFE(0)->GetDof();
  const int trial_elem_dof = trial_fes.GetFE(0)->GetDof();
  const int test_vdim      = test_fes.GetVDim();
  const int trial_vdim     = trial_fes.GetVDim();
  const int test_ndofs     = test_fes.GetNDofs();

  auto I = ReadWriteI();
  for (int i = 0; i < test_vdim * test_ndofs; i++) {
    I[i] = 0;
  }

  for (int test_vdof = 0; test_vdof < test_fes.GetNDofs(); test_vdof++) {
    // Look through each element corresponding to a test_vdof
    const int test_row_offset = test_offsets[test_vdof];
    const int nrow_elems      = test_offsets[test_vdof + 1] - test_row_offset;

    // Build temporary array to get rid of duplicates
    mfem::Array<int> trial_vdofs(nrow_elems * trial_elem_dof);
    trial_vdofs = -1;
    int nnz_row = 0;
    for (int e_index = 0; e_index < nrow_elems; e_index++) {
      const int                  test_offset = test_indices[test_row_offset + e_index];
      const int                  e           = test_offset / test_elem_dof;
      [[maybe_unused]] const int test_i_elem = test_offset % test_elem_dof;

      // find corresponding trial_vdofs
      mfem::Array<int> trial_elem_vdofs(trial_elem_dof);
      for (int j = 0; j < trial_elem_dof; j++) {
        const auto trial_j_vdof = trial_gatherMap[trial_elem_dof * e + j];
        trial_elem_vdofs[j]     = trial_j_vdof;
        if (trial_vdofs.Find(trial_j_vdof) == -1) {
          // we haven't seen this before
          trial_vdofs[nnz_row] = trial_j_vdof;
          nnz_row++;
        }
      }
    }

    // add entries to I
    for (int vi = 0; vi < test_vdim; vi++) {
      I[test_fes.DofToVDof(test_vdof, vi)] = nnz_row * trial_vdim;
    }
  }

  // Perform inclusive scan on all entries
  int nnz = 0;
  for (int i = 0; i < test_ndofs * trial_vdim; i++) {
    int temp = I[i];
    I[i]     = nnz;
    nnz += temp;
  }
  I[test_ndofs * trial_vdim] = nnz;

  return nnz;
}

void AssembledSparseMatrix::FillJ()
{
  auto I = ReadWriteI();
  auto J = WriteJ();

  [[maybe_unused]] auto& test_offsets = __get_offsets(test_restriction);  // offsets for rows.. each row is a test_vdof
  [[maybe_unused]] auto& test_indices =
      __get_indices(test_restriction);  // returns (test_elem_dof , ne) id corresponding to a test_vdof_offset
  [[maybe_unused]] auto& test_gatherMap  = __get_gatherMap(test_restriction);  // returns test_vdof
  [[maybe_unused]] auto& trial_offsets   = __get_offsets(trial_restriction);
  [[maybe_unused]] auto& trial_indices   = __get_indices(trial_restriction);
  [[maybe_unused]] auto& trial_gatherMap = __get_gatherMap(trial_restriction);

  const int                  test_elem_dof  = test_fes.GetFE(0)->GetDof();
  const int                  trial_elem_dof = trial_fes.GetFE(0)->GetDof();
  const int                  test_vdim      = test_fes.GetVDim();
  const int                  trial_vdim     = trial_fes.GetVDim();
  [[maybe_unused]] const int test_ndofs     = test_fes.GetNDofs();

  const int ne = trial_fes.GetNE();
  ea_map.SetSize(test_elem_dof * test_vdim * trial_elem_dof * trial_vdim * ne);
  auto map_ea = Reshape(ea_map.ReadWrite(), test_elem_dof * test_vdim, trial_elem_dof * trial_vdim, ne);

  // initialize J
  for (int j = 0; j < this->J.Capacity(); j++) {
    this->J[j] = -1;
  }

  for (int test_vdof = 0; test_vdof < test_fes.GetNDofs(); test_vdof++) {
    // Look through each element corresponding to a test_vdof
    const int test_row_offset = test_offsets[test_vdof];
    const int nrow_elems      = test_offsets[test_vdof + 1] - test_row_offset;

    // here we assume all the components have the same number of columns
    const int        nnz_row = I[test_fes.DofToVDof(test_vdof, 0) + 1] - I[test_fes.DofToVDof(test_vdof, 0)];
    mfem::Array<int> trial_vdofs(nnz_row);
    trial_vdofs      = -1;
    int j_vdof_index = 0;

    // Build temporary array for assembled J
    for (int e_index = 0; e_index < nrow_elems; e_index++) {
      const int                  test_offset = test_indices[test_row_offset + e_index];
      const int                  e           = test_offset / test_elem_dof;
      [[maybe_unused]] const int test_i_elem = test_offset % test_elem_dof;

      // find corresponding trial_vdofs
      mfem::Array<int> trial_elem_vdofs(trial_elem_dof);
      for (int j_elem = 0; j_elem < trial_elem_dof; j_elem++) {
        const auto trial_j_vdof  = trial_gatherMap[trial_elem_dof * e + j_elem];
        trial_elem_vdofs[j_elem] = trial_j_vdof;

        auto find_index = trial_vdofs.Find(trial_j_vdof);
        if (find_index == -1) {
          // we haven't seen this before
          trial_vdofs[j_vdof_index] = trial_j_vdof;

          // we can add this entry to J
          for (int vi = 0; vi < test_vdim; vi++) {
            const auto i_dof_offset = I[test_fes.DofToVDof(test_vdof, vi)];

            // this access pattern corresnponds to j_vdof_index + vj * nnz_row
            for (int vj = 0; vj < trial_vdim; vj++) {
              const auto column_index = j_vdof_index + vj * nnz_row / trial_vdim;
              const auto j_nnz_index  = i_dof_offset + column_index;
              J[j_nnz_index]          = trial_fes.DofToVDof(trial_vdofs[j_vdof_index], vj);
            }
          }

          // write mapping from ea to csr_nnz_index (can probably optimize this)
          for (int vi = 0; vi < test_vdim; vi++) {
            const auto i_dof_offset = I[test_fes.DofToVDof(test_vdof, vi)];
            for (int vj = 0; vj < trial_vdim; vj++) {
              const auto column_index = j_vdof_index + vj * nnz_row / trial_vdim;
              map_ea(test_i_elem + test_elem_dof * vi, j_elem + trial_elem_dof * vj, e) = i_dof_offset + column_index;
            }
          }

          j_vdof_index++;
        } else {
          // this is a duplicate entry
          // write mapping from ea to csr_nnz_index (can probably optimize this)
          for (int vi = 0; vi < test_vdim; vi++) {
            const auto i_dof_offset = I[test_fes.DofToVDof(test_vdof, vi)];
            for (int vj = 0; vj < trial_vdim; vj++) {
              const auto column_index = find_index + vj * nnz_row / trial_vdim;
              map_ea(test_i_elem + test_elem_dof * vi, j_elem + trial_elem_dof * vj, e) = i_dof_offset + column_index;
            }
          }
        }
      }
    }
  }
}
void AssembledSparseMatrix::FillData(const mfem::Vector& ea_data)
{
  auto Data = WriteData();

  [[maybe_unused]] auto& test_offsets = __get_offsets(test_restriction);  // offsets for rows.. each row is a test_vdof
  [[maybe_unused]] auto& test_indices =
      __get_indices(test_restriction);  // returns (test_elem_dof , ne) id corresponding to a test_vdof_offset
  [[maybe_unused]] auto& test_gatherMap  = __get_gatherMap(test_restriction);  // returns test_vdof
  [[maybe_unused]] auto& trial_offsets   = __get_offsets(trial_restriction);
  [[maybe_unused]] auto& trial_indices   = __get_indices(trial_restriction);
  [[maybe_unused]] auto& trial_gatherMap = __get_gatherMap(trial_restriction);

  const int                  test_elem_dof  = test_fes.GetFE(0)->GetDof();
  const int                  trial_elem_dof = trial_fes.GetFE(0)->GetDof();
  const int                  test_vdim      = test_fes.GetVDim();
  const int                  trial_vdim     = trial_fes.GetVDim();
  [[maybe_unused]] const int test_ndofs     = test_fes.GetNDofs();

  const int ne = trial_fes.GetNE();

  auto map_ea = Reshape(ea_map.Read(), test_elem_dof * test_vdim, trial_elem_dof * trial_vdim, ne);

  auto mat_ea = Reshape(ea_data.Read(), test_elem_dof * test_vdim, trial_elem_dof * trial_vdim, ne);

  // Use map_ea to take ea_data directly to CSR entry
  for (int e = 0; e < ne; e++) {
    for (int i_elem = 0; i_elem < test_elem_dof; i_elem++) {
      for (int vi = 0; vi < test_vdim; vi++) {
        for (int j_elem = 0; j_elem < trial_elem_dof; j_elem++) {
          for (int vj = 0; vj < trial_vdim; vj++) {
            Data[map_ea(i_elem + vi * test_elem_dof, j_elem + vj * trial_elem_dof, e)] +=
                mat_ea(i_elem + vi * test_elem_dof, j_elem + vj * trial_elem_dof, e);
          }
        }
      }
    }
  }
}

}  // namespace mfem_ext

}  // namespace serac


// this test sets up a toy "thermal" problem where the residual includes contributions
// from a temperature-dependent source term and a temperature-gradient-dependent flux
// 
// the same problem is expressed with mfem and weak_form, and their residuals and gradient action
// are compared to ensure the implementations are in agreement.

template <int p, int dim>
void weak_form_test(mfem::ParMesh& mesh, H1<p> test, H1<p> trial, Dimension<dim>)
{
  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = H1_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec);

  ParBilinearForm A(&fespace);

  ConstantCoefficient a_coef(a);
  A.AddDomainIntegrator(new MassIntegrator(a_coef));

  ConstantCoefficient b_coef(b);
  A.AddDomainIntegrator(new DiffusionIntegrator(b_coef));
  A.Assemble(0);
  A.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> J(A.ParallelAssemble());

  ParLinearForm       f(&fespace);
  FunctionCoefficient load_func([&](const Vector& coords) { return 100 * coords(0) * coords(1); });

  f.AddDomainIntegrator(new DomainLFIntegrator(load_func));
  f.Assemble();
  std::unique_ptr<mfem::HypreParVector> F(f.ParallelAssemble());

  ParGridFunction u_global(&fespace);
  u_global.Randomize();

  Vector U(fespace.TrueVSize());
  u_global.GetTrueDofs(U);

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  residual.AddDomainIntegral(
      Dimension<dim>{},
      [&](auto x, auto temperature) {
        auto [u, du_dx] = temperature;
        auto source = a * u - (100 * x[0] * x[1]);
        auto flux = b * du_dx;
        return std::tuple{source, flux};
      },
      mesh);

  mfem::Vector r1 = (*J) * U - (*F);
  mfem::Vector r2 = residual(U);

  if (verbose) {
    std::cout << "||r1||: " << r1.Norml2() << std::endl;
    std::cout << "||r2||: " << r2.Norml2() << std::endl;
    std::cout << "||r1-r2||/||r1||: " << mfem::Vector(r1 - r2).Norml2() / r1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), tol);

  mfem::Operator& grad2 = residual.GetGradient(U);

  mfem::Vector g1 = (*J) * U;
  mfem::Vector g2 = grad2 * U;

  if (verbose) {
    std::cout << "||g1||: " << g1.Norml2() << std::endl;
    std::cout << "||g2||: " << g2.Norml2() << std::endl;
    std::cout << "||g1-g2||/||g1||: " << mfem::Vector(g1 - g2).Norml2() / g1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), tol);
}

// this test sets up a toy "elasticity" problem where the residual includes contributions
// from a displacement-dependent body force term and an isotropically linear elastic stress response
// 
// the same problem is expressed with mfem and weak_form, and their residuals and gradient action
// are compared to ensure the implementations are in agreement.
template <int p, int dim>
void weak_form_test(mfem::ParMesh& mesh, H1<p, dim> test, H1<p, dim> trial, Dimension<dim>)
{
  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = H1_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec, dim);

  ParBilinearForm A(&fespace);

  ConstantCoefficient a_coef(a);
  A.AddDomainIntegrator(new VectorMassIntegrator(a_coef));

  ConstantCoefficient lambda_coef(b);
  ConstantCoefficient mu_coef(b);
  A.AddDomainIntegrator(new ElasticityIntegrator(lambda_coef, mu_coef));
  A.Assemble(0);
  A.Finalize();

  std::unique_ptr<mfem::HypreParMatrix> J(A.ParallelAssemble());

  ParLinearForm             f(&fespace);
  VectorFunctionCoefficient load_func(dim, [&](const Vector& /*coords*/, Vector& force) {
    force    = 0.0;
    force(0) = -1.0;
  });

  f.AddDomainIntegrator(new VectorDomainLFIntegrator(load_func));
  f.Assemble();
  std::unique_ptr<mfem::HypreParVector> F(f.ParallelAssemble());

  ParGridFunction u_global(&fespace);
  u_global.Randomize();

  Vector U(fespace.TrueVSize());
  u_global.GetTrueDofs(U);

  [[maybe_unused]] static constexpr auto I = Identity<dim>();

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  residual.AddDomainIntegral(
      Dimension<dim>{},
      [&](auto /*x*/, auto displacement) {
        auto [u, du_dx] = displacement;
        auto body_force = a * u + I[0];
        auto strain = 0.5 * (du_dx + transpose(du_dx));
        auto stress = b * tr(strain) * I + 2.0 * b * strain;
        return std::tuple{body_force, stress};
      },
      mesh);

  mfem::Vector r1 = (*J) * U - (*F);
  mfem::Vector r2 = residual(U);

  if (verbose) {
    std::cout << "||r1||: " << r1.Norml2() << std::endl;
    std::cout << "||r2||: " << r2.Norml2() << std::endl;
    std::cout << "||r1-r2||/||r1||: " << mfem::Vector(r1 - r2).Norml2() / r1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), tol);

  mfem::Operator& grad = residual.GetGradient(U);

  mfem::Vector g1 = (*J) * U;
  mfem::Vector g2 = grad * U;

  if (verbose) {
    std::cout << "||g1||: " << g1.Norml2() << std::endl;
    std::cout << "||g2||: " << g2.Norml2() << std::endl;
    std::cout << "||g1-g2||/||g1||: " << mfem::Vector(g1 - g2).Norml2() / g1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), tol);
}

// this test sets up part of a toy "magnetic diffusion" problem where the residual includes contributions
// from a vector-potential-proportional J and an isotropically linear H
// 
// the same problem is expressed with mfem and weak_form, and their residuals and gradient action
// are compared to ensure the implementations are in agreement.
template <int p, int dim>
void weak_form_test(mfem::ParMesh& mesh, Hcurl<p> test, Hcurl<p> trial, Dimension<dim>)
{
  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = ND_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec);

  ParBilinearForm B(&fespace);

  ConstantCoefficient a_coef(a);
  B.AddDomainIntegrator(new VectorFEMassIntegrator(a_coef));

  ConstantCoefficient b_coef(b);
  B.AddDomainIntegrator(new CurlCurlIntegrator(b_coef));
  B.Assemble(0);
  B.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> J(B.ParallelAssemble());

  ParLinearForm             f(&fespace);
  VectorFunctionCoefficient load_func(dim, [&](const Vector& coords, Vector& output) {
    double x  = coords(0);
    double y  = coords(1);
    output    = 0.0;
    output(0) = 10 * x * y;
    output(1) = -5 * (x - y) * y;
  });

  f.AddDomainIntegrator(new VectorFEDomainLFIntegrator(load_func));
  f.Assemble();
  std::unique_ptr<mfem::HypreParVector> F(f.ParallelAssemble());

  ParGridFunction u_global(&fespace);
  u_global.Randomize();

  Vector U(fespace.TrueVSize());
  u_global.GetTrueDofs(U);

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  residual.AddDomainIntegral(
      Dimension<dim>{},
      [&](auto x, auto vector_potential) {
        auto [A, curl_A] = vector_potential;
        auto J = a * A - tensor<double, dim>{10 * x[0] * x[1], -5 * (x[0] - x[1]) * x[1]};
        auto H = b * curl_A;
        return std::tuple{J, H};
      },
      mesh);

  mfem::Vector r1 = (*J) * U - (*F);
  mfem::Vector r2 = residual(U);

  if (verbose) {
    std::cout << "||r1||: " << r1.Norml2() << std::endl;
    std::cout << "||r2||: " << r2.Norml2() << std::endl;
    std::cout << "||r1-r2||/||r1||: " << mfem::Vector(r1 - r2).Norml2() / r1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), tol);

  mfem::Operator& grad = residual.GetGradient(U);

  mfem::Vector g1 = (*J) * U;
  mfem::Vector g2 = grad * U;

  if (verbose) {
    std::cout << "||g1||: " << g1.Norml2() << std::endl;
    std::cout << "||g2||: " << g2.Norml2() << std::endl;
    std::cout << "||g1-g2||/||g1||: " << mfem::Vector(g1 - g2).Norml2() / g1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), tol);
}

const mfem::IntegrationRule& DiffusionIntegrator_GetRule(const mfem::FiniteElement& trial_fe,
                                                         const mfem::FiniteElement& test_fe)
{
  int order;
  if (trial_fe.Space() == mfem::FunctionSpace::Pk) {
    order = trial_fe.GetOrder() + test_fe.GetOrder() - 2;
  } else {
    // order = 2*el.GetOrder() - 2;  // <-- this seems to work fine too
    order = trial_fe.GetOrder() + test_fe.GetOrder() + trial_fe.GetDim() - 1;
  }

  if (trial_fe.Space() == mfem::FunctionSpace::rQk) {
    return RefinedIntRules.Get(trial_fe.GetGeomType(), order);
  }
  return mfem::IntRules.Get(trial_fe.GetGeomType(), order);
}

// Copy of DiffusionIntegrator
void DiffusionIntegrator_AssembleElementMatrix(mfem::Coefficient& Q, const mfem::FiniteElement& el,
                                               mfem::ElementTransformation& Trans, mfem::DenseMatrix& elmat)
{
  int    nd       = el.GetDof();
  int    dim      = el.GetDim();
  int    spaceDim = Trans.GetSpaceDim();
  bool   square   = (dim == spaceDim);
  double w;

  mfem::DenseMatrix dshape(nd, dim), dshapedxt(nd, spaceDim), invdfdx(dim, spaceDim);
  mfem::Vector      D(0);
  elmat.SetSize(nd);

  const mfem::IntegrationRule* ir = &DiffusionIntegrator_GetRule(el, el);

  elmat = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcDShape(ip, dshape);

    Trans.SetIntPoint(&ip);
    w = Trans.Weight();
    w = ip.weight / (square ? w : w * w * w);
    // AdjugateJacobian = / adj(J),         if J is square
    //                    \ adj(J^t.J).J^t, otherwise
    Mult(dshape, Trans.AdjugateJacobian(), dshapedxt);

    double q = Q.Eval(Trans, ip);

    w *= q;

    AddMult_a_AAt(w, dshapedxt, elmat);
  }
}

/* Start of dark arts to get at Bilinearform.mat*/
template <typename From, auto V, typename Result>
struct forbidden {
  friend Result _get_mat(From& from) { return from.*V; }
};

mfem::SparseMatrix* _get_mat(mfem::ParBilinearForm&);

template struct forbidden<mfem::ParBilinearForm, &mfem::ParBilinearForm::mat, mfem::SparseMatrix*>;

/* End of dark arts */

template <int p, int dim>
void weak_form_matrix_test(mfem::ParMesh& mesh, H1<p> test, H1<p> trial, Dimension<dim>)
{
  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = H1_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec);

  ParBilinearForm A(&fespace);

  ConstantCoefficient a_coef(a);
  A.AddDomainIntegrator(new MassIntegrator(a_coef));

  ConstantCoefficient   b_coef(b);
  [[maybe_unused]] auto diff_integ = new DiffusionIntegrator(b_coef);

  // modify integration rule temporarily to match weak_form?
  A.AddDomainIntegrator(diff_integ);
  constexpr int skip_zeros = 0;
  A.Assemble(skip_zeros);
  A.Finalize();

  // Save a deep-copy of rank local assembled sparse matrix
  mfem::SparseMatrix A_spmat_mfem(*_get_mat(A));

  std::unique_ptr<mfem::HypreParMatrix> J(A.ParallelAssemble());

  ParLinearForm       f(&fespace);
  FunctionCoefficient load_func([&](const Vector& coords) { return 100 * coords(0) * coords(1); });

  f.AddDomainIntegrator(new DomainLFIntegrator(load_func));
  f.Assemble();
  std::unique_ptr<mfem::HypreParVector> F(f.ParallelAssemble());

  ParGridFunction u_global(&fespace);
  u_global.Randomize();

  Vector U(fespace.TrueVSize());
  u_global.GetTrueDofs(U);

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  residual.AddDomainIntegral(
      Dimension<dim>{},
      [&]([[maybe_unused]] auto x, auto temperature) {
        auto [u, du_dx] = temperature;
        auto f0         = a * u - (100 * x[0] * x[1]);
        auto f1         = b * du_dx;
        return std::tuple{f0, f1};
      },
      mesh);

  mfem::Vector r1 = (*J) * U - (*F);
  mfem::Vector r2 = residual(U);

  if (verbose) {
    std::cout << "||r1||: " << r1.Norml2() << std::endl;
    std::cout << "||r2||: " << r2.Norml2() << std::endl;
    std::cout << "||r1-r2||/||r1||: " << mfem::Vector(r1 - r2).Norml2() / r1.Norml2() << std::endl;
  }
  //  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), tol);

  mfem::Operator& grad2 = residual.GetGradient(U);

  mfem::Vector g1 = (*J) * U;
  mfem::Vector g2 = grad2 * U;

  if (verbose) {
    std::cout << "||g1||: " << g1.Norml2() << std::endl;
    std::cout << "||g2||: " << g2.Norml2() << std::endl;
    std::cout << "||g1-g2||/||g1||: " << mfem::Vector(g1 - g2).Norml2() / g1.Norml2() << std::endl;
  }
  // EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), tol);

  Array<int> dofs;
  fespace.GetElementDofs(0, dofs);
  mfem::Vector K_e(mesh.GetNE() * dofs.Size() * fespace.GetVDim() * dofs.Size() * fespace.GetVDim());
  K_e = 0.;
  residual.GradientMatrix(K_e);
  std::cout << "K_e: (" << K_e.Size() << ")" << std::endl << std::endl;
  ;

  // // Reprocess each element from LEXICOGRAPHIC -> NATIVE for the tensorbasiscase
  // auto inv_dof_map = dynamic_cast<const mfem::TensorBasisElement*>(fespace.GetFE(0))
  //                        ->GetDofMap();  // for quads the change from NATIVE -> lexicographic is the same

  // mfem::SparseMatrix A_spmat_weak(A_spmat_mfem.Height());
  // constexpr auto     ordering_type = mfem::Ordering::byNODES;
  // {
  //   auto dk =
  //       mfem::Reshape(K_e.ReadWrite(), dofs.Size() * fespace.GetVDim(), dofs.Size() * fespace.GetVDim(),
  //       mesh.GetNE());
  //   for (int e = 0; e < mesh.GetNE(); e++) {
  //     DenseMatrix mat(dofs.Size() * fespace.GetVDim());
  //     for (int i = 0; i < dofs.Size(); i++) {
  //       for (int id = 0; id < fespace.GetVDim(); id++) {
  //         for (int j = 0; j < dofs.Size(); j++) {
  //           for (int jd = 0; jd < fespace.GetVDim(); jd++) {
  //             int inv_i_vdof = mfem::Ordering::Map<ordering_type>(dofs.Size(), fespace.GetVDim(), inv_dof_map[i],
  //             id); int inv_j_vdof = mfem::Ordering::Map<ordering_type>(dofs.Size(), fespace.GetVDim(),
  //             inv_dof_map[j], jd); int i_vdof     = mfem::Ordering::Map<ordering_type>(dofs.Size(),
  //             fespace.GetVDim(), i, id); int j_vdof     = mfem::Ordering::Map<ordering_type>(dofs.Size(),
  //             fespace.GetVDim(), j, id); mat(inv_i_vdof, inv_j_vdof) = dk(i_vdof, j_vdof, e);
  //           }
  //         }
  //       }
  //     }
  //     // Copy back
  //     for (int i = 0; i < dofs.Size(); i++) {
  //       for (int id = 0; id < fespace.GetVDim(); id++) {
  //         for (int j = 0; j < dofs.Size(); j++) {
  //           for (int jd = 0; jd < fespace.GetVDim(); jd++) {
  //             int i_vdof            = mfem::Ordering::Map<ordering_type>(dofs.Size(), fespace.GetVDim(), i, id);
  //             int j_vdof            = mfem::Ordering::Map<ordering_type>(dofs.Size(), fespace.GetVDim(), j, id);
  //             dk(i_vdof, j_vdof, e) = mat(i_vdof, j_vdof);
  //           }
  //         }
  //       }
  //     }

  //     Array<int> elem_vdofs;
  //     fespace.GetElementVDofs(e, elem_vdofs);
  //     A_spmat_weak.AddSubMatrix(elem_vdofs, elem_vdofs, mat, skip_zeros);
  //   }
  // }
  // A_spmat_weak.Finalize();

  mfem::SparseMatrix       A_spmat_weak(A_spmat_mfem.Height());
  mfem::ElementRestriction elemRestriction(fespace, mfem::ElementDofOrdering::LEXICOGRAPHIC);
  elemRestriction.FillSparseMatrix(K_e, A_spmat_weak);
  A_spmat_weak.Finalize();

  for (int r = 0; r < A_spmat_weak.Height(); r++) {
    auto columns = A_spmat_weak.GetRowColumns(r);
    for (int c = 0; c < A_spmat_weak.RowSize(r); c++) {
      EXPECT_NEAR(A_spmat_mfem(r, columns[c]), A_spmat_weak(r, columns[c]), 1.e-10);
    }
  }

  // test AssembledSparseMatrix
  serac::mfem_ext::AssembledSparseMatrix A_serac_mat(fespace, fespace, mfem::ElementDofOrdering::LEXICOGRAPHIC);
  A_serac_mat.FillData(K_e);
  A_serac_mat.Finalize();

  for (int r = 0; r < A_serac_mat.Height(); r++) {
    auto columns = A_serac_mat.GetRowColumns(r);
    std::cout << "row " << r << " : " << A_spmat_mfem.RowSize(r) << std::endl;
    for (int c = 0; c < A_serac_mat.RowSize(r); c++) {
      EXPECT_NEAR(A_spmat_mfem(r, columns[c]), A_serac_mat(r, columns[c]), 1.e-10);
    }
  }
}

template <int p, int dim>
void weak_form_matrix_test(mfem::ParMesh& mesh, H1<p, dim> test, H1<p, dim> trial, Dimension<dim>)
{
  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = H1_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec, dim);

  ParBilinearForm A(&fespace);

  ConstantCoefficient a_coef(a);
  auto                VMI = new VectorMassIntegrator(a_coef);
  A.AddDomainIntegrator(VMI);

  ConstantCoefficient lambda_coef(b);
  ConstantCoefficient mu_coef(b);
  auto                EI = new ElasticityIntegrator(lambda_coef, mu_coef);
  A.AddDomainIntegrator(EI);
  constexpr int skip_zeros = 0;
  A.Assemble(skip_zeros);
  A.Finalize();

  // Save a deep-copy of rank local assembled sparse matrix
  mfem::SparseMatrix                    A_spmat_mfem(*_get_mat(A));
  std::unique_ptr<mfem::HypreParMatrix> J(A.ParallelAssemble());

  ParLinearForm             f(&fespace);
  VectorFunctionCoefficient load_func(dim, [&](const Vector& /*coords*/, Vector& force) {
    force    = 0.0;
    force(0) = -1.0;
  });

  f.AddDomainIntegrator(new VectorDomainLFIntegrator(load_func));
  f.Assemble();
  std::unique_ptr<mfem::HypreParVector> F(f.ParallelAssemble());

  ParGridFunction u_global(&fespace);
  u_global.Randomize();

  Vector U(fespace.TrueVSize());
  u_global.GetTrueDofs(U);

  [[maybe_unused]] static constexpr auto I = Identity<dim>();

  using test_space  = decltype(test);
  using trial_space = decltype(trial);

  WeakForm<test_space(trial_space)> residual(&fespace, &fespace);

  residual.AddDomainIntegral(
      Dimension<dim>{},
      [&](auto /*x*/, auto displacement) {
        auto [u, du_dx] = displacement;
        auto f0         = a * u + I[0];
        auto strain     = 0.5 * (du_dx + transpose(du_dx));
        auto f1         = b * tr(strain) * I + 2.0 * b * strain;
        return std::tuple{f0, f1};
      },
      mesh);

  mfem::Vector r1 = (*J) * U - (*F);
  mfem::Vector r2 = residual(U);

  if (verbose) {
    std::cout << "||r1||: " << r1.Norml2() << std::endl;
    std::cout << "||r2||: " << r2.Norml2() << std::endl;
    std::cout << "||r1-r2||/||r1||: " << mfem::Vector(r1 - r2).Norml2() / r1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), tol);

  mfem::Operator& grad = residual.GetGradient(U);

  mfem::Vector g1 = (*J) * U;
  mfem::Vector g2 = grad * U;

  if (verbose) {
    std::cout << "||g1||: " << g1.Norml2() << std::endl;
    std::cout << "||g2||: " << g2.Norml2() << std::endl;
    std::cout << "||g1-g2||/||g1||: " << mfem::Vector(g1 - g2).Norml2() / g1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), tol);

  Array<int> dofs;
  fespace.GetElementDofs(0, dofs);
  mfem::Vector K_e(mesh.GetNE() * dofs.Size() * fespace.GetVDim() * dofs.Size() * fespace.GetVDim());
  K_e = 0.;
  residual.GradientMatrix(K_e);
  std::cout << "K_e: (" << K_e.Size() << ")" << std::endl << std::endl;

  serac::mfem_ext::AssembledSparseMatrix A_serac_mat(fespace, fespace, mfem::ElementDofOrdering::LEXICOGRAPHIC);
  A_serac_mat.FillData(K_e);
  A_serac_mat.Finalize();

  for (int r = 0; r < A_spmat_mfem.Height(); r++) {
    auto columns = A_spmat_mfem.GetRowColumns(r);
    std::cout << "row " << r << " : " << A_spmat_mfem.RowSize(r) << std::endl;
    for (int c = 0; c < A_spmat_mfem.RowSize(r); c++) {
      EXPECT_NEAR(A_spmat_mfem(r, columns[c]), A_serac_mat(r, columns[c]), 1.e-10);
    }
  }

  A_serac_mat.Print();

  A_spmat_mfem.Print();
}

template <int p, int dim>
void weak_form_matrix_test(mfem::ParMesh& mesh, Hcurl<p> test, Hcurl<p> trial, Dimension<dim>)
{
  static constexpr double a = 1.7;
  static constexpr double b = 2.1;

  auto                  fec = ND_FECollection(p, dim);
  ParFiniteElementSpace fespace(&mesh, &fec);

  ParBilinearForm B(&fespace);

  ConstantCoefficient a_coef(a);
  B.AddDomainIntegrator(new VectorFEMassIntegrator(a_coef));

  ConstantCoefficient b_coef(b);
  B.AddDomainIntegrator(new CurlCurlIntegrator(b_coef));
  B.Assemble(0);
  B.Finalize();

  mfem::SparseMatrix                    B_spmat_mfem(*_get_mat(B));
  std::unique_ptr<mfem::HypreParMatrix> J(B.ParallelAssemble());

  ParLinearForm             f(&fespace);
  VectorFunctionCoefficient load_func(dim, [&](const Vector& coords, Vector& output) {
    double x  = coords(0);
    double y  = coords(1);
    output    = 0.0;
    output(0) = 10 * x * y;
    output(1) = -5 * (x - y) * y;
  });

  f.AddDomainIntegrator(new VectorFEDomainLFIntegrator(load_func));
  f.Assemble();
  std::unique_ptr<mfem::HypreParVector> F(f.ParallelAssemble());

  ParGridFunction u_global(&fespace);
  u_global.Randomize();

  Vector U(fespace.TrueVSize());
  u_global.GetTrueDofs(U);

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

  mfem::Vector r1 = (*J) * U - (*F);
  mfem::Vector r2 = residual(U);

  if (verbose) {
    std::cout << "||r1||: " << r1.Norml2() << std::endl;
    std::cout << "||r2||: " << r2.Norml2() << std::endl;
    std::cout << "||r1-r2||/||r1||: " << mfem::Vector(r1 - r2).Norml2() / r1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(r1 - r2).Norml2() / r1.Norml2(), tol);

  mfem::Operator& grad = residual.GetGradient(U);

  mfem::Vector g1 = (*J) * U;
  mfem::Vector g2 = grad * U;

  if (verbose) {
    std::cout << "||g1||: " << g1.Norml2() << std::endl;
    std::cout << "||g2||: " << g2.Norml2() << std::endl;
    std::cout << "||g1-g2||/||g1||: " << mfem::Vector(g1 - g2).Norml2() / g1.Norml2() << std::endl;
  }
  EXPECT_NEAR(0., mfem::Vector(g1 - g2).Norml2() / g1.Norml2(), tol);

  Array<int> dofs;
  fespace.GetElementDofs(0, dofs);
  mfem::Vector K_e(mesh.GetNE() * dofs.Size() * fespace.GetVDim() * dofs.Size() * fespace.GetVDim());
  K_e = 0.;
  residual.GradientMatrix(K_e);
  std::cout << "K_e: (" << K_e.Size() << ")" << std::endl << std::endl;

  serac::mfem_ext::AssembledSparseMatrix B_serac_mat(fespace, fespace, mfem::ElementDofOrdering::LEXICOGRAPHIC);
  B_serac_mat.FillData(K_e);
  B_serac_mat.Finalize();

  for (int r = 0; r < B_spmat_mfem.Height(); r++) {
    auto columns = B_spmat_mfem.GetRowColumns(r);
    std::cout << "row " << r << " : " << B_spmat_mfem.RowSize(r) << std::endl;
    for (int c = 0; c < B_spmat_mfem.RowSize(r); c++) {
      EXPECT_NEAR(B_spmat_mfem(r, columns[c]), B_serac_mat(r, columns[c]), 1.e-10);
    }
  }

  B_serac_mat.Print();

  B_spmat_mfem.Print();
}

// TEST(thermal, 2D_linear) { weak_form_test(*mesh2D, H1<1>{}, H1<1>{}, Dimension<2>{}); }
// TEST(thermal, 2D_quadratic) { weak_form_test(*mesh2D, H1<2>{}, H1<2>{}, Dimension<2>{}); }
// TEST(thermal, 2D_cubic) { weak_form_test(*mesh2D, H1<3>{}, H1<3>{}, Dimension<2>{}); }

TEST(thermal, 2D_linear_mat) { weak_form_matrix_test(*mesh2D, H1<1>{}, H1<1>{}, Dimension<2>{}); }
TEST(thermal, 2D_quadratic_mat) { weak_form_matrix_test(*mesh2D, H1<2>{}, H1<2>{}, Dimension<2>{}); }
// TEST(thermal, 2D_cubic_mat) { weak_form_test(*mesh2D, H1<3>{}, H1<3>{}, Dimension<2>{}); }

// TEST(thermal, 3D_linear) { weak_form_test(*mesh3D, H1<1>{}, H1<1>{}, Dimension<3>{}); }
// TEST(thermal, 3D_quadratic) { weak_form_test(*mesh3D, H1<2>{}, H1<2>{}, Dimension<3>{}); }
// TEST(thermal, 3D_cubic) { weak_form_test(*mesh3D, H1<3>{}, H1<3>{}, Dimension<3>{}); }

TEST(thermal, 3D_linear_mat) { weak_form_matrix_test(*mesh3D, H1<1>{}, H1<1>{}, Dimension<3>{}); }
// TEST(thermal, 3D_quadratic_mat) { weak_form_matrix_test(*mesh3D, H1<2>{}, H1<2>{}, Dimension<3>{}); }
// TEST(thermal, 3D_cubic_mat) { weak_form_matrix_test(*mesh3D, H1<3>{}, H1<3>{}, Dimension<3>{}); }

// TEST(hcurl, 2D_linear) { weak_form_test(*mesh2D, Hcurl<1>{}, Hcurl<1>{}, Dimension<2>{}); }
// TEST(hcurl, 2D_quadratic) { weak_form_test(*mesh2D, Hcurl<2>{}, Hcurl<2>{}, Dimension<2>{}); }
// TEST(hcurl, 2D_cubic) { weak_form_test(*mesh2D, Hcurl<3>{}, Hcurl<3>{}, Dimension<2>{}); }

// TEST(hcurl, 2D_linear_mat) { weak_form_matrix_test(*mesh2D, Hcurl<1>{}, Hcurl<1>{}, Dimension<2>{}); }

// TEST(hcurl, 3D_linear) { weak_form_test(*mesh3D, Hcurl<1>{}, Hcurl<1>{}, Dimension<3>{}); }
// TEST(hcurl, 3D_quadratic) { weak_form_test(*mesh3D, Hcurl<2>{}, Hcurl<2>{}, Dimension<3>{}); }
// TEST(hcurl, 3D_cubic) { weak_form_test(*mesh3D, Hcurl<3>{}, Hcurl<3>{}, Dimension<3>{}); }

// TEST(elasticity, 2D_linear) { weak_form_test(*mesh2D, H1<1, 2>{}, H1<1, 2>{}, Dimension<2>{}); }
// TEST(elasticity, 2D_quadratic) { weak_form_test(*mesh2D, H1<2, 2>{}, H1<2, 2>{}, Dimension<2>{}); }
// TEST(elasticity, 2D_cubic) { weak_form_test(*mesh2D, H1<3, 2>{}, H1<3, 2>{}, Dimension<2>{}); }

TEST(elasticity, 2D_linear_mat) { weak_form_matrix_test(*mesh2D, H1<1, 2>{}, H1<1, 2>{}, Dimension<2>{}); }
TEST(elasticity, 2D_quadratic_mat) { weak_form_matrix_test(*mesh2D, H1<2, 2>{}, H1<2, 2>{}, Dimension<2>{}); }

// TEST(elasticity, 3D_linear) { weak_form_test(*mesh3D, H1<1, 3>{}, H1<1, 3>{}, Dimension<3>{}); }
// TEST(elasticity, 3D_quadratic) { weak_form_test(*mesh3D, H1<2, 3>{}, H1<2, 3>{}, Dimension<3>{}); }
// TEST(elasticity, 3D_cubic) { weak_form_test(*mesh3D, H1<3, 3>{}, H1<3, 3>{}, Dimension<3>{}); }

TEST(elasticity, 3D_linear_mat) { weak_form_matrix_test(*mesh3D, H1<1, 3>{}, H1<1, 3>{}, Dimension<3>{}); }
// TEST(elasticity, 3D_quadratic_mat) { weak_form_matrix_test(*mesh3D, H1<2, 3>{}, H1<2, 3>{}, Dimension<3>{}); }

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  axom::slic::SimpleLogger logger;

  int serial_refinement   = 0;
  int parallel_refinement = 0;

  // std::string meshfile2D = SERAC_REPO_DIR "/data/meshes/star.mesh";
  // mesh2D = mesh::refineAndDistribute(buildMeshFromFile(meshfile2D), serial_refinement, parallel_refinement);
  mesh2D = mesh::refineAndDistribute(serac::buildRectangleMesh(2, 1), serial_refinement, parallel_refinement);

  std::string meshfile3D = SERAC_REPO_DIR "/data/meshes/beam-hex.mesh";
  mesh3D = mesh::refineAndDistribute(buildMeshFromFile(meshfile3D), serial_refinement, parallel_refinement);

  int result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
