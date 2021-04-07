#include <chrono>
#include <iostream>

#include "mfem.hpp"

#include "serac/serac_config.hpp"
#include "serac/numerics/expr_template_ops.hpp"

#include "serac/physics/utilities/variational_form/tensor.hpp"
#include "serac/physics/utilities/variational_form/integral.hpp"
#include "serac/physics/utilities/variational_form/quadrature.hpp"
#include "serac/physics/utilities/variational_form/finite_element.hpp"
#include "serac/physics/utilities/variational_form/tuple_arithmetic.hpp"

class timer {
  typedef std::chrono::high_resolution_clock::time_point time_point;
  typedef std::chrono::duration<double>                  duration_type;

public:
  void start() { then = std::chrono::high_resolution_clock::now(); }

  void stop() { now = std::chrono::high_resolution_clock::now(); }

  double elapsed() { return std::chrono::duration_cast<duration_type>(now - then).count(); }

private:
  time_point then, now;
};

namespace mfem {

template<int T_D1D = 0, int T_Q1D = 0>
static void PADiffusionApply3D(const int NE,
                               const bool symmetric,
                               const Array<double> &b,
                               const Array<double> &g,
                               const Array<double> &bt,
                               const Array<double> &gt,
                               const Vector &d_,
                               const Vector &x_,
                               Vector &y_,
                               int d1d = 0, int q1d = 0)
{
   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;
   MFEM_VERIFY(D1D <= MAX_D1D, "");
   MFEM_VERIFY(Q1D <= MAX_Q1D, "");
   auto B = Reshape(b.Read(), Q1D, D1D);
   auto G = Reshape(g.Read(), Q1D, D1D);
   auto Bt = Reshape(bt.Read(), D1D, Q1D);
   auto Gt = Reshape(gt.Read(), D1D, Q1D);
   auto D = Reshape(d_.Read(), Q1D*Q1D*Q1D, symmetric ? 6 : 9, NE);
   auto X = Reshape(x_.Read(), D1D, D1D, D1D, NE);
   auto Y = Reshape(y_.ReadWrite(), D1D, D1D, D1D, NE);
   for (int e = 0; e < NE; e++) {
      constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
      constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
      double grad[max_Q1D][max_Q1D][max_Q1D][3];
      for (int qz = 0; qz < Q1D; ++qz)
      {
         for (int qy = 0; qy < Q1D; ++qy)
         {
            for (int qx = 0; qx < Q1D; ++qx)
            {
               grad[qz][qy][qx][0] = 0.0;
               grad[qz][qy][qx][1] = 0.0;
               grad[qz][qy][qx][2] = 0.0;
            }
         }
      }
      for (int dz = 0; dz < D1D; ++dz)
      {
         double gradXY[max_Q1D][max_Q1D][3];
         for (int qy = 0; qy < Q1D; ++qy)
         {
            for (int qx = 0; qx < Q1D; ++qx)
            {
               gradXY[qy][qx][0] = 0.0;
               gradXY[qy][qx][1] = 0.0;
               gradXY[qy][qx][2] = 0.0;
            }
         }
         for (int dy = 0; dy < D1D; ++dy)
         {
            double gradX[max_Q1D][2];
            for (int qx = 0; qx < Q1D; ++qx)
            {
               gradX[qx][0] = 0.0;
               gradX[qx][1] = 0.0;
            }
            for (int dx = 0; dx < D1D; ++dx)
            {
               const double s = X(dx,dy,dz,e);
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  gradX[qx][0] += s * B(qx,dx);
                  gradX[qx][1] += s * G(qx,dx);
               }
            }
            for (int qy = 0; qy < Q1D; ++qy)
            {
               const double wy  = B(qy,dy);
               const double wDy = G(qy,dy);
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  const double wx  = gradX[qx][0];
                  const double wDx = gradX[qx][1];
                  gradXY[qy][qx][0] += wDx * wy;
                  gradXY[qy][qx][1] += wx  * wDy;
                  gradXY[qy][qx][2] += wx  * wy;
               }
            }
         }
         for (int qz = 0; qz < Q1D; ++qz)
         {
            const double wz  = B(qz,dz);
            const double wDz = G(qz,dz);
            for (int qy = 0; qy < Q1D; ++qy)
            {
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  grad[qz][qy][qx][0] += gradXY[qy][qx][0] * wz;
                  grad[qz][qy][qx][1] += gradXY[qy][qx][1] * wz;
                  grad[qz][qy][qx][2] += gradXY[qy][qx][2] * wDz;
               }
            }
         }
      }
      // Calculate Dxyz, xDyz, xyDz in plane
      for (int qz = 0; qz < Q1D; ++qz)
      {
         for (int qy = 0; qy < Q1D; ++qy)
         {
            for (int qx = 0; qx < Q1D; ++qx)
            {
               const int q = qx + (qy + qz * Q1D) * Q1D;
               const double O11 = D(q,0,e);
               const double O12 = D(q,1,e);
               const double O13 = D(q,2,e);
               const double O21 = symmetric ? O12 : D(q,3,e);
               const double O22 = symmetric ? D(q,3,e) : D(q,4,e);
               const double O23 = symmetric ? D(q,4,e) : D(q,5,e);
               const double O31 = symmetric ? O13 : D(q,6,e);
               const double O32 = symmetric ? O23 : D(q,7,e);
               const double O33 = symmetric ? D(q,5,e) : D(q,8,e);
               const double gradX = grad[qz][qy][qx][0];
               const double gradY = grad[qz][qy][qx][1];
               const double gradZ = grad[qz][qy][qx][2];
               grad[qz][qy][qx][0] = (O11*gradX)+(O12*gradY)+(O13*gradZ);
               grad[qz][qy][qx][1] = (O21*gradX)+(O22*gradY)+(O23*gradZ);
               grad[qz][qy][qx][2] = (O31*gradX)+(O32*gradY)+(O33*gradZ);
            }
         }
      }
      for (int qz = 0; qz < Q1D; ++qz)
      {
         double gradXY[max_D1D][max_D1D][3];
         for (int dy = 0; dy < D1D; ++dy)
         {
            for (int dx = 0; dx < D1D; ++dx)
            {
               gradXY[dy][dx][0] = 0;
               gradXY[dy][dx][1] = 0;
               gradXY[dy][dx][2] = 0;
            }
         }
         for (int qy = 0; qy < Q1D; ++qy)
         {
            double gradX[max_D1D][3];
            for (int dx = 0; dx < D1D; ++dx)
            {
               gradX[dx][0] = 0;
               gradX[dx][1] = 0;
               gradX[dx][2] = 0;
            }
            for (int qx = 0; qx < Q1D; ++qx)
            {
               const double gX = grad[qz][qy][qx][0];
               const double gY = grad[qz][qy][qx][1];
               const double gZ = grad[qz][qy][qx][2];
               for (int dx = 0; dx < D1D; ++dx)
               {
                  const double wx  = Bt(dx,qx);
                  const double wDx = Gt(dx,qx);
                  gradX[dx][0] += gX * wDx;
                  gradX[dx][1] += gY * wx;
                  gradX[dx][2] += gZ * wx;
               }
            }
            for (int dy = 0; dy < D1D; ++dy)
            {
               const double wy  = Bt(dy,qy);
               const double wDy = Gt(dy,qy);
               for (int dx = 0; dx < D1D; ++dx)
               {
                  gradXY[dy][dx][0] += gradX[dx][0] * wy;
                  gradXY[dy][dx][1] += gradX[dx][1] * wDy;
                  gradXY[dy][dx][2] += gradX[dx][2] * wy;
               }
            }
         }
         for (int dz = 0; dz < D1D; ++dz)
         {
            const double wz  = Bt(dz,qz);
            const double wDz = Gt(dz,qz);
            for (int dy = 0; dy < D1D; ++dy)
            {
               for (int dx = 0; dx < D1D; ++dx)
               {
                  Y(dx,dy,dz,e) +=
                     ((gradXY[dy][dx][0] * wz) +
                      (gradXY[dy][dx][1] * wz) +
                      (gradXY[dy][dx][2] * wDz));
               }
            }
         }
      }
   }
}

}

template < ::Geometry g, int P, int Q, typename lambda > 
void H1_kernel(const mfem::Vector & U, mfem::Vector & R, const mfem::Vector & J_, int num_elements, lambda && qf) {

  using trial = H1<P>;
  using test = H1<P>;
  using test_element = finite_element< g, trial >;
  using trial_element = finite_element< g, test >;
  using element_residual_type = typename trial_element::residual_type;
  static constexpr int dim = dimension_of(g);
  static constexpr int test_ndof = test_element::ndof;
  static constexpr int trial_ndof = trial_element::ndof;
  static constexpr auto rule = GaussQuadratureRule< g, Q >();

  auto J = mfem::Reshape(J_.Read(), rule.size(), dim, dim, num_elements);
  auto u = impl::Reshape<trial>(U.Read(), trial_ndof, num_elements);
  auto r = impl::Reshape<test>(R.ReadWrite(), test_ndof, num_elements);

  for (int e = 0; e < num_elements; e++) {
    tensor u_elem = impl::Load<trial_element>(u, e);

    element_residual_type r_elem{};

    for (int q = 0; q < static_cast<int>(rule.size()); q++) {
      auto xi = rule.points[q];
      auto dxi = rule.weights[q];
      auto J_q = make_tensor< dim, dim >([&](int i, int j){ return J(q, i, j, e); });
      double dx = impl::Measure(J_q) * dxi;

      auto dN = trial_element::shape_function_gradients(xi);
      auto inv_J = inv(J_q);

      auto grad_u = dot(dot(u_elem, dN), inv_J);

      auto qf_output = qf(grad_u);

      r_elem += dot(dN, dot(inv_J, qf_output)) * dx;
    }

    impl::Add(r, r_elem, e);
  }

}

template < ::Geometry g, int P, int Q, typename lambda > 
void H1_kernel_constexpr(const mfem::Vector & U, mfem::Vector & R, const mfem::Vector & J_, int num_elements, lambda && qf) {

  using trial = H1<P>;
  using test = H1<P>;
  using test_element = finite_element< g, trial >;
  using trial_element = finite_element< g, test >;
  using element_residual_type = typename trial_element::residual_type;
  static constexpr int dim = dimension_of(g);
  static constexpr int test_ndof = test_element::ndof;
  static constexpr int trial_ndof = trial_element::ndof;
  static constexpr auto rule = GaussQuadratureRule< g, Q >();

  auto J = mfem::Reshape(J_.Read(), rule.size(), dim, dim, num_elements);
  auto u = impl::Reshape<trial>(U.Read(), trial_ndof, num_elements);
  auto r = impl::Reshape<test>(R.ReadWrite(), test_ndof, num_elements);

  for (int e = 0; e < num_elements; e++) {
    tensor u_elem = impl::Load<trial_element>(u, e);

    element_residual_type r_elem{};

    for_constexpr< rule.size() >([&](auto q){
      static constexpr auto xi = rule.points[q];
      static constexpr auto dxi = rule.weights[q];
      auto J_q = make_tensor< dim, dim >([&](int i, int j){ return J(q, i, j, e); });
      double dx = impl::Measure(J_q) * dxi;

      static constexpr auto dN = trial_element::shape_function_gradients(xi);
      auto inv_J = inv(J_q);

      auto grad_u = dot(dot(u_elem, dN), inv_J);

      auto qf_output = qf(grad_u);

      r_elem += dot(dN, dot(inv_J, qf_output)) * dx;
    });

    impl::Add(r, r_elem, e);
  }

}


int main(int argc, char* argv[]) {

  int p = 1;
  int refinements = 0;
  const char * mesh_file = SERAC_REPO_DIR"/data/meshes/beam-hex.mesh";

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p, "-p", "--polynomial_order", "");
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&refinements, "-r", "--ref", "");

  args.Parse();
  if (!args.Good()) {
    exit(1);
  }
  args.PrintOptions(std::cout);

  mfem::Mesh mesh(mesh_file, 1, 1);
  for (int l = 0; l < refinements; l++) {
    mesh.UniformRefinement();
  }

  auto fec = mfem::H1_FECollection(p, 3);
  mfem::FiniteElementSpace fespace(&mesh, &fec);

  auto num_elements = mesh.GetNE();

  std::cout << "mesh contains " << num_elements << " elements" << std::endl;

  const mfem::FiniteElement & el = *(fespace.GetFE(0));
  const mfem::IntegrationRule & ir = mfem::IntRules.Get(el.GetGeomType(), el.GetOrder() * 2);
  auto geom = mesh.GetGeometricFactors(ir, mfem::GeometricFactors::JACOBIANS);

  auto maps = &el.GetDofToQuad(ir, mfem::DofToQuad::TENSOR);

  auto G = fespace.GetElementRestriction(mfem::ElementDofOrdering::LEXICOGRAPHIC);

  mfem::Vector U_L(fespace.GetTrueVSize(), mfem::Device::GetMemoryType());
  mfem::Vector U_E(G->Height(), mfem::Device::GetMemoryType());

  mfem::Vector J_Q = geom->J;

  mfem::Vector R_E(G->Height(), mfem::Device::GetMemoryType());
  mfem::Vector R_E2(G->Height(), mfem::Device::GetMemoryType());
  mfem::Vector R_E3(G->Height(), mfem::Device::GetMemoryType());

  mfem::Vector R_L(fespace.GetTrueVSize(), mfem::Device::GetMemoryType());

  U_L.Randomize();

  timer stopwatch[3];

  stopwatch[0].start();
  G->Mult(U_L, U_E); 
  stopwatch[0].stop();
  std::cout << "U_L -> U_E time: " << stopwatch[0].elapsed() << std::endl;

  static constexpr double k = 1.0;
  constexpr auto diffusion_qfunc = [](auto grad_u){
    return k * grad_u; // heat_flux
  };

  bool symmetric = false;
  if (p == 1) {
    stopwatch[0].start();
    mfem::PADiffusionApply3D< 2, 2 >(num_elements, symmetric, maps->B, maps->G, maps->Bt, maps->Gt, J_Q, U_E, R_E);
    stopwatch[0].stop();

    stopwatch[1].start();
    H1_kernel<Geometry::Hexahedron, 1, 2 >(U_E, R_E2, J_Q, num_elements, diffusion_qfunc);
    stopwatch[1].stop();

    stopwatch[2].start();
    H1_kernel_constexpr<Geometry::Hexahedron, 1, 2 >(U_E, R_E3, J_Q, num_elements, diffusion_qfunc);
    stopwatch[2].stop();

  }

  if (p == 2) {
    stopwatch[0].start();
    mfem::PADiffusionApply3D< 3, 3 >(num_elements, symmetric, maps->B, maps->G, maps->Bt, maps->Gt, J_Q, U_E, R_E);
    stopwatch[0].stop();
    stopwatch[1].start();
    H1_kernel<Geometry::Hexahedron, 2, 3 >(U_E, R_E2, J_Q, num_elements, diffusion_qfunc);
    stopwatch[1].stop();
    stopwatch[2].start();
    H1_kernel_constexpr<Geometry::Hexahedron, 2, 3 >(U_E, R_E3, J_Q, num_elements, diffusion_qfunc);
    stopwatch[2].stop();
  }

  if (p == 3) {
    stopwatch[0].start();
    mfem::PADiffusionApply3D< 4, 4 >(num_elements, symmetric, maps->B, maps->G, maps->Bt, maps->Gt, J_Q, U_E, R_E);
    stopwatch[0].stop();
    stopwatch[1].start();
    H1_kernel<Geometry::Hexahedron, 3, 4 >(U_E, R_E2, J_Q, num_elements, diffusion_qfunc);
    stopwatch[1].stop();
    stopwatch[2].start();
    H1_kernel_constexpr<Geometry::Hexahedron, 3, 4 >(U_E, R_E3, J_Q, num_elements, diffusion_qfunc);
    stopwatch[2].stop();
  }

  std::cout << "mfem::PADiffusionApply3D time: " << stopwatch[0].elapsed() << " seconds" << std::endl;
  std::cout << "H1_kernel< diffusion > time: " << stopwatch[1].elapsed() << " seconds" << std::endl;
  std::cout << "H1_kernel_constexpr< diffusion > time: " << stopwatch[2].elapsed() << " seconds" << std::endl;
  std::cout << "||R_1||: " << R_E.Norml2() << std::endl;
  std::cout << "||R_2||: " << R_E2.Norml2() << std::endl;
  std::cout << "||R_3||: " << R_E3.Norml2() << std::endl;

}