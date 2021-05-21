// specialization of finite_element for L2 on segment geometry
//
// this specialization defines shape functions (and their gradients) that
// interpolate at Gauss-Lobatto nodes for the appropriate polynomial order
// 
// note: mfem assumes the parent element domain is [0,1]
template <int p, int c>
struct finite_element<Geometry::Segment, L2< p, c > > {

  static constexpr auto geometry = Geometry::Segment;
  static constexpr auto family = Family::L2;
  static constexpr int components = c;
  static constexpr int dim = 1;
  static constexpr int ndof = (p + 1);

  static constexpr tensor<double, ndof> shape_functions(double xi) {
    return GaussLegendreInterpolation<ndof>(xi);
  }

  static constexpr tensor<double, ndof> shape_function_gradients(double xi) {
    return GaussLegendreInterpolationDerivative<ndof>(xi);
  }
};
