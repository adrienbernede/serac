#include "tuple_arithmetic.hpp"

#include <iostream>

int main() {

  constexpr auto abs = [](auto x){ return (x < 0) ? -x : x; };

  static constexpr auto I = Identity<3>();
  static constexpr double rho = 3.0;
  static constexpr double mu = 2.0;
  constexpr auto sigma = [](auto p, auto v, auto L) {      
    return rho * outer(v, v) + 2.0 * mu * sym(L) - p * I;
  };

  constexpr auto dsigma_dp = [](auto /*p*/, auto /*v*/, auto /*L*/) { return -1.0 * I; };

  constexpr auto dsigma_dv = [&](auto /*p*/, auto v, auto /*L*/) { 
    return make_tensor<3,3,3>([&](int i, int j, int k){
      return rho * ((i == k) * v[j] + (j == k) * v[i]);
    }); 
  };

  constexpr auto dsigma_dL = [&](auto /*p*/, auto /*v*/, auto /*L*/) { 
    return make_tensor<3,3,3,3>([&](int i, int j, int k, int l){
      return mu * ((i == k) * (j == l) + (i == l) * (j == k));
    });
  };

  constexpr double p = 3.14;
  constexpr tensor v = {{1.0, 2.0, 3.0}};
  constexpr tensor L = {{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}};

  constexpr auto exact_dsigma_dp = dsigma_dp(p, v, L);
  constexpr auto exact_dsigma_dv = dsigma_dv(p, v, L);
  constexpr auto exact_dsigma_dL = dsigma_dL(p, v, L);

  auto stress = std::apply(sigma, make_dual(p, v, L));

  for_constexpr<3,3>([&](auto i, auto j){
    std::cout << abs(       exact_dsigma_dp[i][j] - std::get<0>(stress[i][j].gradient)) << std::endl;
    std::cout << abs(sqnorm(exact_dsigma_dv[i][j] - std::get<1>(stress[i][j].gradient))) << std::endl;
    std::cout << abs(sqnorm(exact_dsigma_dL[i][j] - std::get<2>(stress[i][j].gradient))) << std::endl;
  });

  {
    auto grad = get_gradient(stress);
    std::cout << std::get<0>(grad) - exact_dsigma_dp << std::endl;
    std::cout << std::get<1>(grad) - exact_dsigma_dv << std::endl;
    std::cout << std::get<2>(grad) - exact_dsigma_dL << std::endl;
  }

  constexpr auto incompressibility = [](auto /*p*/, auto /*v*/, auto L) {      
    return tr(L);
  };

  auto dilatation = std::apply(incompressibility, make_dual(p, v, L));

  {
    std::cout << "dilatation derivatives:" << std::endl;
    auto grad = get_gradient(dilatation);
    std::cout << std::get<0>(grad) << std::endl;
    std::cout << std::get<1>(grad) << std::endl;
    std::cout << std::get<2>(grad) << std::endl;
  }


  constexpr auto tuple_func = [](auto p, auto v, auto L) {      
    return std::tuple{
      rho * outer(v, v) + 2.0 * mu * sym(L) - p * I,
      tr(L)
    };
  };



  auto stress_and_dilatation = std::apply(tuple_func, make_dual(p, v, L));

  {
    auto grad = get_gradient(stress_and_dilatation);

    //std::tuple<
    //  std::tuple<tensor<double, 3, 3>, tensor<double, 3, 3, 3>, tensor<double, 3, 3, 3, 3> >, 
    //  std::tuple<zero, zero, tensor<double, 3, 3> > 
    //>

    std::cout << "stress derivatives:" << std::endl;
    std::cout << std::get<0>(std::get<0>(grad)) << std::endl;
    std::cout << std::get<1>(std::get<0>(grad)) << std::endl;
    std::cout << std::get<2>(std::get<0>(grad)) << std::endl;

    std::cout << "stress derivative errors:" << std::endl;
    std::cout << std::get<0>(std::get<0>(grad)) - exact_dsigma_dp << std::endl;
    std::cout << std::get<1>(std::get<0>(grad)) - exact_dsigma_dv << std::endl;
    std::cout << std::get<2>(std::get<0>(grad)) - exact_dsigma_dL << std::endl;

    std::cout << "dilatation derivatives:" << std::endl;
    std::cout << std::get<0>(std::get<1>(grad)) << std::endl;
    std::cout << std::get<1>(std::get<1>(grad)) << std::endl;
    std::cout << std::get<2>(std::get<1>(grad)) << std::endl;
  }

}