#include "mfem.hpp"

void foo(std::vector<int>&);

void bar(std::vector<int>& vec)
{
  foo(vec);
}

#include "serac/serac_config.hpp" // for SERAC_REPO_DIR

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  std::ifstream imesh(SERAC_REPO_DIR"/data/meshes/beam-hex.mesh");
  mfem::Mesh serial(imesh, 1, 1, true);
  mfem::ParMesh parallel(MPI_COMM_WORLD, serial);
  std::cout << parallel.GetNE() << "\n";
  MPI_Finalize();
}
