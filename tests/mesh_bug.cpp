#include "mfem.hpp"

void foo(std::vector<int>&);

void bar(std::vector<int>& vec)
{
  foo(vec);
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  std::ifstream imesh("/home/essman1/sandbox/serac/data/meshes/beam-hex.mesh");
  mfem::Mesh serial(imesh, 1, 1, true);
  mfem::ParMesh parallel(MPI_COMM_WORLD, serial);
  std::cout << parallel.GetNE() << "\n";
  MPI_Finalize();
}
