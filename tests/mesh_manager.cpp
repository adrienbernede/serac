#include <optional>
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"
#include "mfem.hpp"

void foo(std::vector<int>&) {}

static std::optional<axom::sidre::MFEMSidreDataCollection> datacoll_;
// static std::optional<mfem::ConduitDataCollection> datacoll_;
