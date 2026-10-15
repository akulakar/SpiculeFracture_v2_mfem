#pragma once
// Result writer: ADIOS2 BP5 (default) or ParaView, plus resolved input and run details.

#include <memory>
#include <mfem.hpp>

#include "spicule/Config.hpp"
#include "spicule/Setup.hpp"

namespace spicule
{
    class OutputWriter
    {
    public:
        OutputWriter(const Config &cfg, mfem::ParMesh &mesh, Fields &fields, MPI_Comm comm);
        // Writes the accepted fields u and d.

        void Save(int cycle, double time);

    private:
        mfem::ParMesh &mesh_;
        Fields &fields_;
        std::unique_ptr<mfem::adios2stream> aos_;
        std::unique_ptr<mfem::ParaViewDataCollection> pv_;
    };

} // namespace spicule
