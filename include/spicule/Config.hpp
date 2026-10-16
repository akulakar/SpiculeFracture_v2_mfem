#pragma once
// Typed view of the JSON input. Every other module takes only its own section.

#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace spicule
{
    using json = nlohmann::json;

    struct MeshParams
    {
        std::string file;
        int order = 1;
        int refLevels = 0;
    };

    struct FEParams
    {
        int order = 1;
        std::string basis = "GaussLobatto"; // GaussLobatto | Positive | ClosedUniform
        std::string ordering = "byNODES";
        int stressOrder = 0; // L2 order for stress/strain output
    };

    struct MaterialParams
    {
        double lame = 0.0, shear = 0.0; // inputs
        double Ey = 0.0, nu = 0.0;      // derived
        double Gc = 0.0, eps = 0.0, Kepsilon = 0.0;
    };

    struct EigenstrainParams
    {
        std::string type = "uniform";
        double beta = 0.0, deltaC = 0.0, scale = 1.0;
    };

    struct BCEntry
    {
        std::vector<int> attrs;      // boundary attributes (0-based)
        std::vector<int> components; // displacement components that are Dirichlet
        bool isNonZero = false;      // true -> value follows the load controller
        double value = 0.0;          // default 0.0 value for Dirichlet boundary unles isNonZero = true
    };

    struct LoadParams
    {
        double finalDisp = 0.0, minInc = 0.0, maxInc = 0.0;
    };

    struct StaggeredParams
    {
        double tol = 1e-6;
        int maxIter = 1000;
    };

    struct LinearSolverParams
    {
        std::string library = "mfem"; // mfem, hypre, petsc
        std::string type = "pcg";     // pcg or gmres
        std::string preconditioner = "amg";
        double tol = 1e-6;
        int maxIter = 300;
    };

    struct OutputParams
    {
        int saveFrequency = 1; // save every file
        std::string format = "adios2";
        std::string path = "../results";
    };

    struct Config
    {
        std::string testName;
        MeshParams mesh;
        FEParams fe;
        MaterialParams material;
        EigenstrainParams eigenstrain;
        std::vector<BCEntry> bcs;
        LoadParams load;
        StaggeredParams staggered;
        LinearSolverParams dispSolver, dmgSolver;
        OutputParams output;
        json inputJson; // resolved JSON, dumped next to results for reproducibility
    };

    // Loads `path`, optionally merges `overridePath` on top (RFC 7386 merge-patch),
    // validates required keys, and fills the typed Config. Throws std::runtime_error.
    Config LoadConfig(const std::string &path, const std::string &overridePath = "");

} // namespace spicule
