#include "spicule/Config.hpp"

#include <fstream>
#include <stdexcept>

namespace spicule
{
    namespace
    {
        // Read input
        json ReadJson(const std::string &path)
        {
            std::ifstream in(path);
            if (!in.is_open())
                throw std::runtime_error("Could not open the input file: " + path);
            json input;
            in >> input;
            return input;
        }

        // Get value. Has to exist in json, no fallback
        template <class T>
        T Get(const json &input, const char *section, const char *key)
        {
            if (!input.contains(section) || !input.at(section).contains(key))
                throw std::runtime_error(std::string("Missing input key: ") + section + "." + key);
            else
                return input.at(section).at(key).get<T>();
        }

        template <class T>
        T Get(const json &input, const char *section, const char *sub_section, const char *key)
        {
            if (!input.contains(section) || !input.at(section).contains(sub_section) || !input.at(section).at(sub_section).contains(key))
                throw std::runtime_error(std::string("Missing input key: ") + section + "." + key);
            else
                return input.at(section).at(sub_section).at(key).get<T>();
        }

        // Get value. Need not exist in json, there is fallback
        template <class T>
        T GetOr(const json &input, const char *section, const char *key, T fallback)
        {
            if (input.contains(section) && input.at(section).contains(key))
                return input.at(section).at(key).get<T>();
            else
                return fallback;
        }

        template <class T>
        T GetOr(const json &input, const char *section, const char *sub_section, const char *key, T fallback)
        {
            if (input.contains(section) && input.at(section).contains(sub_section) && input.at(section).at(sub_section).contains(key))
                return input.at(section).at(sub_section).at(key).get<T>();
            else
                return fallback;
        }
    } // unnamed namespace

    Config LoadConfig(const std::string &path, const std::string &overridePath)
    {
        json input = ReadJson(path);
        if (!overridePath.empty())
            input.merge_patch(ReadJson(overridePath));

        Config config;
        config.inputJson = input;
        config.testName = input.at("Test_Name").get<std::string>();

        config.mesh.file = Get<std::string>(input, "Mesh", "file");
        config.mesh.order = GetOr(input, "Mesh", "order", 1);
        config.mesh.refLevels = GetOr(input, "Mesh", "refinement_levels", 0);

        config.fe.order = GetOr(input, "Finite_Element", "order", config.mesh.order);
        config.fe.basis = GetOr<std::string>(input, "Finite_Element", "basis", "GaussLobatto");
        config.fe.ordering = GetOr<std::string>(input, "Finite_Element", "ordering", "ByNODES");
        config.fe.stressOrder = GetOr(input, "Finite_Element", "stress_order", 0);

        config.material.Ey = Get<double>(input, "Material", "youngs_modulus");
        config.material.nu = Get<double>(input, "Material", "poisson_ratio");
        config.material.Gc = Get<double>(input, "Material", "fracture_toughness");
        config.material.eps = Get<double>(input, "Material", "length_scale");
        config.material.Kepsilon = Get<double>(input, "Material", "K_epsilon");

        config.eigenstrain.type = GetOr<std::string>(input, "Eigenstrain", "type", "uniform");
        config.eigenstrain.beta = GetOr(input, "Eigenstrain", "swelling_coefficient", 0.0);
        config.eigenstrain.deltaC = GetOr(input, "Eigenstrain", "Delta_concentration", 0.0);

        // boundary conditions
        // for each unique dirichlet boundary
        for (const auto &dirichlet_bc : input.at("Boundary_Conditions").at("dirichlet_attributes"))
        {
            BCEntry entry;

            entry.attrs = dirichlet_bc.at("attributes").get<std::vector<int>>();
            entry.components = dirichlet_bc.at("components").get<std::vector<int>>();
            if (dirichlet_bc.contains("value") && dirichlet_bc.at("value").is_string())
                entry.isNonZero = (dirichlet_bc.at("value").get<std::string>() == "nonzero_displacement");
            else
                entry.value = dirichlet_bc.value("value", 0.0);

            config.bcs.push_back(std::move(entry));
        }

        config.load.finalDisp = Get<double>(input, "Boundary_Conditions", "nonzero_displacement", "final_displacement");
        config.load.minInc = Get<double>(input, "Boundary_Conditions", "nonzero_displacement", "min_increment");
        config.load.maxInc = Get<double>(input, "Boundary_Conditions", "nonzero_displacement", "max_increment");

        // staggered iterations settings
        config.staggered.tol = GetOr(input, "Staggered_Solve", "tolerance", 1e-6);
        config.staggered.maxIter = GetOr(input, "Staggered_Solve", "max_iterations", 1000);

        // linear solver settings
        // elasticity (disp) solver
        config.dispSolver.library = GetOr<std::string>(input, "Solvers", "elasticity_solver", "library", "mfem");
        config.dispSolver.type = GetOr<std::string>(input, "Solvers", "elasticity_solver", "type", "pcg");
        config.dispSolver.preconditioner = GetOr<std::string>(input, "Solvers", "elasticity_solver", "preconditioner", "amg_systems");
        config.dispSolver.tol = GetOr(input, "Solvers", "elasticity_solver", "tolerance", 1.0e-6);
        config.dispSolver.maxIter = GetOr(input, "Solvers", "elasticity_solver", "max_iterations", 1000);

        // damage solver
        config.dmgSolver.library = GetOr<std::string>(input, "Solvers", "damage_solver", "library", "mfem");
        config.dmgSolver.type = GetOr<std::string>(input, "Solvers", "damage_solver", "type", "pcg");
        config.dmgSolver.preconditioner = GetOr<std::string>(input, "Solvers", "damage_solver", "preconditioner", "amg");
        config.dmgSolver.tol = GetOr(input, "Solvers", "damage_solver", "tolerance", 1.0e-6);
        config.dmgSolver.maxIter = GetOr(input, "Solvers", "damage_solver", "max_iterations", 1000);

        // output settings
        config.output.saveFrequency = GetOr<int>(input, "Output", "save_frequency", 1);
        config.output.format = GetOr<std::string>(input, "Output", "format", "adios2");
        config.output.path = GetOr<std::string>(input, "Output", "path", "../results");

        return config;
    }

} // namespace spicule
