#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pyrbd_plusplus/mcs.hpp>
#include <pyrbd_plusplus/pathset.hpp>
#include <pyrbd_plusplus/sdp.hpp>

namespace py = pybind11;
using namespace pyrbdpp;
using namespace pyrbdpp::mcs;
using namespace pyrbdpp::pathset;
using namespace pyrbdpp::sdp;

PYBIND11_MODULE(pyrbd_plusplus, m)
{
    m.doc() = "PyRBD++ - Reliability Block Diagram analysis library";

    // Bind Classes
    py::class_<ProbabilityMap>(m, "ProbabilityMap")
        .def(py::init<const std::map<int, double>&>())
        .def("print", &ProbabilityMap::print)
        .def("__getitem__", &ProbabilityMap::operator[]);

    py::class_<SDP>(m, "SDP")
        .def(py::init<bool, std::vector<int>>())
        .def("isComplementary", &SDP::isComplementary)
        .def("getSet", &SDP::getSet, py::return_value_policy::reference_internal);

    // MCS Algorithm
    auto mcs_mod = m.def_submodule("mcs", "Module for MCS algorithm");
    mcs_mod.doc() = "Module for MCS algorithm";

    mcs_mod.def("to_probaset", &mcs::toProbaSet,
                "Convert minimal cut sets to probability sets",
                py::arg("src"), py::arg("dst"), py::arg("min_cut_sets"));

    mcs_mod.def("to_probaset_debug", &mcs::toProbaSetDebug,
                "Debug Version: Convert minimal cut sets to probability sets",
                py::arg("src"), py::arg("dst"), py::arg("min_cut_sets"));

    mcs_mod.def("eval_avail", 
                [](NodeID src, NodeID dst, const std::map<int, double>& probabilities, const MinCutSets& min_cut_sets) {
                    ProbabilityMap probMap(probabilities); 
                    return mcs::evalAvail(src, dst, probMap, min_cut_sets);
                },
                "Evaluate availability for single source destination pair using MCS approach",
                py::arg("src"), py::arg("dst"), py::arg("probabilities"), py::arg("min_cut_sets"));

    mcs_mod.def("eval_avail_topo", 
                [](const std::vector<std::pair<NodeID, NodeID>>& node_pairs, 
                   const std::map<int, double>& probabilities, 
                   const std::vector<MinCutSets>& min_cut_sets_list) {
                    ProbabilityMap probMap(probabilities); 
                    return mcs::evalAvailTopo(node_pairs, probMap, min_cut_sets_list);
                },
                "Evaluate availability for each node pairs in topology using MCS (serial)",
                py::arg("node_pairs"), py::arg("probabilities"), py::arg("min_cut_sets_list"));
    
    mcs_mod.def("eval_avail_topo_parallel", 
                [](const std::vector<std::pair<NodeID, NodeID>>& node_pairs, 
                   const std::map<int, double>& probabilities, 
                   const std::vector<MinCutSets>& min_cut_sets_list) {
                    ProbabilityMap probMap(probabilities); 
                    return mcs::evalAvailTopoParallel(node_pairs, probMap, min_cut_sets_list);
                },
                "Evaluate availability for each node pairs in topology using MCS (parallel)",
                py::arg("node_pairs"), py::arg("probabilities"), py::arg("min_cut_sets_list"),
                py::call_guard<py::gil_scoped_release>());

    // PathSet Algorithm
    auto pathset_mod = m.def_submodule("pathset", "Module for PathSet algorithm");
    pathset_mod.doc() = "Module for PathSet algorithm";

    pathset_mod.def("to_probaset", &pathset::toProbaSet,
                "Convert path sets to probability sets",
                py::arg("src"), py::arg("dst"), py::arg("path_sets"));

    pathset_mod.def("to_probaset_debug", &pathset::toProbaSetDebug,
                "Debug Version: Convert path sets to probability sets",
                py::arg("src"), py::arg("dst"), py::arg("path_sets"));

    pathset_mod.def("eval_avail", 
                [](NodeID src, NodeID dst, const std::map<int, double>& probabilities, const PathSets& path_sets) {
                    ProbabilityMap probMap(probabilities); 
                    return pathset::evalAvail(src, dst, probMap, path_sets);
                },
                "Evaluate availability for single source destination pair using PathSet approach",
                py::arg("src"), py::arg("dst"), py::arg("probabilities"), py::arg("path_sets"));
    
    pathset_mod.def("eval_avail_topo", 
                [](const std::vector<std::pair<NodeID, NodeID>>& node_pairs, 
                   const std::map<int, double>& probabilities, 
                   const std::vector<PathSets>& pathsets_list) {
                    ProbabilityMap probMap(probabilities); 
                    return pathset::evalAvailTopo(node_pairs, probMap, pathsets_list);
                },
                "Evaluate availability for each node pairs in topology using PathSet (serial)",
                py::arg("node_pairs"), py::arg("probabilities"), py::arg("pathsets_list"));
    
    pathset_mod.def("eval_avail_topo_parallel", 
                [](const std::vector<std::pair<NodeID, NodeID>>& node_pairs, 
                   const std::map<int, double>& probabilities, 
                   const std::vector<PathSets>& pathsets_list) {
                    ProbabilityMap probMap(probabilities); 
                    return pathset::evalAvailTopoParallel(node_pairs, probMap, pathsets_list);
                },
                "Evaluate availability for each node pairs in topology using PathSet (parallel)",
                py::arg("node_pairs"), py::arg("probabilities"), py::arg("pathsets_list"),
                py::call_guard<py::gil_scoped_release>());
    
    // SDP Algorithm
    auto sdp_mod = m.def_submodule("sdp", "Module for SDP algorithm");
    sdp_mod.doc() = "Module for SDP algorithm";

    sdp_mod.def("to_sdp_set", &sdp::toSDPSet,
                "Convert path sets to SDP sets (serial)",
                py::arg("src"), py::arg("dst"), py::arg("path_sets"));
    
    sdp_mod.def("to_sdp_set_debug", &sdp::toSDPSetDebug,
                "Debug Version: Convert path sets to SDP sets",
                py::arg("src"), py::arg("dst"), py::arg("path_sets"));

    // sdp_mod.def("to_sdp_set_parallel", &sdp::toSDPSetParallel,
    //             "Convert path sets to SDP sets (parallel)",
    //             py::arg("src"), py::arg("dst"), py::arg("path_sets"));

    sdp_mod.def("eval_avail", 
                [](NodeID src, NodeID dst, const std::map<int, double>& probabilities, PathSets& path_sets) {
                    ProbabilityMap probMap(probabilities); 
                    return sdp::evalAvail(src, dst, probMap, path_sets);
                },
                "Evaluate availability for single source destination pair using SDP approach",
                py::arg("src"), py::arg("dst"), py::arg("probabilities"), py::arg("path_sets"));

    sdp_mod.def("eval_avail_parallel", 
                [](NodeID src, NodeID dst, const std::map<int, double>& probabilities, PathSets& path_sets) {
                    ProbabilityMap probMap(probabilities); 
                    return sdp::evalAvailParallel(src, dst, probMap, path_sets);
                },
                "Evaluate availability for single source destination pair using SDP approach (parallel)",
                py::arg("src"), py::arg("dst"), py::arg("probabilities"), py::arg("path_sets"),
                py::call_guard<py::gil_scoped_release>());
    
    sdp_mod.def("eval_avail_topo", 
                [](const std::vector<std::pair<NodeID, NodeID>>& node_pairs, 
                   const std::map<int, double>& probabilities, 
                   std::vector<PathSets>& pathsets_list) {
                    ProbabilityMap probMap(probabilities); 
                    return sdp::evalAvailTopo(node_pairs, probMap, pathsets_list);
                },
                "Evaluate availability for each node pairs in topology using SDP (serial)",
                py::arg("node_pairs"), py::arg("probabilities"), py::arg("pathsets_list"));

    sdp_mod.def("eval_avail_topo_parallel", 
                [](const std::vector<std::pair<NodeID, NodeID>>& node_pairs, 
                   const std::map<int, double>& probabilities, 
                   std::vector<PathSets>& pathsets_list) {
                    ProbabilityMap probMap(probabilities); 
                    return sdp::evalAvailTopoParallel(node_pairs, probMap, pathsets_list);
                },
                "Evaluate availability for each node pairs in topology using SDP (parallel)",
                py::arg("node_pairs"), py::arg("probabilities"), py::arg("pathsets_list"),
                py::call_guard<py::gil_scoped_release>());
}