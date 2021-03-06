#include "masb_nodes.hpp"

using namespace geoflow::nodes::mat;

void register_nodes(geoflow::NodeRegister& node_register) {
    node_register.register_node<ComputeMedialAxisNode>("ComputeMedialAxisNode");
    node_register.register_node<ComputeNormalsNode>("ComputeNormalsNode");
    node_register.register_node<SegmentMakerNode>("SegmentMaker");
    node_register.register_node<TestPointsNode>("TestPoints");
    node_register.register_node<RegionGrowMedialAxisNode>("RegionGrowMedialAxis");
    node_register.register_node<MATCSVWriterNode>("MATCSVWriter");
    node_register.register_node<MATCSVLoaderNode>("MATCSVLoader");
    node_register.register_node<PLYReaderNode>("PLYReader");
    node_register.register_node<PLYWriterNode>("PLYWriter");
    node_register.register_node<SplitMATInteriorExteriorNode>("SplitMATInteriorExterior");
}