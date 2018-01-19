import po_reader as po_reader
from subprocess import run
from Errors import NoConsensusFound
import toolkit as t

def process_tree_node(poagraph, poagraphref, hbmin, cutoff_search_range, multiplier, output_dir):
    children_nodes = []
    while True:
        consensus = get_top_consensus(poagraph, poagraphref.sources_IDs, hbmin, output_dir)
        break
        # compatibilities_to_tree_nodes_sources = self.calc_compatibility(consensus, tree_node)
        # cutoff_for_max = self.find_cutoff_for_max(compatibilities_to_tree_nodes_sources, comp_range)
        # max_compatible_sources_IDs = self.get_compatible(tree_node, compatibilities_to_tree_nodes_sources, cutoff_for_max)
        # max_consensus = self.get_top_consensus(max_compatible_sources_IDs, hbmin)
        # compatibilities_to_tree_nodes_sources = self.calc_compatibility(max_consensus, tree_node)
        # cutoff_for_node = self.find_cutoff_for_node(compatibilities_to_tree_nodes_sources, multiplier)
        #
        # compatible_sources_IDs = self.get_compatible(tree_node, compatibilities_to_tree_nodes_sources, cutoff_for_node)
        #
        # poagraph_ref = POAGraphRef(compatible_sources_IDs, max_consensus, min(compatibilities_to_tree_nodes_sources))
        #
        # if children_nodes:
        #     decide_with = self.get_min_treshold(children_nodes) #children_nodes[-1].min_compatibility wwszystkie children nodes
        #     if decide_with < poagraph_ref.min_compatibility:
        #         poagraph_ref = POAGraphRef(tree_node.sources_IDs, max_consensus, poagraph_ref.min_compatibility)
        #         children_nodes.append(poagraph_ref)
        #         break
        #     else:
        #         tree_node.sources_IDs = np.array([src_ID for src_ID in tree_node.sources_IDs if src_ID not in compatible_sources_IDs])
        # else:
        #     children_nodes.append(poagraph_ref)


def get_top_consensus(poagraph, sources_IDs, hbmin, output_dir):
    po_file_path, nodes_map = poagraph.save_as_po(output_dir, sources_IDs)
    hb_file_path = t.change_file_extension(po_file_path, '.hb')
    run(['../bin/poa', '-read_msa', po_file_path, '-hb', '-po', hb_file_path, '../bin/blosum80.mat', '-hbmin',
         str(hbmin)])
    try:
        consensus0 = po_reader.read_consensus(hb_file_path, consensusID=0)
    except NoConsensusFound:
        raise NoConsensusFound()

    consensus0.nodes_IDs = [*map(lambda nodeID: nodes_map['org_nodeID'][nodes_map['temp_nodeID'] == nodeID][0], consensus0.nodes_IDs)]
    #consensus0.nodes_IDs = [*map(lambda nodeID: new_to_original_nodes_IDs[consensus0.nodes_IDs[i]])]

    return consensus0


def calc_compatibility(self, consensus, tree_node):
    pass


def get_min_treshold(self, children_nodes):
    pass


def get_compatible(self, tree_node, compatibilities_to_tree_nodes_sources, cutoff_for_node):
    pass


def find_cutoff_for_node(self, compatibilities_to_tree_nodes_sources, multiplier):
    pass


def get_compatible(self, tree_node, compatibilities_to_tree_nodes_sources, cutoff_for_max):
    pass


def find_cutoff_for_max(self, compatibilities_to_tree_nodes_sources, comp_range):
    pass
