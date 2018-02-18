import toolkit as t
import maf_reader as maf_reader
import po_reader as po_reader
import consensus as cons
import POAGraphViz as viz
from Errors import CloseProgram

class Multialignment(object):
    def __init__(self, data_type='ebola'):
        self.name = None
        self.output_dir = None
        self.poagraphs = None
        self.data_type = data_type

    def build_multialignment_from_maf(self, maf_file_path, merge_option):
        """Build multialignment structrure from maf file.

        Keyword arguments:
        maf_file_path -- path to the maf file
        merge_option -- #todo uzupelnic poprawnie
        """
        print("Building multialignment from " + maf_file_path + "...") #todo logging

        self.name = self._get_multialignment_name(maf_file_path)
        self.output_dir = self._get_output_dir(maf_file_path)
        self.poagraphs = maf_reader.parse_to_poagraphs(file_path=maf_file_path,
                                                       merge_option=merge_option,
                                                       multialignment_name=self.name,
                                                       output_dir=self.output_dir)
        if self.poagraphs is None:
            raise CloseProgram("No poagraph was built.")

        for p in self.poagraphs:
            p.data_type = self.data_type

    def build_multialignment_from_po(self, po_file_name):
        """Build multialignment structrure with one poagraph from po file.

                Keyword arguments:
                po_file_name -- path to the po file
                """
        print("Building multialignment from " + po_file_name + "...") #todo logging

        self.name = self._get_multialignment_name(po_file_name)
        self.output_dir = self._get_output_dir(po_file_name)
        self.poagraphs = [po_reader.parse_to_poagraph(file_path=po_file_name,
                                                      output_dir=self.output_dir)]

    def _get_multialignment_name(self, input_file_name):
        """Get name for multialignment based on its filename.

            Keyword arguments:
            input_file_name -- path to the file with multialignment
            """
        return t.get_file_name_without_extension(input_file_name)

    def _get_output_dir(self, input_file_name):
        """Get directory path for operations on this multiaglinemnt output.

            Keyword arguments:
            input_file_name -- path to the file with multialignment
            """
        return t.create_next_sibling_dir(input_file_name, "result")

    def generate_consensus(self, option, hbmin, min_comp, cutoff_search_range, multiplier, stop, re_consensus):
        for i, p in enumerate(self.poagraphs):
            if option is 1:
                print('Generate consensuses (hbmin=', str(hbmin) + ') in one iteration...') # todo logging
            elif option is 3:
                print('Generate tree based consensus (hbmin=', str(hbmin),
                      ', min_comp=', str(min_comp),
                      ', range=', str(cutoff_search_range),
                      ', multiplier=', str(multiplier))

                root_poagraphref_ID = p.create_root_poagraphref()
                poagraph_IDs_to_process = [root_poagraphref_ID]
                while poagraph_IDs_to_process:
                    new_poagraph_IDs_to_process = []
                    for tree_node_ID in poagraph_IDs_to_process:
                        children_nodes_IDs = cons.process_tree_node(p,
                                                                    tree_node_ID,
                                                                    cutoff_search_range,
                                                                    multiplier,
                                                                    re_consensus,
                                                                    stop)
                        new_poagraph_IDs_to_process.extend(children_nodes_IDs)

                        poagraph_IDs_to_process = new_poagraph_IDs_to_process

    def generate_visualization(self, consensus_option, draw_poagraph_option, processing_time):
        vizualization_output_dir = t.create_child_dir(self.output_dir, "visualization")
        viz.generate_visualization(self,
                                   vizualization_output_dir,
                                   consensus_option,
                                   draw_poagraph_option,
                                   processing_time)



