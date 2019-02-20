from Bio import AlignIO
from graph import nucleotides
from graph.Node import Node
from graph.PangraphBuilder.PangraphBuilderBase import PangraphBuilderBase
from metadata.MultialignmentMetadata import MultialignmentMetadata


class PangraphBuilderFromMAF(PangraphBuilderBase):
    def __init__(self, genomes_info: MultialignmentMetadata):
        super().__init__(genomes_info)

    def build(self, input, pangraph):
        input = [*AlignIO.parse(input, "maf")]
        nodes_count = PangraphBuilderFromMAF.get_nodes_count(input)
        pangraph._nodes = [None] * nodes_count
        pangraph._pathmanager.init_paths(self.sequences_names, nodes_count)

        sequence_last_node_id = {seq_id: None for seq_id in self.sequences_names}
        current_node_id = -1
        for block in input:
            block_width = len(block[0].seq)

            for col in range(block_width):
                sequence_name_to_nucleotide = {seq.id: seq[col] for seq in block}
                nodes_codes = sorted([*(
                    set([nucleotide for nucleotide in sequence_name_to_nucleotide.values()])).difference({'-'})])
                column_nodes_ids = [current_node_id + i + 1 for i, _ in enumerate(nodes_codes)]

                for i, nucl in enumerate(nodes_codes):
                    current_node_id += 1
                    node = Node(id=current_node_id,
                                base=nucleotides.code(nucl),
                                in_nodes=set(),
                                aligned_to=PangraphBuilderFromMAF.get_next_aligned_node_id(i, column_nodes_ids))

                    for sequence, nucleotide in sequence_name_to_nucleotide.items():
                        if nucleotide == nucl:
                            pangraph.add_path_to_node(path_name=sequence, node_id=current_node_id)
                            last_node_id = sequence_last_node_id[sequence]
                            if last_node_id is not None:
                                node.in_nodes.add(last_node_id)
                            sequence_last_node_id[sequence] = current_node_id
                    node.in_nodes = list(node.in_nodes)
                    pangraph._nodes[current_node_id] = node
        pangraph.remove_empty_paths()

    @staticmethod
    def get_nodes_count(mafalignment) -> int:
        nodes_count = 0
        for block in mafalignment:
            number_of_columns = len(block[0].seq)

            for col_id in range(number_of_columns):
                letters_in_columns = set([block[i].seq[col_id] for i in range(len(block))]).difference(set('-'))
                nodes_count += len(letters_in_columns)

        return nodes_count

    @staticmethod
    def get_next_aligned_node_id(current_column_i, column_nodes_ids):
        if len(column_nodes_ids) > 1:
            return column_nodes_ids[(current_column_i + 1) % len(column_nodes_ids)]
        return None
