import numpy as np
import toolkit as t

def save_as_po(poagraph, sources_IDs):
    def write_introduction_data(output_po_file, active_nodes_count):
        output_po_file.writelines('VERSION=' + poagraph.version + '\n')
        output_po_file.writelines('NAME=' + poagraph.name + '\n')
        output_po_file.writelines('TITLE=' + poagraph.title + '\n')
        output_po_file.writelines('LENGTH=' + str(active_nodes_count) + '\n')
        output_po_file.writelines('SOURCECOUNT=' + str(len(sources_IDs)) + '\n')

    def write_source_sequences(output_po_file, nodes, sources_weights):
        def get_source_info_lines(source, source_weight):
            source_nodes_IDs = poagraph.ns[source.ID] == True
            return ['SOURCENAME=' + source.name,
                    '\nSOURCEINFO=',
                    " ".join([str(sum(source_nodes_IDs)),
                              str(nodes['temp_ID'][nodes['orig_ID']==np.argmax(source_nodes_IDs == True)][0]),
                              str(source_weight),
                              str(-1),
                              str(source.title)]),
                    '\n']

        for i, orig_src_ID in enumerate(sources_IDs):
            output_po_file.writelines(get_source_info_lines(poagraph.sources[orig_src_ID], sources_weights[i]))

    def write_nodes(output_po_file, nodes, sources):
        def get_node_info(node, active_nodes):
            L_to_return = ['L' + str((nodes['temp_ID'][nodes['orig_ID'] == in_nodeID])[0]) for in_nodeID in
                           node.in_nodes if in_nodeID in active_nodes]
            node_active_sources_IDs = np.where(poagraph.ns[sources_IDs, node.ID]==True)[0]
            S_to_return = ['S' + str(src_ID) for src_ID in node_active_sources_IDs]
            n = node
            all_aligned = []
            while n.aligned_to != None and n.aligned_to != node.ID:
                if nodes['active'][nodes['orig_ID']==n.aligned_to]:
                    all_aligned.append(n.aligned_to)
                n = poagraph.nodes[n.aligned_to]
            new_aligned = all_aligned[0] if all_aligned else None

            if new_aligned in active_nodes:
                A_to_return = "A" + str(nodes['temp_ID'][new_aligned])
            else:
                A_to_return = ""

            return [node.base, ':',
                    "".join(L_to_return),
                    "".join(S_to_return),
                    A_to_return,
                    "\n"]

        active_nodes_IDs = nodes['orig_ID'][nodes['active'] == True]
        for node in poagraph.nodes:
            if node.ID in active_nodes_IDs:
                output_po_file.writelines(get_node_info(node, active_nodes_IDs))

    def get_partial_sources_weights(sources_IDs, nodes):
        def get_source_weight(source):
            src_nodes_IDs = np.argwhere(poagraph.ns[source] == True)
            return np.mean(np.array([nodes['sources_count'][nodes['orig_ID'] == node_ID][0] for node_ID in src_nodes_IDs]))

        def normalize_weight(weight, max_weight, min_weight):
            # if weight == -1:
            #     return -1
            if max_weight - min_weight == 0:
                return 100
            return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)

        weights = [*map(lambda source: get_source_weight(source), sources_IDs)]
        max_weight = max(weights)
        min_weight = min(set(weights) - set([-1]))

        normalized_weights = [*map(lambda weight: normalize_weight(weight, max_weight, min_weight), weights)]
        return normalized_weights

    nodes = np.zeros(shape=(len(poagraph.nodes)), dtype=[('orig_ID', np.uint32),
                                                         ('temp_ID', np.uint32),
                                                         ('active', np.bool),
                                                         ('sources_count', np.uint16)])

    for i, src_ID in enumerate(sources_IDs):
        src_nodes_IDs = poagraph.ns[src_ID] == True
        nodes['active'][src_nodes_IDs] = True
        nodes['sources_count'][src_nodes_IDs] = nodes['sources_count'][src_nodes_IDs] + 1

    active_nodes_count = len(nodes[nodes['active'] == True])
    # nodes['orig_ID'][nodes['active'] == True] = np.where(nodes['active'] == True)[0]
    nodes['orig_ID'][:] = range(len(poagraph.nodes))
    nodes['temp_ID'][nodes['active'] == True] = range(len(nodes[nodes['active'] == True]))

    sources = np.zeros(shape=(len(poagraph.sources)), dtype=[('orig_ID', np.uint32),
                                                             ('temp_ID', np.uint32),
                                                             ('active', np.bool)])

    sources['active'][sources_IDs] = True
    sources['orig_ID'][:] = range(len(poagraph.sources))# np.where(sources['active'] == True)[0]
    sources['temp_ID'][sources['active']] = range(len(sources_IDs))

    po_file_name = t.get_next_child_file_name(poagraph.path, poagraph.name + str(".po"))

    sources_weights = get_partial_sources_weights(sources_IDs, nodes)
    with open(po_file_name, 'w') as output_po_file:
        write_introduction_data(output_po_file, active_nodes_count)
        write_source_sequences(output_po_file, nodes, sources_weights)
        write_nodes(output_po_file, nodes, sources)

    return po_file_name, nodes