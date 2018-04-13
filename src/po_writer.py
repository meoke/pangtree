import numpy as np
import toolkit as t
import math

# def save_as_po_old(poagraph, sources_IDs):
#     def write_introduction_data(output_po_file, active_nodes_count):
#         output_po_file.writelines('VERSION=' + poagraph.version + '\n')
#         output_po_file.writelines('NAME=' + poagraph.name + '\n')
#         output_po_file.writelines('TITLE=' + poagraph.title + '\n')
#         output_po_file.writelines('LENGTH=' + str(active_nodes_count) + '\n')
#         output_po_file.writelines('SOURCECOUNT=' + str(len(sources_IDs)) + '\n')
#
#     def write_source_sequences(output_po_file, nodes, sources_weights):
#         def get_source_info_lines(source, source_weight):
#             source_nodes_IDs = poagraph.ns[source.ID] == True
#             return ['SOURCENAME=' + source.name,
#                     '\nSOURCEINFO=',
#                     " ".join([str(np.sum(source_nodes_IDs)),
#                               str(nodes['temp_ID'][nodes['orig_ID']==np.argmax(source_nodes_IDs == True)][0]),
#                               str(source_weight),
#                               str(-1),
#                               str(source.title)]),
#                     '\n']
#
#         for i, orig_src_ID in enumerate(sources_IDs):
#             output_po_file.writelines(get_source_info_lines(poagraph.sources[orig_src_ID], sources_weights[i]))
#
#     def write_nodes(output_po_file, nodes, sources):
#         def get_node_info(node, active_nodes):
#             # L_to_return = ['L' + str((nodes['temp_ID'][nodes['orig_ID'] == in_nodeID])[0]) for in_nodeID in
#             #                node.in_nodes if in_nodeID in active_nodes]
#
#             append_letter = lambda letter, ID : letter + str(ID)
#             append_letter_vec = np.vectorize(append_letter, otypes=[str])
#             active_in_nodes = np.intersect1d(active_nodes, node.in_nodes)
#             L_to_return = append_letter_vec('L', nodes['temp_ID'][active_in_nodes])
#
#             # todo czemu sources nie potrzebują przemapowania???
#             node_active_sources_IDs = np.where(poagraph.ns[sources_IDs, node.ID]==True)[0]
#             # S_to_return = ['S' + str(src_ID) for src_ID in node_active_sources_IDs]
#             S_to_return = append_letter_vec('S', node_active_sources_IDs)
#
#             n = node
#             all_aligned = []
#             while n.aligned_to != None and n.aligned_to != node.ID:
#                 if nodes['active'][nodes['orig_ID']==n.aligned_to]:
#                     all_aligned.append(n.aligned_to)
#                 n = poagraph.nodes[n.aligned_to]
#             new_aligned = all_aligned[0] if all_aligned else None
#
#             if new_aligned in active_nodes:
#                 A_to_return = "A" + str(nodes['temp_ID'][new_aligned])
#             else:
#                 A_to_return = ""
#
#             return [node.base.lower(), ':',
#                     "".join(L_to_return),
#                     "".join(S_to_return),
#                     A_to_return,
#                     "\n"]
#
#         active_nodes_IDs = nodes['orig_ID'][nodes['active']]
#         for node in poagraph.nodes:
#             if node.ID in active_nodes_IDs:
#                 output_po_file.writelines(get_node_info(node, active_nodes_IDs))
#
#     def get_partial_sources_weights(sources_IDs, nodes):
#         # return np.ones(shape=(len(sources_IDs)))
#         #todo prawdopodobnie przesunąć do poagraph - bo jest też wykorzystywane gdzie indziej
#         def get_source_weight(source):
#             #todo efficiency src_nodes_IDs = np.argwhere(poagraph.ns[source] == True)
#             src_nodes_IDs = np.where(poagraph.ns[source])[0]
#             # old = np.mean(np.array([nodes['sources_count'][nodes['orig_ID'] == node_ID][0] for node_ID in src_nodes_IDs]))
#             # tab = nodes['sources_count'][src_nodes_IDs]
#             # new = np.mean(nodes['sources_count'][src_nodes_IDs])
#             # return np.mean(np.array([nodes['sources_count'][nodes['orig_ID'] == node_ID][0] for node_ID in src_nodes_IDs]))
#             return np.mean(nodes['sources_count'][src_nodes_IDs])
#
#         def normalize_weight(weight, min_weight, diff):
#             # if weight == -1:
#             #     return -1
#             if diff == 0:
#                 return 100
#             return int((weight - min_weight)/diff*100)
#             # return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)
#
#         v_get = np.vectorize(get_source_weight)
#         weights = v_get(np.array(sources_IDs))
#         # weights = [*map(lambda source: get_source_weight(source), sources_IDs)]
#         max_weight = np.max(weights)
#         # min_weight = np.min(set(weights) - set([-1]))
#         min_weight = np.min(weights)
#         diff = max_weight - min_weight
#
#         v_norm = np.vectorize(normalize_weight)
#         normalized_weights = v_norm(weights, min_weight, diff)
#         # normalized_weights = [*map(lambda weight: normalize_weight(weight, ml, diff), weights)]
#         return normalized_weights
#
#     nodes = np.zeros(shape=(len(poagraph.nodes)), dtype=[('orig_ID', np.int32),
#                                                          ('temp_ID', np.int32),
#                                                          ('active', np.bool),
#                                                          ('sources_count', np.uint16)])
#
#     nodes['temp_ID'] = -1
#     for i, src_ID in enumerate(sources_IDs):
#         src_nodes_IDs = poagraph.ns[src_ID] == True
#         nodes['active'][src_nodes_IDs] = True
#         nodes['sources_count'][src_nodes_IDs] = nodes['sources_count'][src_nodes_IDs] + 1
#
#     active_nodes_count = len(nodes[nodes['active'] == True])
#     nodes['orig_ID'][:] = range(len(poagraph.nodes))
#     nodes['temp_ID'][nodes['active'] == True] = range(len(nodes[nodes['active'] == True]))
#
#     sources = np.zeros(shape=(len(poagraph.sources)), dtype=[('orig_ID', np.uint32),
#                                                              ('temp_ID', np.uint32),
#                                                              ('active', np.bool)])
#
#     sources['active'][sources_IDs] = True
#     sources['orig_ID'][:] = range(len(poagraph.sources))
#     sources['temp_ID'][sources['active']] = range(len(sources_IDs))
#
#     po_file_name = t.get_next_child_file_name(poagraph.path, poagraph.name + str(".po"))
#
#     sources_weights = get_partial_sources_weights(sources_IDs, nodes)
#     with open(po_file_name, 'w') as output_po_file:
#         write_introduction_data(output_po_file, active_nodes_count)
#         write_source_sequences(output_po_file, nodes, sources_weights)
#         write_nodes(output_po_file, nodes, sources)
#
#     return po_file_name, nodes

def save_as_po(poagraph, sources_IDs):
    def get_introduction_data(active_nodes_count):
        introduction_data_lines = []
        introduction_data_lines.append('VERSION=' + poagraph.version + '\n')
        introduction_data_lines.append('NAME=' + poagraph.name + '\n')
        introduction_data_lines.append('TITLE=' + poagraph.title + '\n')
        introduction_data_lines.append('LENGTH=' + str(active_nodes_count) + '\n')
        introduction_data_lines.append('SOURCECOUNT=' + str(len(sources_IDs)) + '\n')
        return introduction_data_lines

    def get_source_sequences(nodes, sources_weights):
        sources_sequences_lines = []
        def get_source_info_lines(source, source_weight):
            source_nodes_IDs = poagraph.ns[source.ID] == True
            return ['SOURCENAME=' + source.name + '\n',
                    'SOURCEINFO=' +
                    " ".join([str(np.sum(source_nodes_IDs)),
                              str(nodes['temp_ID'][nodes['orig_ID']==np.argmax(source_nodes_IDs == True)][0]),
                              str(source_weight),
                              str(-1),
                              str(source.title)]) + '\n']

        for i, orig_src_ID in enumerate(sources_IDs):
            sources_sequences_lines.extend(get_source_info_lines(poagraph.sources[orig_src_ID], sources_weights[i]))

        return sources_sequences_lines

    def get_nodes(nodes, sources):
        def get_node_info(node, active_nodes):
            active_in_nodes = []
            for in_node_ID in node.in_nodes:
                if nodes[in_node_ID]['active']:
                    active_in_nodes.append(in_node_ID)

            if len(active_in_nodes):
                L_to_return = 'L' + 'L'.join(nodes['temp_ID'][np.array(active_in_nodes)].astype('str', copy=False))
            else:
                L_to_return = ''

            # todo czemu sources nie potrzebują przemapowania???
            node_active_sources_IDs = np.where(poagraph.ns[sources_IDs, node.ID]==True)[0]
            S_to_return = 'S' + 'S'.join(node_active_sources_IDs.astype('str', copy = False))

            n = node
            all_aligned = []
            while n.aligned_to != None and n.aligned_to != node.ID:
                if nodes[n.aligned_to]['active']:
                    all_aligned.append(n.aligned_to)
                n = poagraph.nodes[n.aligned_to]
            A_to_return = "A" + str(nodes['temp_ID'][all_aligned[0]]) if all_aligned else ''

            # return "".join([node.base.lower(), ':',
            return "".join([node.base, ':',
                            L_to_return,
                            S_to_return,
                            A_to_return,
                            "\n"])

        active_nodes_IDs = nodes['orig_ID'][nodes['active']]
        nodes_lines = [None] * len(active_nodes_IDs)

        for i, orig_node_ID in enumerate(active_nodes_IDs):
            nodes_lines[i] = get_node_info(poagraph.nodes[orig_node_ID], active_nodes_IDs)
        return nodes_lines

    def get_partial_sources_weights(sources_IDs, nodes):
        #todo prawdopodobnie przesunąć do poagraph - bo jest też wykorzystywane gdzie indziej
        def get_source_weight(source):
            src_nodes_IDs = np.where(poagraph.ns[source])[0]
            return np.mean(nodes['sources_count'][src_nodes_IDs])

        def normalize_weight(weight, min_weight, diff):
            if diff == 0:
                return 100
            return int((weight - min_weight)/diff*100)

        v_get = np.vectorize(get_source_weight)
        weights = v_get(np.array(sources_IDs))
        max_weight = np.max(weights)
        min_weight = np.min(weights)
        diff = max_weight - min_weight

        v_norm = np.vectorize(normalize_weight)
        normalized_weights = v_norm(weights, min_weight, diff)
        return normalized_weights

    nodes = np.zeros(shape=(len(poagraph.nodes)), dtype=[('orig_ID', np.int32),
                                                         ('temp_ID', np.int32),
                                                         ('active', np.bool),
                                                         ('sources_count', np.uint16)])

    nodes['temp_ID'] = -1
    for i, src_ID in enumerate(sources_IDs):
        src_nodes_IDs = poagraph.ns[src_ID] == True
        nodes['active'][src_nodes_IDs] = True
        nodes['sources_count'][src_nodes_IDs] = nodes['sources_count'][src_nodes_IDs] + 1

    active_nodes_count = len(nodes[nodes['active'] == True])
    nodes['orig_ID'][:] = range(len(poagraph.nodes))
    nodes['temp_ID'][nodes['active'] == True] = range(len(nodes[nodes['active'] == True]))

    sources = np.zeros(shape=(len(poagraph.sources)), dtype=[('orig_ID', np.uint32),
                                                             ('temp_ID', np.uint32),
                                                             ('active', np.bool)])

    sources['active'][sources_IDs] = True
    sources['orig_ID'][:] = range(len(poagraph.sources))
    sources['temp_ID'][sources['active']] = range(len(sources_IDs))

    po_file_name = t.get_next_child_file_name(poagraph.path, poagraph.name + ".po")

    sources_weights = get_partial_sources_weights(sources_IDs, nodes)
    po_lines = []
    po_lines.extend(get_introduction_data(active_nodes_count))
    po_lines.extend(get_source_sequences(nodes, sources_weights))
    n = get_nodes(nodes, sources)
    po_lines.extend(n)

    with open(po_file_name, 'w') as out:
        out.writelines(po_lines)

    return po_file_name, nodes