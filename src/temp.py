    column_nodes = {}
    nodes_count = 0
    source_to_its_last_node_ID = {source_ID: -1 for source_ID, _ in enumerate(sources)}

    for column_ID in range(nucleotides_matrix.shape[0]):
        for row_ID in range(nucleotides_matrix.shape[1]):

            nucl = nucleotides_matrix[column_ID,row_ID]
            if nucl == nucleotides.code('-'):
                 continue

            if nucl not in column_nodes:
                column_nodes[nucl] = Node(ID=nodes_count, base=nucleotides.decode(nucl))
                nodes_count += 1

            node_id = column_nodes[nucl].ID

            if source_to_its_last_node_ID[row_ID] != -1:
                column_nodes[nucl].add_in_node([source_to_its_last_node_ID[row_ID]])
            source_to_its_last_node_ID[row_ID] = node_id

            try:
                poagraph.ns[row_ID, node_id] = True
            except:
                # pass
                raise ValueError("Wyjątek w maf readerze - kiedy to leci???")

            #progress = round(100 * (column_ID*nucleotides_matrix.shape[1]+row_ID+1) / all_nucles_count, 2) # todo logging
            #progress = column_ID #/ nucleotides_matrix.shape[0];
            #print("\r\tMultialignment ready in {0}".format(str(progress)), end='') # todo logging

        sorted_column_nodes = sorted([*column_nodes.values()], key=lambda node: node.ID)
        for i, node in enumerate(sorted_column_nodes):
            if not len(sorted_column_nodes) == 1:
                node.aligned_to = sorted_column_nodes[(i+1) % len(column_nodes)].ID
            poagraph.nodes[node.ID] = node

        column_nodes = {}