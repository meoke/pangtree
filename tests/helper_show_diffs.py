def show_pangraph_differences(actual_pangraph, expected_pangraph):
    print("Nodes differences: ")
    show_nodes_differences(actual_pangraph._nodes, expected_pangraph._nodes)

    print("Paths differences: ")
    show_paths_differences(actual_pangraph._pathmanager.paths,
                                expected_pangraph._pathmanager.paths)


def show_nodes_differences(actual_nodes, expected_nodes):
    if len(actual_nodes) != len(expected_nodes):
        print(f"Actual graph has {len(actual_nodes)} nodes while expected graph: {len(expected_nodes)}")
    for node_id in range(len(actual_nodes)):
        if actual_nodes[node_id] != expected_nodes[node_id]:
            print(f"Nodes {node_id} differ:")
            print(f"Actual: {actual_nodes[node_id]}")
            print(f"Expected: {expected_nodes[node_id]}")


def show_paths_differences(actual_paths, expected_paths):
    if actual_paths.shape != expected_paths.shape:
        print(f"Actual pm has shape {actual_paths.shape} while expected pm: {expected_paths.shape}")

    for i, (actual_row, expected_row) in enumerate(zip(actual_paths, expected_paths)):
        if (actual_row != expected_row).any():
            print(f"Rows {i} differ.")
            print(f"Actual row: {actual_row}")
            print(f"Expected row: {expected_row}")