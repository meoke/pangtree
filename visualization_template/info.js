var info =
{  name: 'ebola_full',
   consensus_algorithm: 'tree',
   running_time: '33min',
   sources_count: 3,
   nodes_count: 13,
   sequences_per_node: 1.84,
   levels:[0.6,0.8]
};

var treeData = [
      {
        name: "Top Level",
        parent: "null",
        children: [
          {
            "name": "Level 2: A",
            "parent": "Top Level",
            "val": "wartosc",
            "children": [
              {
                "name": "Son of A",
                "parent": "Level 2: A"
              },
              {
                "name": "Daughter of A",
                "parent": "Level 2: A"
              }
            ]
          },
          {
            "name": "Level 2: B",
            "parent": "Top Level"
          }
        ]
      }
    ];