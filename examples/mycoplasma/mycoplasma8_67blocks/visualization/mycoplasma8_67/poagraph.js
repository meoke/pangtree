var poagraph = {
  nodes: [
    {
      data: {
        id: 0,
        nucleobase: 'A',
        source: [0],
        weight: 0.33,
      },
      position: { x: 60, y: 15 },
      // classes: 'czerwony zielony'
    }, {
      data: {
        id: 1,
        nucleobase: 'T',
        source: [1,2],
        weight: 0.66,
      },
      position: { x: 60, y: -15 },
    },{
      data: {
        id: 2,
        nucleobase: 'G',
        source: [0],
        weight: 0.33
      },
      position: { x: 120, y: 15 },
    },{
      data: {
        id: 3,
        nucleobase: 'C',
        source: [1,2],
        weight: 0.66
      },
      position: { x: 120, y: -15 },
    },
    {
      data: {
        id: 4,
        nucleobase: 'C',
        source: [0,1,2],
        weight: 1
      },
      position: { x: 180, y: 0.5 },
    },{
      data: {
        id: 5,
        nucleobase: 'A',
        source: [0],
        weight: 0.33
      },
      position: { x: 240, y: 15 },
    },{
      data: {
        id: 6,
        nucleobase: 'T',
        source: [1],
        weight: 0.33
      },
      position: { x: 240, y: 0.5 },
    },{
      data: {
        id: 7,
        nucleobase: 'G',
        source: [2],
        weight: 0.33
      },
      position: { x: 240, y: -15 },
    },{
      data: {
        id: 8,
        nucleobase: 'G',
        source: [0, 1],
        weight: 0.66
      },
      position: { x: 300, y: 15 },
    },{
      data: {
        id: 9,
        nucleobase: 'C',
        source: [2],
        weight: 0.33
      },
      position: { x: 300, y: -15 },
    },{
      data: {
        id: 10,
        nucleobase: 'A',
        source: [0,1,2],
        weight: 1
      },
      position: { x: 360, y: 0 },
    },{
      data: {
        id: 11,
        nucleobase: 'A',
        source: [0,1,2],
        weight: 1
      },
      position: { x: 420, y: 0 },
    },{
      data: {
        id: 12,
        nucleobase: 'A',
        source: [0, 1, 2],
        weight: 1
      },
      position: { x: 480, y: 0 },
    }
  ],
  edges: [
    {
      data: {
        id: 13,
        source: 0,
        target: 2,
        weight: 0.33,
        consensus: -1,
        level: -1
      },
      classes: 'edge'
    },{
      data: {
        id: 14,
        source: 1,
        target: 3,
        consensus: -1,
        weight: 0.66,
        level: -1
      },
      classes: 'edge'
    },{
      data: {
        id: 15,
        source: 2,
        target: 4,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    },{
      data: {
        id: 16,
        source: 3,
        target: 4,
        consensus: -1,
        weight: 0.66,
        level: -1
      },
      classes: 'edge'
    },{
      data: {
        id: 17,
        source: 4,
        target: 5,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 18,
        source: 4,
        target: 6,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 19,
        source: 4,
        target: 7,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 20,
        source: 5,
        target: 8,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 21,
        source: 6,
        target: 8,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 22,
        source: 7,
        target: 9,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 23,
        source: 8,
        target: 10,
        consensus: -1,
        weight: 0.66,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 24,
        source: 9,
        target: 10,
        consensus: -1,
        weight: 0.33,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 25,
        source: 10,
        target: 11,
        consensus: -1,
        weight: 1,
        level: -1
      },
      classes: 'edge'
    }, {
      data: {
        id: 26,
        source: 11,
        target: 12,
        consensus: -1,
        weight: 1,
        level: -1
      },
      classes: 'edge'
    },
//    level 1
//    consensus 0
    {
      data: {
        id: 27,
        source: 0,
        target: 2,
        consensus: 0,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 28,
        source: 2,
        target: 4,
        consensus: 0,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 29,
        source: 4,
        target: 5,
        consensus: 0,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 30,
        source: 5,
        target: 8,
        consensus: 0,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 31,
        source: 8,
        target: 10,
        consensus: 0,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 32,
        source: 10,
        target: 11,
        consensus: 0,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 33,
        source: 11,
        target: 12,
        consensus: 0,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },
//    consensus 1
     {
      data: {
        id: 34,
        source: 1,
        target: 3,
        consensus: 1,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 35,
        source: 3,
        target: 4,
        consensus: 1,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 36,
        source: 4,
        target: 6,
        consensus: 1,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 37,
        source: 6,
        target: 8,
        consensus: 1,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 38,
        source: 8,
        target: 10,
        consensus: 1,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 39,
        source: 10,
        target: 11,
        consensus: 1,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 40,
        source: 11,
        target: 12,
        consensus: 1,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },
//    consensus 2
    {
      data: {
        id: 41,
        source: 1,
        target: 3,
        consensus: 2,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },
    {
      data: {
        id: 42,
        source: 3,
        target: 4,
        consensus: 2,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 43,
        source: 4,
        target: 7,
        consensus: 2,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 44,
        source: 7,
        target: 9,
        consensus: 2,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    },{
      data: {
        id: 45,
        source: 9,
        target: 10,
        consensus: 2,
        weight: 0.33,
        level: 1
      },
      classes: 'edge consensus'
    }
  ]
};
