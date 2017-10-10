var PoaGraph_elements = {
  nodes: [
    {
      data: {
        id: 0,
        nucleobase: 'T',
        source: [1],
        weight: 0.5,
      },
      position: { x: 10, y: 0 },
      // classes: 'czerwony zielony'
    }, {
      data: {
        id: 1,
        nucleobase: 'H',
        source: [1],
        weight: 0.5,
      },
      position: { x: 20, y: 100 },
    },{
      data: {
        id: 2,
        nucleobase: 'P',
        source: [0],
        weight: 0.5
      },
      position: { x: 30, y: 100 },
    },{
      data: {
        id: 3,
        nucleobase: 'K',
        source: [0, 1],
        weight: 1
      },
      position: { x: 40, y: 100 },
    },
    {
      data: {
        id: 4,
        nucleobase: 'M',
        source: [0, 1],
        weight: 1
      },
      position: { x: 50, y: 0 },
    },{
      data: {
        id: 5,
        nucleobase: 'I',
        source: [0],
        weight: 0.5
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 6,
        nucleobase: 'L',
        source: [1],
        weight: 0.5
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 7,
        nucleobase: 'V',
        source: [0, 1],
        weight: 1
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 8,
        nucleobase: 'R',
        source: [0, 1],
        weight: 1
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 9,
        nucleobase: 'P',
        source: [0],
        weight: 0.5
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 10,
        nucleobase: 'Q',
        source: [0],
        weight: 0.5
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 11,
        nucleobase: 'K',
        source: [0],
        weight: 0.5
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 12,
        nucleobase: 'N',
        source: [0, 1],
        weight: 1
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 13,
        nucleobase: 'E',
        source: [0, 1],
        weight: 1
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 14,
        nucleobase: 'T',
        source: [0, 1],
        weight: 1
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 15,
        nucleobase: 'V',
        source: [0],
        weight: 0.5
      },
      position: { x: 60, y: 100 },
    },{
      data: {
        id: 16,
        nucleobase: 'I',
        source: [ 1],
        weight: 0.5
      }
    },{
      data: {
        id: 17,
        nucleobase: 'M',
        source: [ 1],
        weight: 0.5
      },
      position: { x: 60, y: 100 },
    },
  ],
  edges: [
    {
      data: {
        id: 18,
        source: 0,
        target: 1,
        consensus: -1,
        weight: 0.5
      }
    },{
      data: {
        id: 19,
        source: 1,
        target: 3,
        consensus: -1,
        weight: 0.5
      }
    },{
      data: {
        id: 20,
        source: 2,
        target: 3,
        consensus: -1,
        weight: 0.5
      }
    },{
      data: {
        id: 21,
        source: 3,
        target: 4,
        consensus: 1,
        weight: 1
      }
    },{
      data: {
        id: 22,
        source: 4,
        target: 5,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 23,
        source: 4,
        target: 6,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 24,
        source: 5,
        target: 7,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 25,
        source: 6,
        target: 7,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 26,
        source: 7,
        target: 8,
        consensus: -1,
        weight: 1
      }
    }, {
      data: {
        id: 27,
        source: 8,
        target: 9,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 28,
        source: 8,
        target: 12,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 29,
        source: 9,
        target: 10,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 30,
        source: 10,
        target: 11,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 31,
        source: 11,
        target: 12,
        consensus: 0,
        weight: 0.5
      }
    }, {
      data: {
        id: 32,
        source: 12,
        target: 13,
        consensus: 1,
        weight: 1
      }
    }, {
      data: {
        id: 33,
        source: 13,
        target: 14,
        consensus: 2,
        weight: 1
      }
    }, {
      data: {
        id: 34,
        source: 14,
        target: 15,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 35,
        source: 14,
        target: 16,
        consensus: -1,
        weight: 0.5
      }
    }, {
      data: {
        id: 36,
        source: 16,
        target: 17,
        consensus: -1,
        weight: 0.5
      }
    }
  ]
};
