importScripts('lalolib/lalolib.js');
importScripts('lodash.min.js');

importScripts('planarmap-em/libplanarmap-em.js');

import LinkShadow from "./lib/shadow.js";
import OrthogonalDiagramEmbedding from "./lib/orthemb.js";
import ForceLinkDiagram from "./lib/forcediagram.js";

function sleep(millis)
{
    var date = new Date();
    var curDate = null;
    do { curDate = new Date(); }
    while(curDate-date < millis);
}

self.randomDiagram = Module.cwrap('randomDiagram', 'number',
                                  ['number', 'number', 'number', 'number', 'number', 'number']);

function randomDiagram(n_verts, n_comps, max_att, type) {
    /* Generate a random diagram using planarmap-emscripten */
    
    let vertPtr = Module._malloc(4);

    let nVerts = self.randomDiagram(n_verts, n_comps, max_att, type,
                                    Math.random()*2**32, vertPtr);
    if (nVerts == 0) {
        // No diagram was successfully sampled
        Module._free(vertPtr);
        return null;
    }

    let vertArray = Module.getValue(vertPtr, "i32*");

    let view = Module.HEAP32.subarray(vertArray/4, vertArray/4+(4*nVerts));

    let pd = [];
    for (let vi = 0; vi < nVerts; vi++) {
        let vert = [];
        for (let pos = 0; pos < 4; pos++) {
            vert.push(view[vi*4+pos]);
        }
        pd.push(vert);
    }

    Module._free(vertArray);
    Module._free(vertPtr);

    return pd;
}

var workerFunctions = {
    setRandomLinkDiagram: function(n_verts, n_comps, max_att, type) {
        let sigma = randomDiagram(n_verts, n_comps, max_att, type);
        if (sigma == null) {
            postMessage({
                function: "raiseError",
                arguments: ["A diagram was not sampled successfully"]
            });
            return;
        }

        postMessage({
            function: "updatePDCode",
            arguments: ["["+sigma.map(v => "["+v+"]")+"]"]
        });
        workerFunctions.setLinkDiagram(sigma);
    },

    setLinkDiagram: function(sigma, crossBend) {
        self.shadow = new LinkShadow(sigma);
        self.orthShadow = new OrthogonalDiagramEmbedding(self.shadow);

        let rep = self.orthShadow.orthogonalRep();

        let spec = self.orthShadow.orthogonalSpec();

        let gridEmb = rep.basicGridEmbedding();

        let verts = [];
        for (let [i, vert] of gridEmb) {
            verts[i] = vert;
        }
        
        let edges = [];
        for (let [source, sink, data] of rep.graph.edgeGen()) {
            //console.log(data);
            if (!rep.dummy.has(data)) {
                edges.push([source, sink]);
            }
        }
        let faces = rep.faces.map(f => f.evPairs.map(ev => ev[1]));
        //faces.forEach(f => f.reverse());

        self.components = self.orthShadow.components.map(c => c.map(a => a.vert));
        self.faces = self.orthShadow.faces.map(f => {
            let face = f.arcs.map(a => a.vert);
            face.exterior = f.exterior;
            face.degree = f.degree;
            return face;
        });
        self.force_shadow = new ForceLinkDiagram(verts, edges, self.faces, self.components);
        
        //console.log(faces);
        

        postMessage({
            function: "setLinkDiagram",
            arguments: [self.force_shadow]
        });

        workerFunctions.embedDiagram();
    },

    embedDiagram: function() {
        let tstart = Date.now();

        let thresh = 5e-10;

        let curDate;
        self.n_steps = 100;//50;
        let max_steps = 100;

        for (let i = 0; i < max_steps; i++) {
            let procStart = Date.now();
            console.assert(!isNaN(self.force_shadow.aExp));
            postMessage({
                function: "setLinkDiagram",
                arguments: [self.force_shadow]
            });

            self.force_shadow.update();
            self.force_shadow.aExp -= (1 - 0.4)/self.n_steps;
            self.force_shadow.erExp += (4 - 2)/self.n_steps;
            self.force_shadow.dbar -= (3*self.force_shadow.delta)/self.n_steps;

            //do { curDate = Date.now(); }
            //while( curDate-procStart < 50);
        }

        postMessage({
            function: "finalizeLinkDiagram",
            arguments: [self.force_shadow, self.force_shadow.faces, self.components]
        });
    },

    stepUpdate: function() {
        self.force_shadow.update();
        self.force_shadow.aExp -= (1 - 0.4)/self.n_steps;
        self.force_shadow.reExp += (4 - 2)/self.n_steps;
        self.force_shadow.dbar -= (3*self.force_shadow.delta)/self.n_steps;

        postMessage({
            function: "finalizeLinkDiagram",
            arguments: [self.force_shadow, self.force_shadow.faces, self.components]
        });
    }
}

onmessage = function(e) {
    workerFunctions[e.data.function](...e.data.arguments);
}
