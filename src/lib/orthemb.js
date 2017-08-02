import DiGraph from "./digraph.js";

class Face {
    constructor (link, arcs, exterior=false) {
        this.link = link;
        this.arcs = arcs;
        this.exterior = exterior;

        // this.edges is array of [arc1index, arc2index] around the face
        this.edges = this.arcs.map(a => this.link.edges[a.edge].map(b => b.index));

        this.turns = this.arcs.map(a => 1);
    }

    edgeOfIntersection(other) {
        return this.edges.filter(e => other.edges.some(oe => oe[0]==e[0] && oe[1]==e[1]));
    }

    sourceCapacity() {
        return (this.exterior ? 0 :
                Math.max(4 - this.arcs.length, 0));
    }

    sinkCapacity() {
        return (this.exterior ? this.arcs.length + 4 :
                Math.max(this.arcs.length - 4, 0));
    }

    bend(arc, turns) {
        let i = this.arcs.indexOf(arc);
        turns.reverse();

        for (let t of turns) {
            let nArc = this.link.verts[arc.vert][(arc.vertpos-1)%2];
            let oArc = this.link.edges[nArc.edge][(nArc.edgepos+1)%2];

            this.arcs.splice(i, 0, oArc);
            this.turns.splice(i, 0, t);
        }
    }

    *iterateFrom(arc) {
        let ai = this.arcs.indexOf(arc);
        for (let i = ai; i < this.arcs.length; i++) {
            yield [this.arcs[i], this.turns[i]];
        }
        for (let i = 0; i < ai; i++) {
            yield [this.arcs[i], this.turns[i]];
        }
    }

    orientEdges(arc, orientation) {
        const dirs = ["left", "up", "right", "down"];
        let dir = dirs.indexOf(orientation);
        let ans = new Map();
        for (let [a, t] in this.iterateFrom(arc)) {
            ans.set(a.index, dirs[dir]);
            ans.set(this.link.edgeOpposite(a).index, dirs[(dir+2)%4]);
            dir = (dir+t+4)%4;
        }
        return ans;
    }

}

export default class OrthogonalDiagramEmbedding {
    constructor (shadow) {
        this.shadow = shadow.copy();

        this.faces = shadow.faces.map(f => new Face(shadow, f));
        let F = this.faces.reduce(
            (bigF, f) => (bigF.length >= f.length) ? bigF : f);
        F.exterior = true;

        this.faceNetwork = this.flowNetwork();
        this.bend();
        this.orientEdges();

        this.edges = this.faces.reduce((all, f) => all.concat(f.arcs), []);
        this.repairComponents();
    }

    flowNetwork() {
        let G = new DiGraph();

        let sourceDemand = this.faces.reduce(
            (tot, f) => (tot + f.sourceCapacity()), 0);
        G.addNode('s', {demand: -sourceDemand});
        for (let fi = 0; fi < this.faces.length; fi++) {
            let sourceCapacity = this.faces[fi].sourceCapacity();
            if (sourceCapacity > 0) {
                G.addEdge('s', fi, {weight: 0, capacity: sourceCapacity});
            }
        }

        let sinkDemand = this.faces.reduce(
            (tot, f) => (tot + f.sinkCapacity()), 0);
        G.addNode('t', {demand: sinkDemand});
        for (let fi = 0; fi < this.faces.length; fi++) {
            let sinkCapacity = this.faces[fi].sinkCapacity();
            if (sinkCapacity > 0) {
                G.addEdge(fi, 't', {weight: 0, capacity: sinkCapacity});
            }
        }

        for (let ai = 0; ai < this.faces.length; ai++) {
            for (let bi = 0; bi < this.faces.length; bi++) {
                if (ai != bi && this.faces[ai].edgeOfIntersection(this.faces[bi])) {
                    G.addEdge(ai, bi, {weight: 1});
                }
            }
        }

        console.log(G);
        return G;
    }

    bend() {
        let flow = this.faceNetwork.minCostFlow()[1];
        console.log(flow);

        for (let [a, flows] of flow) {
            for (let [b, w_a] of flows) {
                //console.log(a, b, w_a);
                if (!w_a || a == 's' || b == 's' || a == 't' || b == 't') {
                    continue;
                }

                let w_b = flow.get(b).get(a);

                let [A, B] = [this.faces[a], this.faces[b]];
                console.log(a, b, A, B);
                let [e_a, e_b] = A.edgeOfIntersection(B);

                let turnsA = (new Array(w_a).map(x => 1)).concat(
                    (new Array(w_b).map(x => -1)));
                let turnsB = (new Array(w_b).map(x => 1)).concat(
                    (new Array(w_a).map(x => -1)));

                this.subdivideEdge(e_a, turnsA.length);

                A.bend(e_a, turnsA);
                B.bend(e_b, turnsB);
            }
        }
    }

    repairComponents() {
        this.arcComponentMap = new Map();
        this.orderedArcs = [];

        for (let i = 0; i < this.shadow.components.length; i++) {
            for (let arc of this.shadow.components[i]) {
                this.arcComponentMap.set(arc.index, i);
                this.orderedArcs.push(arc);
            }
        }
    }

    orientEdges() {
        let orientations = new Map();
        orientations.set(this.faces[0].arcs[0].index, 'right');

        let G = this.faceNetwork.copy();
        G.removeNode('s');
        G.removeNode('t');

        for (let node of G.depthFirstSearch(0)) {
            let F = this.faces[node];
            for (let arc of F.arcs) {
                if (orientations.has(arc.index)) {
                    let newOrientations = F.orientEdges(
                        arc, orientations[arc.index]);
                    for (let [a, dir] of newOrientations) {
                        if (orientations.has(a.index)) {
                            console.assert(orientations[a.index] == dir, orientations, a.index, dir);
                        } else {
                            orientations[a.index] = dir;
                        }
                    }
                    break;
                }
            }
        }

        //console.assert(orientations.length, orientations);
        this.orientations = orientations;
    }

    orthogonalSpec() {
        return ['right', 'up'].map(
            dir => this.arcs.filter(a => orientations[a.index] == dir).map(
                a => [a.vert, a.opposite.vert] // TODO
            ));
    }

    orthogonalRep() {
        
    }
}
