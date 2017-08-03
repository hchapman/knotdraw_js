import DiGraph from "./digraph.js";
import OrthogonalRep from "./orthrep.js";

class Face {
    constructor (link, arcs, exterior=false) {
        this.link = link;
        this.arcs = arcs;
        this.exterior = exterior;

        // this.edges is array of [arc1index, arc2index] around the face
        // this.edges = this.arcs.map(a => this.link.edges[a.edge].map(b => b.index));
        this.edges = this.arcs.map(a => a.edge);
        this.edgeMap = []; for (let arc of this.arcs) { this.edgeMap[arc.edge] = arc; }

        this.turns = this.arcs.map(a => 1);
    }

    edgeOfIntersection(other) {
        let commonEdges = this.edges.filter(e => other.edges.some(oe => e == oe));
        if (commonEdges.length > 0) {
            let e = commonEdges.pop();
            return [this.edgeMap[e], other.edgeMap[e]];
        }
        return null;
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
        //turns.reverse();

        for (let t of turns) {
            arc = arc.edgeOpposite.vertNext;
            //arc = arc.vertPrev.edgeOpposite;
            this.arcs.splice(i+1, 0, arc);
            this.turns.splice(i+1, 0, t);
            i += 1;
        }
    }

    *iterateFrom(arc) {
        let ai = this.arcs.indexOf(arc);
        for (let i = ai+1; i < this.arcs.length; i++) {
            yield [this.arcs[i], this.turns[i]];
        }
        for (let i = 0; i < ai+1; i++) {
            yield [this.arcs[i], this.turns[i]];
        }
    }

    orientEdges(arc, orientation) {
        const dirs = ["left", "up", "right", "down"];
        console.assert(this.isValid(), this);
        let dir = dirs.indexOf(orientation);
        let ans = new Map();
        for (let [a, t] of this.iterateFrom(arc)) {
            dir = (dir+t+4)%4;
            ans.set(a.index, dirs[dir]);
            ans.set(a.edgeOpposite.index, dirs[(dir+2)%4]);
        }
        return ans;
    }

    isValid() {
        let face = [];
        let start = this.arcs[0];
        let arc = start;
        while (!face.includes(arc.index)) {
            face.push(arc.index);
            arc = arc.edgeOpposite.vertNext;
        }

        return this.arcs.reduce((val, a, i) => val && (a.index == face[i]), true);
    }

}

export default class OrthogonalDiagramEmbedding {
    constructor (shadow) {
        this.shadow = shadow.copy();

        this.faces = this.shadow.faces.map(f => new Face(this.shadow, f));
        let F = this.faces.reduce(
            (bigF, f) => (bigF.arcs.length >= f.arcs.length) ? bigF : f);
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
                if (ai != bi &&
                    this.faces[ai].edgeOfIntersection(this.faces[bi]) != null) {
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
                let [arc_ai, arc_bi] = A.edgeOfIntersection(B);

                let arc_a = this.shadow.arcs[arc_ai];
                let arc_b = this.shadow.arcs[arc_bi];

                let turnsA = (new Array(w_a).fill(1)).concat(
                    (new Array(w_b).fill(-1)));
                let turnsB = (new Array(w_b).fill(1)).concat(
                    (new Array(w_a).fill(-1)));

                this.subdivideEdge(arc_a, turnsA.length);

                A.bend(arc_a, turnsA);
                B.bend(arc_b, turnsB);
            }
        }
    }

    subdivideEdge(arc, n) {
        let head = arc;
        let tail;
        let backwards = !this.shadow.components.some(c => c.includes(arc));
        if (backwards) {
            tail = head;
            head = tail.edgeOpposite;
        } else {
            tail = head.edgeOpposite;
        }

        let strands = (new Array(2*n).fill(0)).map(i => this.shadow.newArc());
        strands.splice(0, 0, head);
        strands.push(tail);

        // Glue edges
        for (let si = 0; si < strands.length; si += 2) {
            strands[si].edgeOpposite = strands[si+1];
            strands[si+1].edgeOpposite = strands[si];

            this.shadow.setEdge(this.shadow.edges.length,
                                [strands[si].index, strands[si+1].index]);
        }

        // Glue degree 2 joints
        for (let si = 1; si < strands.length-1; si += 2) {
            strands[si].vertNext = strands[si+1];
            strands[si].vertPrev = strands[si+1];

            strands[si+1].vertNext = strands[si];
            strands[si+1].vertPrev = strands[si];

            this.shadow.setVert(this.shadow.verts.length,
                                [strands[si].index, undefined, strands[si+1].index, undefined]);
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
        let orientations = [];
        orientations[this.faces[0].arcs[0].index] = 'right';

        let G = this.faceNetwork.copy();
        G.removeNode('s');
        G.removeNode('t');

        for (let node of G.depthFirstSearch(0)) {
            let F = this.faces[node];
            for (let arc of F.arcs) {
                if (arc.index in orientations) {
                    let newOrientations = F.orientEdges(
                        arc, orientations[arc.index]);
                    for (let [a, dir] of newOrientations) {
                        if (a in orientations) {
                            console.assert(orientations[a] == dir, orientations, a, dir);
                        } else {
                            orientations[a] = dir;
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
            dir => this.edges.filter(a => this.orientations[a.index] == dir).map(
                a => [a.vert, a.edgeOpposite.vert] // TODO
            ));
    }

    orthogonalRep() {
        return new OrthogonalRep(...this.orthogonalSpec());
    }
}
