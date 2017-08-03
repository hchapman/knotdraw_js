import DiGraph from "./digraph.js";

class OrthogonalFace {
    constructor(graph, edgeAndVert) {
        let [edge, vertex] = edgeAndVert;
        this.evPairs = [];
        this.evPairs.push(edgeAndVert);

        while (true) {
            [edge, vertex] = this.evPairs[this.evPairs.length-1];
            edge = graph.nextEdgeAtVertex(edge, vertex);
            vertex = (vertex === edge.source) ? edge.sink : edge.source;
            if (this.evPairs[0][0] == edge && this.evPairs[0][1] == vertex) {
                break;
            } else {
                this.evPairs.push([edge, vertex]);
            }
        }

        this.addTurns();
    }

    addTurns() {
        this.turns = [];
        for (let i = 0; i < this.evPairs.length; i++) {
            let [e0, v0] = this.evPairs[i];
            let [e1, v1] = this.evPairs[(i+1)%this.evPairs.length];

            if (e0.kind == e1.kind) {
                this.turns.push(0);
            } else {
                let t = (e0.tail == v0) ^ (e1.head == v0) ^ (e0.kind == 'horizontal');
                this.turns.push(t ? -1 : 1);
            }
        }

        let rotation = this.turns.reduce((s, x) => s+x);
        console.assert( Math.abs(rotation) == 4, rotation );
        this.exterior = (rotation == -4);
    }
}

export default class OrthogonalRep {
    constructor(horiz_pairs, vert_pairs) {
        this.graph = new DiGraph();
        this.nedges = 0;
        for (let [a, b] of horiz_pairs) {
            this.graph.addEdge(a, b, {kind: "horizontal", index: this.nedges++});
        }
        for (let [a, b] of vert_pairs) {
            this.graph.addEdge(a, b, {kind: "vertical", index: this.nedges++});
        }

        this.buildFaces();
        this.makeTurnRegular();
    }

    buildFaces() {
        this.faces = [];
        let unseen = new Map();
        for (let [source, sink, data] of this.graph.edgeGen()) {
            if (!unseen.has(data.index)) { unseen.set(data.index, new Map()); }
            unseen.get(data.index).set(source, data);
            unseen.get(data.index).set(sink, data);
        }

        while (unseen.size > 0) {
            let [vert, edge] = unseen.values().next().value.entries().next().value;
            unseen.get(edge.index).delete(vert);
            if (unseen.get(edge.index).size <= 0) { unseen.delete(edge.index); }

            let face = OrthogonalFace(this, es);
            face.evPairs.forEach(
                x => {
                    let [edge, vert] = x;
                    unseen.get(edge.index).delete(vert);
                    if (unseen.get(edge.index).size <= 0) {
                        unseen.delete(edge.index);
                    }
                });
            this.faces.push(face);
        }
    }

    makeTurnRegular() {
        let dummy = new Set();
        let regular = self.faces.filter(F => F.isTurnRegular());
        let irregular = self.faces.filter(F => !F.isTurnRegular());

        let i = 0;
        while (irregular.length > 0) {
            let F = irregular.pop();
            let [i, j] = F.kittyCorner();
            let [v0, v1] = [F.evPairs[i][1], F.evPairs[j][1]];

            let kind = ['vertical', 'horizontal'][i%2];
            let e;
            if (this.graph.incoming(v0).some(e => e.kind == kind)) {
                e = this.graph.addEdge(v0, v1, {kind: kind, index: this.nedges++});
            } else {
                e = this.graph.addEdge(v1, v0, {kind: kind, index: this.nedges++});
            }
            dummy.add(e);

            for (let v of [v0, v1]) {
                F = OrthogonalFace(this, (e, v));
                if (F.isTurnRegular()) {
                    regular.push(F);
                } else {
                    irregular.push(F);
                }
            }

            i++;
        }

        [this.faces, this.dummy] = [regular, dummy];
    }

    DagFromDirection(kind) {
        let H = new DiGraph();
        this.graph.nodes.forEach((data, v) =>
                                 H.addNode(v));
        for (let [source, sink, data] of this.graph.edgeGen()) {
            
        }
    }

    chainCoordinates(kind) {
        let D = this.DagFromDirection(kind);
        let chainCoords = topologicalNumbering(D);

        let coords = new Map();
        this.graph.nodes.forEach((data, v) =>
                                 coords.set(v, chainCoords[D.vertexToChain(v)]));
    }

    basicGridEmbedding() {
        let V = this.chainCoordinates('horizontal');
        let H = this.chainCoordinates('vertical');

        let emb = new Map();
        this.graph.nodes.forEach((data, v) => emb.set(v, [H[v], V[v]]));
        return emb;
    }
}
