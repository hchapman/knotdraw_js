import DiGraph from "./digraph.js";

function partialSums(L) {
    let ans = [0];
    for (let x of L) {
        ans.push(ans[ans.length-1] + x);
    }
    return ans;
}

function elementMap(part) {
    let ans = [];
    for (let p of part) {
        for (let x of p) {
            ans[x] = p;
        }
    }
    return ans;
}

function basicTopologicalNumbering(G) {
    let inValences = new Map();
    G.nodes.forEach(
        (data, v) => inValences.set(v, G.indegree(v)));

    let numbering = [];

    let currSources = [];
    for (let [v, deg] of inValences) {
        if (deg == 0) { currSources.push(v); }
    }

    let currNumber = 0;

    while (inValences.size > 0) {
        let newSources = [];
        for (let v of currSources) {
            inValences.delete(v);
            numbering[v] = currNumber;

            for (let e of G.outgoing(v)) {
                let w = e.sink;
                inValences.set(w, inValences.get(w) - 1);
                if (inValences.get(w) == 0) {
                    newSources.push(w);
                }
            }
            currSources = newSources;
        }
        currNumber += 1;
    }

    return numbering;
}

function topologicalNumbering(G) {
    let n = basicTopologicalNumbering(G);
    let success = true;

    while (success) {
        success = false;
        for (let [v, d] of G.nodes) {
            let below = G.incoming(v).filter(e => e.dummy == false).length;
            let above = G.outgoing(v).filter(e => e.dummy == false).length;

            if (above != below) {
                let newPos;
                if (above > below) {
                    newPos = Math.min(...G.outgoing(v).map(e => n[e.sink])) - 1;
                } else {
                    newPos = Math.max(...G.incoming(v).map(e => n[e.source])) + 1;
                }

                if (newPos != n[v]) {
                    n[v] = newPos;
                    success = true;
                }
            }
        }
    }

    return n;
}

function kittyCorner(turns) {
    let rotations = partialSums(turns);
    //let reflexCorners = turns.filter(t => t == -1).map((t, i) => i);
    let reflexCorners = turns.map((t, i) => i).filter(i => turns[i] == -1);

    for (let r0 of reflexCorners) {
        for (let r1 of reflexCorners.filter(r => r > r0)) {
            if (rotations[r1] - rotations[r0] == 2) {
                return [r0, r1];
            }
        }
    }
    return null;
}

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
                let t = (e0.source == v0) ^ (e1.sink == v0) ^ (e0.kind == 'horizontal');
                this.turns.push(t ? -1 : 1);
            }
        }

        let rotation = this.turns.reduce((s, x) => s+x);
        console.assert( Math.abs(rotation) == 4, rotation );
        this.exterior = (rotation == -4);
    }

    kittyCorner() {
        if (!this.exterior) {
            return kittyCorner(this.turns);
        }
        return null;
    }

    isTurnRegular() {
        return this.kittyCorner() === null;
    }

    switches(swap) {
        function edgeToEndpoints(e) {
            if (swap && e.kind == 'horizontal') {
                return [e.sink, e.source];
            }
            return [e.source, e.sink];
        }

        let ans = [];
        for (let i = 0; i < this.evPairs.length; i++) {
            let [e0, v0] = this.evPairs[i];
            let [t0, h0] = edgeToEndpoints(e0);
            let [t1, h1] = edgeToEndpoints(this.evPairs[(i+1)%this.evPairs.length][0]);

            if (t0 == t1 && t1 == v0) {
                ans.push({index: i, kind: 'source', turn: this.turns[i]});
            } else if (h0 == h1 && h1 == v0) {
                ans.push({index: i, kind: 'sink', turn: this.turns[i]});
            }
        }

        return ans;
    }

    saturationEdges(swap) {
        function saturateFace(faceInfo) {
            for (let i = 0; i < faceInfo.length; i++) {
                if (faceInfo[i].turn == -1) {
                    faceInfo = faceInfo.slice(i).concat(faceInfo.slice(0, i));
                    break;
                }
            }

            for (let i = 0; i < faceInfo.length-2; i++) {
                let [x, y, z] = faceInfo.slice(i, i+3);
                if (x.turn == -1 && y.turn == z.turn && z.turn == 1) {
                    let [a, b] = x.kind == 'sink' ? [x, z] : [z, x];
                    let remaining = faceInfo.slice(0, i)
                        .concat([{index: z.index, kind: z.kind, turn: 1}])
                        .concat(faceInfo.slice(i+3));
                    return [[a.index, b.index]].concat(saturateFace(remaining));
                }
            }
            return [];
        }

        let newEdges = saturateFace(this.switches(swap));
        return newEdges.map(([a,b]) => [this.evPairs[a][1], this.evPairs[b][1]]);
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

            let face = new OrthogonalFace(this, [edge, vert]);
            face.evPairs.forEach(
                x => {
                    let [edge, vert] = x;
                    if (unseen.has(edge.index)) {
                        unseen.get(edge.index).delete(vert);
                        if (unseen.get(edge.index).size <= 0) {
                            unseen.delete(edge.index);
                        }
                    }
                });
            this.faces.push(face);
        }
    }

    link(vert) {
        let link = this.graph.incidentEdges(vert);
        function score(e) {
            return (e.source == vert)*2 + (e.kind == 'vertical');
        }
        link.sort((a, b) => score(a) > score(b));
        return link;
    }

    nextEdgeAtVertex(edge, vert) {
        let link = this.link(vert);
        return link[ (link.indexOf(edge)+1) % link.length ];
    }

    makeTurnRegular() {
        let dummy = new Set();
        let regular = this.faces.filter(F => F.isTurnRegular());
        let irregular = this.faces.filter(F => !F.isTurnRegular());

        let ell = 0;
        while (irregular.length > 0) {
            let F = irregular.pop();
            let [i, j] = F.kittyCorner();
            let [v0, v1] = [F.evPairs[i][1], F.evPairs[j][1]];

            let kind = ['vertical', 'horizontal'][ell%2];
            let e;
            if (this.graph.incoming(v0).some(e => e.kind == kind)) {
                e = this.graph.addEdge(v0, v1, {kind: kind, index: this.nedges++});
            } else {
                e = this.graph.addEdge(v1, v0, {kind: kind, index: this.nedges++});
            }
            dummy.add(e);

            for (let v of [v0, v1]) {
                F = new OrthogonalFace(this, [e, v]);
                if (F.isTurnRegular()) {
                    regular.push(F);
                } else {
                    irregular.push(F);
                }
            }

            ell++;
        }

        [this.faces, this.dummy] = [regular, dummy];
    }

    saturationEdges(swap) {
        return this.faces.reduce(
            (arr, f) => arr.concat(f.saturationEdges(swap)), []);
    }

    dagFromDirection(kind) {
        let H = new DiGraph();
        this.graph.nodes.forEach((data, v) =>
                                 H.addNode(v));
        for (let [source, sink, data] of this.graph.edgeGen()) {
            if (data.kind == kind) { H.addEdge(source, sink); }
        }

        let maximalChains = H.weakComponents();
        let vertexToChain = elementMap(maximalChains);

        let D = new DiGraph();
        maximalChains.forEach((c, i) => D.addNode(c));

        for (let [source, sink, data] of this.graph.edgeGen()) {
            if (data.kind != kind) {
                let d = D.addEdge(vertexToChain[source], vertexToChain[sink]);
                d.dummy = (this.dummy.has(data));
            }
        }

        for (let [u, v] of this.saturationEdges(false)) {
            let d = D.addEdge(vertexToChain[u], vertexToChain[v]);
            d.dummy = true;
        }
        for (let [u, v] of this.saturationEdges(true)) {
            if (kind == 'vertical') {
                let t = u;
                u = v;
                v = t;
            }

            let d = D.addEdge(vertexToChain[u], vertexToChain[v]);
            d.dummy = true;
        }

        D.vertexToChain = vertexToChain;
        return D;
    }

    chainCoordinates(kind) {
        let D = this.dagFromDirection(kind);
        let chainCoords = topologicalNumbering(D);

        let coords = new Map();
        this.graph.nodes.forEach((data, v) =>
                                 coords.set(v, chainCoords[D.vertexToChain[v]]));
        return coords;
    }

    basicGridEmbedding() {
        let V = this.chainCoordinates('horizontal');
        let H = this.chainCoordinates('vertical');

        let emb = new Map();
        this.graph.nodes.forEach((data, v) => emb.set(v, [H.get(v), V.get(v)]));
        return emb;
    }
}
