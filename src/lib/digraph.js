export default class DiGraph {
    /* Lifted from networkx's network_simplex routines */

    constructor() {
        this.nodes = new Map();
        this.edges = new Map();

        this.redges = new Map();
    }

    copy() {
        let G = new DiGraph();
        G.nodes = new Map(this.nodes);
        for (let [u, sinks] of this.edges) {
            G.edges.set(u, new Map(this.edges.get(u)));
        }
        for (let [v, sources] of this.redges) {
            G.redges.set(v, new Map(this.redges.get(v)));
        }
        return G;
    }

    addNode(name, attrs=undefined) {
        if (attrs === undefined) { attrs = {}; }
        if (this.nodes.has(name)) {
            for (let [k, v] of Object.entries(attrs)) {
                this.nodes.get(name)[k] = v;
            }
        } else {
            this.nodes.set(name, {name: name});
            for (let [k, v] of Object.entries(attrs)) {
                this.nodes.get(name)[k] = v;
            }
        }
        if (!this.edges.has(name)) {
            this.edges.set(name, new Map());
        }
        if (!this.redges.has(name)) {
            this.redges.set(name, new Map());
        }
    }

    removeNode(name) {
        this.nodes.delete(name);
        this.edges.delete(name);
        for (let [source, sinks] of this.edges) {
            sinks.delete(name);
        }
        for (let [sink, sources] of this.redges) {
            sources.delete(name);
        }
    }

    addEdge(source, sink, attrs=undefined) {
        if (!this.nodes.has(source)) {
            this.addNode(source);
        }
        if (!this.nodes.has(sink)) {
            this.addNode(sink);
        }
        this.edges.get(source).set(sink, {source: source, sink: sink});

        let edge = this.edges.get(source).get(sink);
        this.redges.get(sink).set(source, edge);

        if (attrs !== undefined) {
            for (let [k, v] of Object.entries(attrs)) {
                edge[k] = v;
            }
        }

        return edge;
    }

    removeEdge(source, sink) {
        if (this.edges.has(source)) {
            this.edges.get(source).delete(sink);
        }
        if (this.redges.has(sink)) {
            this.redges.get(sink).delete(source);
        }
    }

    getEdge(source, sink) {
        return this.edges.get(source).get(sink);
    }

    hasEdge(source, sink) {
        if (!this.edges.has(source)) {
            return false;
        }
        if (!this.edges.get(source).has(sink)) {
            return false;
        }
        return true;
    }

    *edgeGen() {
        for (let [source, sinks] of this.edges) {
            for (let [sink, d] of sinks) {
                yield [source, sink, d];
            }
        }
    }

    getReverseEdges() {
        return this.redges;
    }

    incidentEdges(vert) {
        let ans = [];
        for (let [sink, data] of this.edges.get(vert)) {
            ans.push(data);
        }
        for (let [source, data] of this.redges.get(vert)) {
            ans.push(data);
        }
        return ans;
    }

    incoming(vert) {
        return Array.from(this.redges.get(vert).values());
    }
    outgoing(vert) {
        return Array.from(this.edges.get(vert).values());
    }

    outdegree(vert) {
        return this.outgoing(vert).length;
    }
    indegree(vert) {
        return this.incoming(vert).length;
    }

    treePath(start, stop) {
        // Requires that this is a tree to avoid real programming...

        if (start == stop) { return [start]; }

        // this.edges is a map source -> sink
        let rev_edges = this.getReverseEdges();

        let search_nodes = new Set();
        for (let u of this.nodes.keys()) {
            search_nodes.add(u);
        }

        // DFS for now, since it's easier
        let path;
        let search_stack = [[[start], start]];

        do {
            let [path, node] = search_stack.pop();
            search_nodes.delete(node);
            if (this.edges.has(node)) {
                for (let [v, d] of this.edges.get(node)) {
                    if (v == stop) {
                        return path.concat([v]);
                    }
                    if (search_nodes.has(v)) {
                        search_stack.push([path.concat([v]), v]);
                    }
                }
            }
            if (rev_edges.has(node)) {
                for (let [v, d] of rev_edges.get(node)) {
                    if (v == stop) {
                        return path.concat([v]);
                    }
                    if (search_nodes.has(v)) {
                        search_stack.push([path.concat([v]), v]);
                    }
                }
            }
        } while (search_nodes.size > 0);

        return path;
    }

    connectedComponentSubgraphs() {
        // Requires that this is a tree to avoid real programming...
        // CAVEAT: We only care about nodes so "meh" to edges

        // this.edges is a map source -> sink
        let rev_edges = this.getReverseEdges();

        let search_nodes = new Set();
        for (let u of this.nodes.keys()) {
            search_nodes.add(u);
        }

        // DFS for now, since it's easier
        let search_stack = [];
        let components = [];

        while (search_nodes.size > 0) {
            search_stack = [search_nodes.values().next().value];

            let G = new DiGraph();
            while (search_stack.length > 0) {
                let node = search_stack.pop();
                search_nodes.delete(node);
                G.addNode(node);
                if (this.edges.has(node)) {
                    for (let [v, d] of this.edges.get(node)) {
                        if (search_nodes.has(v)) {
                            search_stack.push(v);
                        }
                    }
                }
                if (rev_edges.has(node)) {
                    for (let [v, d] of rev_edges.get(node)) {
                        if (search_nodes.has(v)) {
                            search_stack.push(v);
                        }
                    }
                }
            }
            components.push(G);
        }

        return components;
    }

    weakComponents() {
        return this.connectedComponentSubgraphs().map(
            c => Array.from(c.nodes.keys())
        );
    }

    *depthFirstSearch(start) {
        let stack = [start];
        let seen = new Set(stack);

        let rev_edges = this.getReverseEdges();

        while (stack.length > 0) {
            let current = stack.pop();
            for (let v of this.edges.get(current).keys()) {
                if (!seen.has(v)) {
                    stack.push(v);
                    seen.add(v);
                }
            }
            for (let v of rev_edges.get(current).keys()) {
                if (!seen.has(v)) {
                    stack.push(v);
                    seen.add(v);
                }
            }

            yield current;
        }
    }

    minCostFlow() {
        let [H, T, y, artificialEdges, flowCost, r] = this._initialTreeSolution();

        let c = [];
        for (let [u, sinks] of H.edges) {
            c[u] = [];
            for (let [v, d] of sinks) {
                c[u][v] = _.get(d, 'weight', 0) + y[u] - y[v];
            }
        }

        let reverse;
        while (true) {
            let newEdge = H._findEnteringEdge(c);
            if (newEdge.length != 2) {
                break;
            }

            let cycleCost = Math.abs(c[newEdge[0]][newEdge[1]]);

            let path1 = T.treePath(r, newEdge[0]);
            let path2 = T.treePath(r, newEdge[1]);

            let join = r;

            for (let i = 1; i < path1.length; i++) {
                let node = path1[i];
                if (i < path2.length && node == path2[i]) {
                    join = node;
                } else {
                    break;
                }
            }

            path1 = path1.slice(path1.indexOf(join));
            path2 = path2.slice(path2.indexOf(join));
            let cycle = [];

            //if (_.get(H.getEdge(...newEdge), "flow", 0) == 0) {
            let hEdgeD = H.getEdge(...newEdge);
            if (hEdgeD.flow === 0 || hEdgeD.flow === undefined) {
                reverse = false;
                path2.reverse();
                cycle = path1.concat(path2);
            } else {
                reverse = true;
                path1.reverse();
                cycle = path2.concat(path1);
            }

            let [leavingEdge, eps] = H._findLeavingEdge(T, cycle, newEdge, reverse);
            console.assert(leavingEdge.length == 2);

            if (eps != 0) {
                flowCost -= cycleCost * eps;
                if (cycle.length == 3) {
                    if (reverse) {
                        eps = -eps;
                    }
                    let [u, v] = newEdge;
                    let uvEdge = H.getEdge(u, v);
                    let vuEdge = H.getEdge(v, u);

                    if (uvEdge.flow === undefined) {
                        uvEdge.flow = eps;
                    } else {
                        uvEdge.flow += eps;
                    }
                    if (vuEdge.flow === undefined) {
                        vuEdge.flow = eps;
                    } else {
                        vuEdge.flow += eps;
                    }
                    //_.set(H.getEdge(u, v), "flow", _.get(H.getEdge(u, v), "flow", 0) + eps);
                    //_.set(H.getEdge(v, u), "flow", _.get(H.getEdge(v, u), "flow", 0) + eps);
                } else {
                    for (let j = 0; j < cycle.length-1; j++) {
                        let u = cycle[j];
                        let v = cycle[j+1];
                        if ((u == newEdge[0] && v == newEdge[1]) ||
                            T.hasEdge(u, v)) {
                            //_.set(H.getEdge(u, v), "flow", _.get(H.getEdge(u, v), "flow", 0) + eps);
                            let uvEdge = H.getEdge(u, v);

                            if (uvEdge.flow === undefined) {
                                uvEdge.flow = eps;
                            } else {
                                uvEdge.flow += eps;
                            }
                        } else {
                            //_.set(H.getEdge(v, u), "flow", _.get(H.getEdge(v, u), "flow", 0) - eps);
                            let vuEdge = H.getEdge(v, u);

                            if (vuEdge.flow === undefined) {
                                vuEdge.flow = eps;
                            } else {
                                vuEdge.flow -= eps;
                            }
                        }
                    }
                }
            }

            T.addEdge(...newEdge);
            T.removeEdge(...leavingEdge);

            if (newEdge[0] != leavingEdge[0] || newEdge[1] != leavingEdge[1]) {
                let forest = T.copy();
                forest.removeEdge(...newEdge);

                let [R, notR] = forest.connectedComponentSubgraphs();

                if (notR.nodes.has(r)) {
                    let tmp = R;
                    R = notR;
                    notR = tmp;
                }

                if (R.nodes.has(newEdge[0])) {
                    for (let v of notR.nodes.keys()) {
                        y[v] += c[newEdge[0]][newEdge[1]];
                    }
                } else {
                    for (let v of notR.nodes.keys()) {
                        y[v] -= c[newEdge[0]][newEdge[1]];
                    }
                }

                for (let [u, sinks] of H.edges) {
                    for (let [v, d] of sinks) {
                        if (notR.nodes.has(u) || notR.nodes.has(v)) {
                            let uvEdge = H.getEdge(u, v);

                            if (uvEdge.weight === undefined) {
                                c[u][v] = y[u] - y[v];
                            } else {
                                c[u][v] = uvEdge.weight + y[u] - y[v];
                            }
                            //c[u][v] = _.get(H.getEdge(u, v), "weight", 0) + y[u] - y[v];
                        }
                    }
                }
            }
        }

        for (let [u, v] of artificialEdges) {
            H.removeEdge(u, v);
        }

        for (let u of H.nodes.keys()) {
            if (!this.nodes.has(u)) {
                H.removeNode(u);
            }
        }

        let flowDict = this._createFlowDict(H);

        return [flowCost, flowDict];
    }

    _createFlowDict(H) {
        let flowDict = new Map();

        for (let [u, sinks] of this.edges) {
            flowDict.set(u, new Map());
            for (let [v, d] of sinks) {
                if (H.hasEdge(u, v)) {
                    flowDict.get(u).set(v, _.get(H.getEdge(u, v), "flow", 0));
                } else {
                    flowDict.get(u).set(v, 0);
                }
            }
        }

        return flowDict;
    }

    _initialTreeSolution() {
        let H = new DiGraph();

        let maxWeight = 0;
        for (let [u, sinks] of this.edges) {
            for (let [v, d] of sinks) {
                if (_.get(d, 'capacity', 1) > 0) {
                    H.addEdge(u, v, d);
                    maxWeight = Math.max(maxWeight, _.get(d, 'weight', 0));
                }
            }
        }

        for (let [n, d] of this.nodes) {
            if (_.get(d, 'demand', 0) != 0) {
                H.addNode(n, d);
            }
        }

        let nodeIter = H.nodes.entries();
        let [r, _d] = nodeIter.next().value;

        let T = new DiGraph();
        let y = []; y[r] = 0;
        let artificialEdges = [];
        let flowCost = 0;

        let n = H.nodes.size;
        let hugeWeight = 1 + n*maxWeight;

        for (let [v, d] of nodeIter) {
            let vDemand = _.get(d, 'demand', 0);
            if (vDemand >= 0) {
                if (!H.hasEdge(r, v)) {
                    H.addEdge(r, v, {weight: hugeWeight, flow: vDemand});
                    artificialEdges.push([r, v]);
                    y[v] = hugeWeight;
                    T.addEdge(r, v);
                    flowCost += vDemand * hugeWeight;
                } else {
                    if (!"capacity" in H.getEdge(r, v) ||
                        vDemand <= H.getEdge(r, v).capacity) {
                        H.getEdge(r, v).flow = vDemand;
                        y[v] = _.get(H.getEdge(r, v), "weight", 0);
                        T.addEdge(r, v);
                        flowCost += vDemand * y[v];
                    } else {
                        let newLabel = -H.nodes.size;
                        H.addEdge(r, newLabel, {weight: hugeWeight, flow: vDemand});
                        H.addEdge(newLabel, v, {weight: hugeWeight, flow: vDemand});
                        artificialEdges.push([r, newLabel]);
                        artificialEdges.push([newLabel, v]);
                        y[v] = 2*hugeWeight;
                        y[newLabel] = hugeWeight;
                        T.addEdge(r, newLabel);
                        T.addEdge(newLabel, v);
                        flowCost += vDemand * y[v];
                    }
                }
            } else { // vDemand < 0
                if (!H.hasEdge(v, r)) {
                    H.addEdge(v, r, {weight: hugeWeight, flow: -vDemand});
                    artificialEdges.push([v, r]);
                    y[v] = -hugeWeight;
                    T.addEdge(v, r);
                    flowCost += vDemand * hugeWeight;
                } else {
                    if (!"capacity" in H.getEdge(v, r) ||
                        -vDemand <= H.getEdge(v, r).capacity) {
                        H.getEdge(v, r).flow = vDemand;
                        y[v] = -_.get(H.getEdge(v, r), "weight", 0);
                        T.add_edge(v, r);
                        flowCost += vDemand * y[v];
                    } else {
                        let newLabel = -H.nodes.size;
                        H.addEdge(v, newLabel, {weight: hugeWeight, flow: -vDemand});
                        H.addEdge(newLabel, r, {weight: hugeWeight, flow: -vDemand});
                        artificialEdges.push([v, newLabel]);
                        artificialEdges.push([newLabel, r]);
                        y[v] = -2*hugeWeight;
                        y[newLabel] = -hugeWeight;
                        T.addEdge(v, newLabel);
                        T.addEdge(newLabel, r);
                        flowCost += vDemand * y[v];
                    }
                }

            }
        }

        return [H, T, y, artificialEdges, flowCost, r];
    }

    _findEnteringEdge(c) {
        let newEdge = [];

        for (let [u, sinks] of this.edges) {
            for (let [v, d] of sinks) {
                //if (_.get(d, 'flow', 0) == 0) {
                if (d.flow === 0 || d.flow === undefined) {
                    if (c[u][v] < 0) {
                        newEdge = [u, v];
                        return newEdge;
                    }
                } else {
                    if ("capacity" in d &&
                        //_.get(d, 'flow', 0) == d.capacity &&
                        ((d.flow === undefined && d.capacity === 0) ||
                         d.flow === d.capacity) &&
                        c[u][v] > 0) {
                        newEdge = [u, v];
                        return newEdge;
                    }
                }
            }
        }
        return newEdge;
    }

    _findLeavingEdge(T, cycle, newEdge, reverse) {
        let eps = false;
        let leavingEdge = [];

        if (cycle.length == 3) {
            let [u, v] = newEdge;

            if (reverse) {
                if (_.get(this.getEdge(u, v), "flow", 0) >
                    _.get(this.getEdge(v, u), "flow", 0)) {
                    return [[v, u], _.get(this.getEdge(v, u), "flow", 0)];
                } else {
                    return [[u, v], _.get(this.getEdge(u, v), "flow", 0)];
                }
            } else {
                let uv_res = (_.get(this.getEdge(u, v), "capacity", 0) -
                              _.get(this.getEdge(u, v), "flow", 0));
                let vu_res = (_.get(this.getEdge(v, u), "capacity", 0) -
                              _.get(this.getEdge(v, u), "flow", 0));

                if (uv_res > vu_res) {
                    return [[v, u], vu_res];
                } else {
                    return [[u, v], uv_res];
                }
            }
        }

        for (let i = 0; i < cycle.length-1; i++) {
            let u = cycle[i];

            let edgeCapacity = false;
            let edge = [];
            let v = cycle[i+1];

            if ((u == newEdge[0] && v == newEdge[1]) ||
                T.hasEdge(u, v)) {
                if ("capacity" in this.getEdge(u, v)) {
                    edgeCapacity = (this.getEdge(u, v).capacity -
                                    _.get(this.getEdge(u, v), "flow", 0));
                    edge = [u, v];
                }
            } else {
                edgeCapacity = _.get(this.getEdge(v, u), "flow", 0);
                edge = [v, u];
            }

            if (edge.length > 0) {
                if (leavingEdge.length > 0) {
                    if (edgeCapacity < eps) {
                        eps = edgeCapacity;
                        leavingEdge = edge;
                    }
                } else {
                    eps = edgeCapacity;
                    leavingEdge = edge;
                }
            }
        }

        return [leavingEdge, eps];
    }

}
