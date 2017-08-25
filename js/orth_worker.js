/******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__webpack_require__.d = function(exports, name, getter) {
/******/ 		if(!__webpack_require__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, {
/******/ 				configurable: false,
/******/ 				enumerable: true,
/******/ 				get: getter
/******/ 			});
/******/ 		}
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__webpack_require__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__webpack_require__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(__webpack_require__.s = 4);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
class Point {

    constructor(x, y){
        this.x = x;
        this.y = y;
    }

    getPosition() {
        return {
            x: this.x,
            y: this.y
        };
    }

    increaseX(amount) {
        this.x += amount;
    }

    decreaseX(amount) {
        this.x -= amount;
    }

    increaseY(amount) {
        this.y += amount;
    }

    decreaseY(amount) {
        this.y -= amount;
    }
}
/* harmony export (immutable) */ __webpack_exports__["a"] = Point;



/***/ }),
/* 1 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
class DiGraph {
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
/* harmony export (immutable) */ __webpack_exports__["a"] = DiGraph;



/***/ }),
/* 2 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__point__ = __webpack_require__(0);


class Quad {

    constructor(origin, width, height) {
        this.origin = origin;
        this.width = width;
        this.height = height;
        this.xmin = origin.x;
        this.xmax = origin.x + width;
        this.ymin = origin.y;
        this.ymax = origin.y + height;
    }

    static copy(quad) {
        let origin = new __WEBPACK_IMPORTED_MODULE_0__point__["a" /* default */](quad.origin.x, quad.origin.y);
        let width = quad.width;
        let height = quad.height;

        return new Quad(origin, width, height);
    }

    intersects(quad) {
        if (this.xmin < quad.xmax && this.xmax > quad.xmin &&
            this.ymin < quad.ymax && this.ymax > quad.ymin) return true;
        return false;
    }

    containsPoint(point) {
        return (point.x <= this.xmax && point.x >= this.xmin) &&
            (point.y <= this.ymax && point.y >= this.ymin);
    }
}
/* harmony export (immutable) */ __webpack_exports__["a"] = Quad;



/***/ }),
/* 3 */,
/* 4 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
Object.defineProperty(__webpack_exports__, "__esModule", { value: true });
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__lib_shadow_js__ = __webpack_require__(5);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__lib_orthemb_js__ = __webpack_require__(6);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__lib_forcediagram_js__ = __webpack_require__(8);
importScripts('lalolib/lalolib.js');
importScripts('lodash.min.js');

importScripts('planarmap-em/libplanarmap-em.js');





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
        self.shadow = new __WEBPACK_IMPORTED_MODULE_0__lib_shadow_js__["a" /* default */](sigma);
        self.orthShadow = new __WEBPACK_IMPORTED_MODULE_1__lib_orthemb_js__["a" /* default */](self.shadow);

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
        self.force_shadow = new __WEBPACK_IMPORTED_MODULE_2__lib_forcediagram_js__["a" /* default */](verts, edges, self.faces, self.components);
        
        //console.log(faces);
        

        postMessage({
            function: "setLinkDiagram",
            arguments: [self.force_shadow]
        });

        workerFunctions.embedDiagram();
    },

    embedDiagram: function(n_steps=0) {
        let tstart = Date.now();

        let thresh = 5e-10;

        let curDate;
        self.n_steps = n_steps;//50;
        let max_steps = self.n_steps;

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
    },

    deleteMonogon: function(fi) {
        self.force_shadow.deleteMonogon(fi);
        self.force_shadow.initParams();
        self.force_shadow.numIter = 30;
        workerFunctions.embedDiagram(10);
    }
}

onmessage = function(e) {
    workerFunctions[e.data.function](...e.data.arguments);
}


/***/ }),
/* 5 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
class Arc {
    constructor(index) {
        this.index = index;
    }

    toString() {
        return this.index.toString();
    }

    valueOf() {
        return this.index;
    }
}

class LinkShadow {
    constructor(verts) {
        this.nv = verts.length;
        this.ne = this.nv*2;
        this.na = this.nv*4;

        this.arcs = [];
        this.edges = [];
        this.verts = [];

        // Edge i is canonically of the form [2i, 2i+1]
        for (let ei = 0; ei < this.ne; ei++) {
            //console.log(ei);
            this.newArc(2*ei);
            this.newArc(2*ei+1);
            //console.log(this.arcs)

            this.setEdge(ei, [2*ei, 2*ei+1]);
        }

        // Set verts by looking through verts
        for (let vi in verts) {
            this.setVert(vi, verts[vi]);
        }

        this.generateFaces();
        this.generateComponents();
        this.connectArcs();
    }

    copy() {
        let link = new LinkShadow(this.verts.map(
            v => v.map(a => a.index)
        ));

        return link;
    }

    generateFaces() {
        let leftArcs = new Set(this.arcs);
        this.faces = [];

        while(leftArcs.size > 0) {
            let startArc = Array.from(leftArcs).pop();

            let face = [];
            let arc = startArc;
            let _failsafe = 0;
            do {
                leftArcs.delete(arc);
                face.push(arc);
                arc.face = this.faces.length;
                //console.log(arc);

                let oArc = this.edges[arc.edge][(arc.edgepos+1)%2];
                arc = this.verts[oArc.vert][(oArc.vertpos+1)%4];
                //console.log(arc);
                _failsafe += 1;
                if (_failsafe > 500) {
                    console.log("Failure");
                    return this.faces;
                }
            } while (arc != startArc)

            //face.reverse();
            this.faces.push(face);
        }
        return this.faces;
    }

    generateComponents(oneOrient=true) {
        let leftArcs = new Set(this.arcs);

        this.components = [];
        while (leftArcs.size > 0) {
            let startArc = Array.from(leftArcs).pop();

            let component = this.component(startArc);
            for (let arc of component) {
                leftArcs.delete(arc);
                if (oneOrient) {
                    // Delete the other arc edge-opposite this one
                    leftArcs.delete(this.edges[arc.edge][(arc.edgepos+1)%2]);
                }
            }

            this.components.push(component);
        }
        return this.components;
    }

    component(arc) {
        let startArc = arc;

        let component = [];
        do {
            component.push(arc);

            let oArc = this.edges[arc.edge][(arc.edgepos+1)%2];
            arc = this.verts[oArc.vert][(oArc.vertpos+2)%4];
        } while (arc != startArc)

        return component;
    }

    outVertI(arc) {
        return arc.vert;
    }

    inVertI(arc) {
        return this.edgeOpposite(arc).vert;
    }

    edgeOpposite(arc) {
        return this.edges[arc.edge][(arc.edgepos+1)%2];
    }
    vertNext(arc) {
        return this.verts[arc.vert][(arc.vertpos+1)%4];
    }
    vertPrev(arc) {
        return this.verts[arc.vert][(arc.vertpos+3)%4];
    }
    vertOppo(arc) {
        return this.verts[arc.vert][(arc.vertpos+2)%4];
    }

    newArc(idx) {
        if (idx === undefined) { idx = this.arcs.length; }

        this.arcs[idx] = new Arc(idx);
        // {index: idx, edge:undefined, edgepos:undefined, vert:undefined, vertpos:undefined};
        return this.arcs[idx];
    }

    connectArcs() {
        /* Hook up the arcs so they know, locally, their linkings */

        for (let arc of this.arcs) {
            arc.edgeOpposite = this.edgeOpposite(arc);
            arc.vertNext = this.vertNext(arc);
            arc.vertPrev = this.vertPrev(arc);
            arc.vertOppo = this.vertOppo(arc);
        }
    }

    setEdge(idx, ais) {
        this.edges[idx] = ais.map((ai) => {return this.arcs[ai];}, this);

        for (let i in ais) {
            this.arcs[ais[i]].edge = idx;
            this.arcs[ais[i]].edgepos = parseInt(i);
        }
    }

    setVert(idx, ais) {
        this.verts[idx] = ais.map((ai) => {return this.arcs[ai];}, this);

        for (let i in ais) {
            if (ais[i] === undefined) { continue; }

            this.arcs[ais[i]].vert = parseInt(idx);
            this.arcs[ais[i]].vertpos = parseInt(i);
        }
    }
}
/* harmony export (immutable) */ __webpack_exports__["a"] = LinkShadow;



/***/ }),
/* 6 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__digraph_js__ = __webpack_require__(1);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__orthrep_js__ = __webpack_require__(7);



class Face {
    constructor (link, arcs, exterior=false) {
        this.link = link;
        this.arcs = arcs;
        this.exterior = exterior;

        // Real degree, i.e. the number of 4-valent vertices around this face
        this.degree = this.arcs.length;

        // Real crossings around this face, as opposed to joints
        this.cross = this.arcs.map(a => a.vert);

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

class OrthogonalDiagramEmbedding {
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
        let G = new __WEBPACK_IMPORTED_MODULE_0__digraph_js__["a" /* default */]();

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

        return G;
    }

    bend() {
        let flow = this.faceNetwork.minCostFlow()[1];

        for (let [a, flows] of flow) {
            for (let [b, w_a] of flows) {
                //console.log(a, b, w_a);
                if (!w_a || a == 's' || b == 's' || a == 't' || b == 't') {
                    continue;
                }

                let w_b = flow.get(b).get(a);

                let [A, B] = [this.faces[a], this.faces[b]];
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
            strands[si].vertOppo = strands[si+1];

            strands[si+1].vertNext = strands[si];
            strands[si+1].vertPrev = strands[si];
            strands[si+1].vertOppo = strands[si];

            this.shadow.setVert(this.shadow.verts.length,
                                [strands[si].index, undefined, strands[si+1].index, undefined]);
        }
    }

    repairComponents() {
        this.arcComponentMap = new Map();
        this.orderedArcs = [];
        this.components = [];

        for (let i = 0; i < this.shadow.components.length; i++) {
            for (let arc of this.shadow.components[i]) {
                this.arcComponentMap.set(arc.index, i);
                this.orderedArcs.push(arc);
            }

            let comp_root_a = this.shadow.components[i][0];
            let arc = comp_root_a;
            let component = [];
            do {
                component.push(arc);
                arc = arc.vertOppo.edgeOpposite;
            } while (arc.index != comp_root_a.index)
            this.components.push(component);
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
        return new __WEBPACK_IMPORTED_MODULE_1__orthrep_js__["a" /* default */](...this.orthogonalSpec());
    }
}
/* harmony export (immutable) */ __webpack_exports__["a"] = OrthogonalDiagramEmbedding;



/***/ }),
/* 7 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__digraph_js__ = __webpack_require__(1);


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

class OrthogonalRep {
    constructor(horiz_pairs, vert_pairs) {
        this.graph = new __WEBPACK_IMPORTED_MODULE_0__digraph_js__["a" /* default */]();
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
        let H = new __WEBPACK_IMPORTED_MODULE_0__digraph_js__["a" /* default */]();
        this.graph.nodes.forEach((data, v) =>
                                 H.addNode(v));
        for (let [source, sink, data] of this.graph.edgeGen()) {
            if (data.kind == kind) { H.addEdge(source, sink); }
        }

        let maximalChains = H.weakComponents();
        let vertexToChain = elementMap(maximalChains);

        let D = new __WEBPACK_IMPORTED_MODULE_0__digraph_js__["a" /* default */]();
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
/* harmony export (immutable) */ __webpack_exports__["a"] = OrthogonalRep;



/***/ }),
/* 8 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__Quadtree_es6_quadtree_js__ = __webpack_require__(9);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__Quadtree_es6_math_point_js__ = __webpack_require__(0);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_2__Quadtree_es6_math_quad_js__ = __webpack_require__(2);




class ForceLinkVertex {
    constructor (idx, coord) {
        this.index = idx;
        this.coord = coord;
        this.edges = new Map();
        this.adj = new Map();
    }
}

class ForceLinkDiagram {
    /* Link diagram embedding improved by ImPrEd */
    constructor (verts, edges, faces, components) {
        this.verts = new Map(verts.map(
            (v, vi) => [vi, new ForceLinkVertex(vi, v)]));

        this.edges = edges;
        for (let edge of this.edges) {
            this.setSymmEdge(edge[0], edge[1], edge);
        }
        this.faces = faces;
        this.components = components;

        this.freeVi = [];
        this.maxVi = Math.max(...verts.keys());

        this.adjMap = {};
        for (let edge of edges) {
            let [a, b] = edge;
            if (a in this.adjMap) {
                this.adjMap[a].push(b);
            } else {
                this.adjMap[a] = [b];
            }

            if (b in this.adjMap) {
                this.adjMap[b].push(a);
            } else {
                this.adjMap[b] = [a];
            }
        }

        let N = this.verts.size;
        function index(edge) {
            return edge[0] * N + edge[1];
        }
        function deindex(index) {
            return [(index - index%N) / N, index % N];
        }
        this.paths = [];
        let edgeSearch = new Set(this.edges.map(e => index(e)));
        while (edgeSearch.size) {
            let path;
            let edge = deindex(edgeSearch.values().next().value);
            edgeSearch.delete(index(edge));

            let [prev, next] = edge;
            path = [prev, next];

            while (this.adjMap[prev].length == 2) {
                let tmp = prev;
                prev = (this.adjMap[prev][0] == next ?
                        this.adjMap[prev][1] : this.adjMap[prev][0]);
                next = tmp;
                path.unshift(prev);

                edgeSearch.delete(index([prev, next]));
                edgeSearch.delete(index([next, prev]));
            }

            [prev, next] = edge;
            while (this.adjMap[next].length == 2) {
                let tmp = next;
                next = (this.adjMap[next][0] == prev ?
                        this.adjMap[next][1] : this.adjMap[next][0]);
                prev = tmp;
                path.push(next);

                edgeSearch.delete(index([prev, next]));
                edgeSearch.delete(index([next, prev]));
            }

            this.paths.push(path);
        }

        this.initParams();
        this.calculateSurroundingEdges();

        this._cacheDistances;
    }

    initParams() {
        this.delta = 2;
        this.gamma = 5;

        this.dbar = 3*this.delta;

        this.alpha = 2*this.delta;
        this.beta = 3*this.delta;

        //this.flexEvery = 2;

        this.aExp = 1;
        this.erExp = 2;

        this.numIter = 0;
    }

    getEdge(ai, bi) {
        return this.verts.get(ai).edges.get(bi);
    }

    setEdge(ai, bi, edge) {
        this.verts.get(ai).edges.set(bi, edge);
    }

    setSymmEdge(ai, bi, edge) {
        this.setEdge(ai, bi, edge); this.setEdge(bi, ai, edge);
    }

    deleteVertex(vi) {
        /* Deletes a vert, breaking paths, components, and faces (you must fix) */
        this.verts.delete(vi);

        this.surrEdges.delete(vi);
        this.freeVi.push(vi);

        // Delete this vertex from any faces it is in
    }

    deleteEdge(edge) {
        /* Deletes an edge, breaking paths, components, and faces (you must fix) */

        let [ai, bi] = edge;

        // Delete this edge from the matrix
        this.verts.get(ai).edges.delete(bi);
        this.verts.get(bi).edges.delete(ai);

        // Remove from edges
        this.edges.splice(this.edges.indexOf(edge), 1);

        // ai is no longer adjacent to bi
        this.adjMap[ai].splice(this.adjMap[ai].indexOf(bi), 1);
        this.adjMap[bi].splice(this.adjMap[bi].indexOf(ai), 1);

        // Update surrEdges; faces is potentially broken after this op
        for (let f of this.faces) {
            // check if this edge is in this face
            let fi = f.findIndex((v, i, f) => f[i]==ai && f[(i+1)%f.length]==bi);
            if (fi == -1) {
                fi = f.findIndex((v, i, f) => f[i]==bi && f[(i+1)%f.length]==ai);
            }
            if (fi == -1) { continue; }

            // Update surrEdges for vertices in f which AREN'T vi
            for (let ui of f) {
                if (this.surrEdges.has(ui))
                    this.surrEdges.get(ui).delete(edge);
            }
        }
    }

    deleteMonogon(fi) {
        let delFace = this.faces[fi].slice(); // copy of face to delete
        let monoXi; // The unique crossing vertex index in this monogon
        let monoPi; // The unique path which consists of this loop, to delete

        monoXi = delFace.find(vi => this.adjMap[vi].length == 4);

        for (let i = 0; i < delFace.length; i++) {
            let ai, bi;
            if (i == delFace.length-1) {
                [ai, bi] = [delFace[i], delFace[0]];
            } else {
                [ai, bi] = delFace.slice(i, i+2);
            }
            let delEdge = this.getEdge(ai, bi);

            // Delete edge ai->bi
            this.deleteEdge(delEdge);
        }

        for (let ai of delFace) {
            if (!this.adjMap[ai].length) {
                // Vertex ai is stranded, delete it
                this.deleteVertex(ai);

                // Find the path for this monogon, to delete later
                if (monoPi === undefined) {
                    monoPi = this.paths.findIndex(p => p.includes(ai));
                }

                // Remove ai from its component
                this.components.forEach(c => {
                    let cai = c.indexOf(ai);
                    if (cai != -1) {
                        c.splice(cai, 1);
                    }
                });
            }
        }

        // Delete this path
        this.paths.splice(monoPi, 1);

        // monoXi is in two paths; we want to combine them.
        let path1_i = this.paths.findIndex(p => p.includes(monoXi));
        let path2_i = path1_i+1 + this.paths.slice(path1_i+1).findIndex(p => p.includes(monoXi));
        let [p1, p2] = [this.paths[path1_i], this.paths[path2_i]];
        let x_p1i = p1.indexOf(monoXi); // either 0 or length-1
        if (x_p1i == 0) { p1.reverse(); } // crossing is now at end of p1
        if (p1.indexOf(monoXi) == 0) {
            // edge case--self loop edge (impossible) or
            // (actual) -> single twist diagram is being untwisted
            // sufficient to just remove path2
            this.paths.splice(path2_i, 1);

        } else {

            let x_p2i = p2.indexOf(monoXi); // either 0 or length-1
            if (x_p2i != 0) { p2.reverse(); } // crossing is not at front of p2

            // remove paths from this.paths
            this.paths.splice(path2_i, 1);
            this.paths.splice(path1_i, 1);

            this.paths.push(p1.concat(p2.slice(1)));

        }

        for (let vi of delFace) {
            this.faces.forEach(f => {
                let fai = f.indexOf(vi);
                if (fai != -1) {
                    if (vi == monoXi) {
                        fai = f.indexOf(vi, fai+1);
                        if (fai == -1) {
                            return;
                        }
                    }
                    f.splice(fai, 1);
                }
            });
        }

        // Remove the crossing---once---from its component
        this.components.forEach(c => {
            let cai = c.indexOf(monoXi);
            if (cai != -1) {
                c.splice(cai, 1);
            }
        });

        // Delete this face
        this.faces.splice(fi, 1);
    }

    newVertIdx() {
        if (this.freeVi.length) {
            return this.freeVi.pop();
        }
        return ++this.maxVi;
    }

    distance(u, v) {
        if (this._cacheDistances !== undefined) {
            let cu = this._cacheDistances.get(u);
            if (cu !== undefined) {
                let cuv = cu.get(v);
                if (cuv !== undefined) {
                    return cuv;
                }
            } else {
                cu = new Map()
                this._cacheDistances.set(u, cu);
            }
            let d = norm(sub(u, v));
            cu.set(v, d);
            return d;
        }
        // this._cacheDistances === undefined --> no caching
        return norm(sub(u, v));
    }

    forceAvert(u, v) {
        return mul(Math.pow(this.distance(u, v)/this.delta, this.aExp),
                   sub(v, u));
    }

    forceRvert(u, v) {
        let d = this.distance(u, v);
        //let a = Math.pow(this.delta/d, this.erExp);
        //return [(u[0]-v[0])*a, (u[1]-v[1])*a];
        return mul(Math.pow(this.delta/d, this.erExp),
                   sub(u, v));
    }

    computeVe(v, a, b) {
        let m = (a[1] - b[1])/(a[0] - b[0]);
        let n = -1 / m;
        let c = a[1] - m*a[0];
        let d = v[1] - n*v[0];

        if (m == 0) {
            return [v[0], a[1]];
        } else if (n == 0) {

            return [a[0], v[1]];
        } else {
            let x = (d - c) / (m - n);
            return [x, m*x + c];
        }
    }

    veOnEdge(ve, a, b) {
        return (((ve[0] <= a[0] && ve[0] >= b[0]) ||
                 (ve[0] <= b[0] && ve[0] >= a[0])) &&
                ((ve[1] <= a[1] && ve[1] >= b[1]) ||
                 (ve[1] <= b[1] && ve[1] >= a[1])));
    }

    forceRedge(u, a, b, ve) {
        let d = this.distance(u, ve);
        if (d >= this.gamma) {
            // node and "virtual edge" too far
            return [0, 0];
        }

        return mul(-Math.pow(this.gamma - d, this.erExp)/d,
                   sub(ve, u));
    }

    surroundingEdges(ui) {
        // calculate the surrounding edges SUi
        let edges = new Set();
        for (let face of this.faces) {
            if (face.includes(ui)) {
                for (let i = 0; i < face.length-1; i++) {
                    console.assert(this.edges.some(
                        e => ((e[0] == face[i] && e[1] == face[i+1]) ||
                              (e[1] == face[i] && e[0] == face[i+1]))), this.edges, face, i, ui);
                    if (face[i] != ui && face[i+1] != ui) {
                        edges.add(this.getEdge(face[i], face[i+1]));
                    }
                    console.assert(!edges.has(undefined));
                }
                if (face[face.length-1] != ui && face[0] != ui) {
                    edges.add(this.getEdge(face[face.length-1], face[0]));
                }
                console.assert(!edges.has(undefined));
            }
        }
        console.assert(!edges.has(undefined));
        return edges;
    }

    calculateSurroundingEdges() {
        this.surrEdges = new Map();
        for (let [i, v] of this.verts) {
            this.surrEdges.set(i, this.surroundingEdges(i));
        }
    }

    octant(x, y) {
        if (x >= 0) {
            if (y >= 0) {
                if (x >= y) {
                    return 0;
                } else {
                    return 1;
                }
            } else {
                if (x >= -y) {
                    return 7;
                } else {
                    return 6;
                }
            }
        } else {
            if (y >= 0) {
                if (-x >= y) {
                    return 3;
                } else {
                    return 2;
                }
            } else {
                if (-x >= -y) {
                    return 4;
                } else {
                    return 5;
                }
            }
        }
    }

    move (ui, FU, MU) {
        console.assert (this.verts.get(ui) !== undefined, this.verts);
        let i = this.octant(FU[0], FU[1]);

        MU[i] = Math.max(0, MU[i]);
        let fU = norm(FU);
        let du;
        if (fU <= MU[i]) {
            du = FU;
        } else {
            du = mul(MU[i]/fU, FU);
        }

        if (norm(du) < 10*Number.EPSILON) { return; }

        this.verts.get(ui).coord[0] += du[0];
        this.verts.get(ui).coord[1] += du[1];
    }

    triangleIncludes(a, b, c, p,vi) {
        let v0, v1, v2;
        let dot00, dot01, dot11, dot20, dot21;

        v0 = sub(b, a);
        v1 = sub(c, a);
        v2 = sub(p, a);

        dot00 = mul(v0, v0);
        dot01 = mul(v0, v1);
        dot11 = mul(v1, v1);
        dot20 = mul(v2, v0);
        dot21 = mul(v2, v1);

        let invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        let u = (dot11 * dot20 - dot01 * dot21) * invDenom;
        let v = (dot00 * dot21 - dot01 * dot20) * invDenom;

        return (u >= 0) && (v >= 0) && (u+v < 1);
    }

    contract() {
        for (let path of this.paths) {
            let vArray = Array.from(this.verts.entries());
            let reContract = true;
            while (reContract) {
                reContract = false;
                for (let i = 0; i < path.length-2; i++) {
                    let [ai, bi, ci] = path.slice(i, i+3);
                    let [a, b, c] = [this.verts.get(ai).coord,
                                     this.verts.get(bi).coord,
                                     this.verts.get(ci).coord];
                    if (norm(sub(a, b)) < this.alpha) {
                        if (this.faces.filter(f => f.includes(ai)).some(f => f.length <= 3)) { continue; }

                        if (!vArray.some(
                            ([vi, v]) => (vi != ai && vi != bi && vi != ci && this.triangleIncludes(a,b,c,v,vi)))) {
                            // contract

                            // remove bi from its adjacent vertices at ai and ci
                            this.adjMap[ai].splice(this.adjMap[ai].indexOf(bi), 1);
                            this.adjMap[ci].splice(this.adjMap[ci].indexOf(bi), 1);
                            this.adjMap[ai].push(ci);
                            this.adjMap[ci].push(ai);
                            let nEdge = [ai, ci];
                            this.setSymmEdge(ai, ci, nEdge);
                            this.edges.push(nEdge);

                            // remove bi from any faces
                            let bFaces = this.faces.filter(f => f.includes(bi));
                            // first, update these faces' surrEdges
                            for (let f of bFaces) {
                                for (let ui of f) {
                                    for (let [ev, bEdge] of this.verts.get(bi).edges) {
                                        if (bEdge === undefined) { continue; }
                                        this.surrEdges.get(ui).delete(bEdge);
                                    }
                                    if (!nEdge.includes(ui)) {
                                        this.surrEdges.get(ui).add(nEdge);
                                    }
                                }
                            }

                            // finally, remove bi from the faces
                            bFaces.forEach(f => f.splice(f.indexOf(bi), 1));

                            // remove bi from any components.
                            this.components.filter(cp => cp.includes(bi)).forEach(
                                cp => {while(cp.includes(bi)) { cp.splice(cp.indexOf(bi), 1); }}
                            );

                            // remove edges including bi
                            this.edges.splice(this.edges.indexOf(this.getEdge(ai, bi), 1));
                            this.edges.splice(this.edges.indexOf(this.getEdge(ci, bi), 1));

                            this.verts.get(ai).edges.delete(bi);
                            this.verts.get(ci).edges.delete(bi);

                            // remove bi from this path
                            path.splice(path.indexOf(bi), 1);
                            reContract = true;

                            // delete the vertex at bi
                            this.verts.delete(bi);
                            delete this.adjMap[bi];
                            break; // Only contract one bend of a path at a time???
                        }
                    }
                }
                //break;
            }
        }
    }

    expand() {
        for (let path of this.paths) {
            for (let i = 0; i < path.length-1; i++) {
                let [ai, bi] = path.slice(i,i+2);
                let [a, b] = [this.verts.get(ai).coord,
                              this.verts.get(bi).coord];
                let elen = norm(sub(a, b));
                if (elen > this.beta) {
                    // Insert a new vertex along this path
                    let vi = this.newVertIdx();
                    let v = new ForceLinkVertex(vi, mul(.5, add(a, b)));
                    this.verts.set(vi, v);

                    // Edge to delete; which was split
                    let splitEdge = this.getEdge(ai, bi);

                    // Add the new edges to the diagram
                    let nEdgeA = [ai, vi];
                    let nEdgeB = [vi, bi];

                    // Update the edgematrix of edges
                    this.setSymmEdge(ai, vi, nEdgeA);
                    this.setSymmEdge(bi, vi, nEdgeB);

                    // Remove splitEdge from edges, then add in these new edges
                    this.edges.splice(this.edges.indexOf(splitEdge), 1);
                    this.edges.push(nEdgeA); this.edges.push(nEdgeB);

                    // vi is adjacent to the two endpoints ai, bi
                    this.adjMap[vi] = [ai, bi];

                    // ai is no longer adjacent to bi
                    this.adjMap[ai].splice(this.adjMap[ai].indexOf(bi), 1, vi);
                    this.adjMap[bi].splice(this.adjMap[bi].indexOf(ai), 1, vi);

                    // Add vi to faces adjacent to splitEdge
                    for (let f of this.faces) {
                        let fi = f.findIndex((v, i, f) => f[i]==ai && f[(i+1)%f.length]==bi);
                        if (fi == -1) {
                            fi = f.findIndex((v, i, f) => f[i]==bi && f[(i+1)%f.length]==ai);
                        }
                        if (fi == -1) { continue; }

                        // Update surrEdges for vertices in f which AREN'T vi
                        for (let ui of f) {
                            this.surrEdges.get(ui).delete(splitEdge);
                            if (ai != ui) {
                                this.surrEdges.get(ui).add(nEdgeA);
                            }
                            if (bi != ui) {
                                this.surrEdges.get(ui).add(nEdgeB);
                            }
                        }

                        f.splice(fi+1,0,vi);
                    }

                    // Add vi to any components
                    for (let comp of this.components) {
                        let fi = comp.findIndex((v, i, f) => f[i]==ai && f[(i+1)%f.length]==bi);
                        if (fi == -1) {
                            fi = comp.findIndex((v, i, f) => f[i]==bi && f[(i+1)%f.length]==ai);
                        }
                        if (fi == -1) { continue; }

                        if (fi != -1) {
                            comp.splice(fi+1,0,vi);
                        }
                    }

                    // Set surrEdges for vi
                    this.surrEdges.set(vi, this.surroundingEdges(vi));

                    // Add vi to this path
                    path.splice(i+1,0,vi);
                    break; // only expand once per path, for now?
                }
            }
        }
    }

    update() {

        // We can't trust distances while we contract/expand, maybe??
        this._cacheDistances = undefined;

        if (true) {
            this.contract();
        }

        // console.log(this.numIter);
        if (this.numIter < 20) {
            this.expand();
        }

        this.numIter++;

        // Distances are fixed until next move
        this._cacheDistances = new Map();

        let F = new Map(Array.from(this.verts.keys()).map(k => [k, [0,0]])); //zeros(this.verts.size);
        let M = new Map();
        for (let [vi, v] of this.verts) {
            M.set(vi, [this.dbar, this.dbar, this.dbar, this.dbar,
                       this.dbar, this.dbar, this.dbar, this.dbar]);
        }

        let barycenter = mul(1/this.verts.size, sum(Array.from(this.verts.values()).map(v => v.coord), 2));

        //for (let ui = 0; ui < this.verts.size; ui++) {
        for (let [ui, u_v] of this.verts) {
            let u = u_v.coord;
            let Fu = F.get(ui);

            // Calculate gravity force
            let db = sub(barycenter, u);
            let nDb = norm(db);
            Fu[0] += db[0]/nDb;
            Fu[1] += db[1]/nDb;

            // Calculate total node-node repulsive force
            for (let [vi, v_v] of this.verts) {
                if (ui != vi) {
                    if (this.paths.some(p => p.includes(ui) && p.includes(vi))) {
                        //console.log("consec");
                        continue;
                    }
                    let v = v_v.coord;
                    //if (this.adjMap[ui].length == 2) {
                    //    if (this.adjMap[ui].includes(vi)) {
                    //        continue;
                    //    }
                    //}

                    let dF = this.forceRvert(u, v);
                    if (!isNaN(dF[0])) {
                        Fu[0] += dF[0];
                        Fu[1] += dF[1];
                    }
                }
            }

            // calculate edge attractive force
            for (let vi of this.adjMap[ui]) {
                let v = this.verts.get(vi).coord;
                let dF = this.forceAvert(u, v);

                Fu[0] += dF[0];
                Fu[1] += dF[1];
            }

            // calculate node-edge repulsive force
            for (let edge of this.surrEdges.get(ui)) {
                let [ai, bi] = edge;
                if (ui == ai || ui == bi) {
                    continue;
                }
                let [a, b] = [this.verts.get(ai).coord, this.verts.get(bi).coord];
                let ve = this.computeVe(
                    u, a, b);

                if (this.veOnEdge(ve, a, b)) {
                    let dF = this.forceRedge(
                        u, a, b, ve);
                    if (!isNaN(dF[0])) {
                        Fu[0] += dF[0];
                        Fu[1] += dF[1];
                    }
                }
            }

            let MU = M.get(ui);

            for (let edge of this.surrEdges.get(ui)) {
                let [ai, bi] = edge;
                if (ui == ai || ui == bi) {
                    continue;
                }
                let [a_v, b_v] = [this.verts.get(ai), this.verts.get(bi)];
                let [a, b] = [a_v.coord, b_v.coord];
                let ve = this.computeVe(u, a, b);

                let cv;

                if (this.veOnEdge(ve, a, b)) {
                    cv = sub(ve, u);

                } else {
                    let va = sub(a, u);
                    let vb = sub(b, u);
                    if (norm(va) < norm(vb)) {
                        cv = va;
                    } else {
                        cv = vb;
                    }
                }

                let i = this.octant(cv[0], cv[1]);

                let maxR = norm(cv)/2.1;
                let cv_angle = Math.atan2(cv[1], cv[0]);

                let ell = (i+4)%8;
                for (let j = 0; j < MU.length; j++) {
                    if ((i-j+8)%8 == 0) {
                        MU[j] = Math.min(MU[j], maxR);
                    } else if ((i-j+8)%8 == 1 || (i-j+8)%8 == 2) {
                        MU[j] = Math.min(MU[j], maxR /
                                         Math.cos(cv_angle - (j+1)*Math.PI/4));
                    } else if ((i-j+8)%8 == 6 || (i-j+8)%8 == 7) {
                        MU[j] = Math.min(MU[j], maxR /
                                         Math.cos(cv_angle - (j)*Math.PI/4));
                    }
                }

                let cw_angle = cv_angle + Math.PI;

                let Ma = M.get(ai);
                for (let j = 0; j < Ma.length; j++) {
                    if ((ell-j+8)%8 == 0) {
                        Ma[j] = Math.min(Ma[j], maxR);
                    } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                        Ma[j] = Math.min(Ma[j], maxR /
                                         Math.cos(cw_angle - (j+1)*Math.PI/4));
                    } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                        Ma[j] = Math.min(Ma[j], maxR /
                                         Math.cos(cw_angle - (j)*Math.PI/4));
                    }
                }

                let Mb = M.get(bi);
                for (let j = 0; j < Mb.length; j++) {
                    if ((ell-j+8)%8 == 0) {
                        Mb[j] = Math.min(Mb[j], maxR);
                    } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                        Mb[j] = Math.min(Mb[j], maxR /
                                         Math.cos(cw_angle - (j+1)*Math.PI/4));
                    } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                        Mb[j] = Math.min(Mb[j], maxR /
                                         Math.cos(cw_angle - (j)*Math.PI/4));
                    }
                }

            }
        }

        // Move all verts based on their maximal movements and forces
        for (let [ui, u] of this.verts) {
            this.move(ui, F.get(ui), M.get(ui));
        }
    }
}
/* harmony export (immutable) */ __webpack_exports__["a"] = ForceLinkDiagram;



/***/ }),
/* 9 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__math_quad__ = __webpack_require__(2);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__math_point__ = __webpack_require__(0);



class QuadTree {
    constructor(bounds, maxObjects, maxLevels, level = 0) {
        this.bounds = bounds;
        this.maxObjects = maxObjects;

        this.maxLevels = maxLevels;
        this.level = level;

        this.objects = [];
        this.children = [];
    }

    set northWest(tree) {
        this.children[0] = tree;
    }

    set northEast(tree) {
        this.children[1] = tree;
    }

    set southWest(tree) {
        this.children[2] = tree;
    }

    set southEast(tree) {
        this.children[3] = tree;
    }

    get northWest() {
        return this.children[0];
    }

    get northEast() {
        return this.children[1];
    }

    get southWest() {
        return this.children[2];
    }

    get southEast() {
        return this.children[3];
    }


    insert(obj) {
        // if we have a sub-tree
        if (this.northWest !== undefined) {
            return this._getChild(obj.point).insert(obj);
        }

        // check to see the point is within bounds of the tree
        if (!this.bounds.containsPoint(obj.point)) {
            return false;
        }

        this.objects.push(obj);

        // split the tree if there are too many objects at this level
        if (this.objects.length > this.maxObjects && this.level < this.maxLevels) {
            if (this.northWest === undefined) {
                this._subdivide();
            }

            // add all the objects in this node to its children
            this.objects.forEach(obj => {
                this._getChild(obj.point).insert(obj);
            });

            this.objects = [];
        }
        return true;
    }

    queryRange(range, selector = null, cb = null) {
        let pointsInRange = [];

        if (!this.bounds.intersects(range)) {
            return pointsInRange;
        }

        this.objects.forEach(obj => {
            if (range.containsPoint(obj.point)) {
                pointsInRange.push(obj[selector] || obj);
            }
        });

        if (this.northWest !== undefined) {
            this.children.forEach(child => {
                pointsInRange = pointsInRange.concat(child.queryRange(range, selector));
            });
        }

        if (cb) {
            return cb(pointsInRange);
        }

        return pointsInRange;
    }

    clear() {
        this.objects = [];

        this.children.forEach(child => {
            if (child !== undefined) {
                child.clear();
            }
        });

        this.children = [];
    }

    _subdivide() {
        let nextLevel = this.level + 1;
        let subWidth = Math.round(this.bounds.width / 2);
        let subHeight = Math.round(this.bounds.height / 2);
        let x = Math.round(this.bounds.xmin);
        let y = Math.round(this.bounds.ymin);

        let nwBounds = new __WEBPACK_IMPORTED_MODULE_0__math_quad__["a" /* default */](new __WEBPACK_IMPORTED_MODULE_1__math_point__["a" /* default */](x, y), subWidth, subHeight);
        let neBounds = new __WEBPACK_IMPORTED_MODULE_0__math_quad__["a" /* default */](new __WEBPACK_IMPORTED_MODULE_1__math_point__["a" /* default */](x + subWidth, y), subWidth, subHeight);
        let swBounds = new __WEBPACK_IMPORTED_MODULE_0__math_quad__["a" /* default */](new __WEBPACK_IMPORTED_MODULE_1__math_point__["a" /* default */](x, y + subHeight), subWidth, subHeight);
        let seBounds = new __WEBPACK_IMPORTED_MODULE_0__math_quad__["a" /* default */](new __WEBPACK_IMPORTED_MODULE_1__math_point__["a" /* default */](x + subWidth, y + subHeight), subWidth, subHeight);

        this.northWest = new QuadTree(nwBounds, this.maxObjects, this.maxLevels, nextLevel);
        this.northEast = new QuadTree(neBounds, this.maxObjects, this.maxLevels, nextLevel);
        this.southWest = new QuadTree(swBounds, this.maxObjects, this.maxLevels, nextLevel);
        this.southEast = new QuadTree(seBounds, this.maxObjects, this.maxLevels, nextLevel);
    }

    _getChild(point) {
        if (this.northWest.bounds.containsPoint(point)) return this.northWest;
        if (this.northEast.bounds.containsPoint(point)) return this.northEast;
        if (this.southWest.bounds.containsPoint(point)) return this.southWest;
        if (this.southEast.bounds.containsPoint(point)) return this.southEast;
    }
}
/* unused harmony export default */



/***/ })
/******/ ]);
//# sourceMappingURL=orth_worker.js.map