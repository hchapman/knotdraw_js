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
/******/ 	return __webpack_require__(__webpack_require__.s = 1);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */,
/* 1 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
Object.defineProperty(__webpack_exports__, "__esModule", { value: true });
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__lib_shadow_js__ = __webpack_require__(2);
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_1__lib_orthemb_js__ = __webpack_require__(3);
importScripts('lalolib/lalolib.js');
importScripts('lodash.min.js');




function leastSquares(X /* : Matrix */, Y /* : Matrix */) /* : leastSquares */ {
    //console.log("X", X, inv(X), det(X), qr(X, true), Y);
    let betaHat = solve(mul(X, transpose(X)), mul(X, Y));

    return betaHat;
}

class ForceLinkDiagram {
    /* Link diagram embedding improved by ImPrEd */
    constructor (verts, edges, faces) {
        this.verts = verts;
        this.edges = edges;
        this.faces = faces;

        //console.log("+++++++");
        //console.log(verts);
        //console.log(edges);
        //console.log(faces);

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

        this.delta = 2;
        this.gamma = 5;

        this.dbar = 3*this.delta;

        this.aExp = 1;
        this.erExp = 2;

        this.calculateSurroundingEdges();
    }

    distance(u, v) {
        return norm(sub(u, v));
    }

    forceAvert(u, v) {
        //console.log(this.distance(u,v));
        return mul(Math.pow(this.distance(u, v)/this.delta, this.aExp),
                   sub(v, u));
    }

    forceRvert(u, v) {
        let d = this.distance(u, v);
        return mul(Math.pow(this.delta/d, this.erExp),
                   sub(u, v));
    }

    computeVe(v, a, b) {
        let m = (a[1] - b[1])/(a[0] - b[0]);
        let n = -1 / m;
        let c = a[1] - m*a[0];
        let d = v[1] - n*v[0];

        let x = (d - c) / (m - n);
        return [x, m*x + c];
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
        let edges = [];
        for (let face of this.faces) {
            if (face.includes(ui)) {
                if (ui == 63) {
                    console.log(face);
                }
                for (let i = 0; i < face.length-1; i++) {
                    console.assert(this.edges.filter(
                        e => ((e[0] == face[i] && e[1] == face[i+1]) ||
                              (e[1] == face[i] && e[0] == face[i+1]))).length > 0);
                    edges.push([face[i], face[i+1]]);
                }
                edges.push([face[face.length-1], face[0]]);
            }
        }
        return edges;
    }

    calculateSurroundingEdges() {
        this.surrEdges = [];
        for (let i = 0; i < this.verts.length; i++) {
            this.surrEdges[i] = this.surroundingEdges(i);
        }
    }

    move (ui, FUx, FUy, MU) {
        let i;
        if (FUx >= 0) {
            if (FUy >= 0) {
                if (FUx >= FUy) {
                    i = 0;
                } else {
                    i = 1;
                }
            } else {
                if (FUx >= -FUy) {
                    i = 7;
                } else {
                    i = 6;
                }
            }
        } else {
            if (FUy >= 0) {
                if (-FUx >= FUy) {
                    i = 3;
                } else {
                    i = 2;
                }
            } else {
                if (-FUx >= -FUy) {
                    i = 4;
                } else {
                    i = 5;
                }
            }
        }

        let FU = [FUx, FUy];

        let fU = norm(FU);
        let du;
        if (fU <= MU[i]) {
            du = FU;
        } else {
            du = mul(MU[i]/fU, FU);
        }

        //if (ui == 9) {
        //    console.log(ui, i, du, FU, MU[i]);
        //}

        //console.log(this.verts[ui], du, FU, MU);
        this.verts[ui][0] += du[0];
        this.verts[ui][1] += du[1];
        //console.log(FU, MU[i], du, this.verts[ui]);
    }

    update() {
        //console.log(this.adjMap);
        let FX = zeros(this.verts.length);
        let FY = zeros(this.verts.length);
        let M = [];
        for (let i = 0; i < this.verts.length; i++) {
            M.push([this.dbar, this.dbar, this.dbar, this.dbar,
                    this.dbar, this.dbar, this.dbar, this.dbar]);
        }

        let barycenter = mul(1/this.verts.length, sum(this.verts, 2));

        for (let ui = 0; ui < this.verts.length; ui++) {
            // Calculate gravity force
            let db = sub(barycenter, this.verts[ui]);
            let nDb = norm(db);
            FX[ui] += db[0]/nDb;
            FY[ui] += db[1]/nDb;

            // Calculate total node-node repulsive force
            for (let vi = 0; vi < this.verts.length; vi++) {
                if (ui != vi) {
                    if (this.distance(ui, vi) >= 3*this.delta) {
                        continue;
                    }

                    if (this.adjMap[ui].length == 2) {
                        if (this.adjMap[ui].includes(vi)) {
                            continue;
                        }
                    }

                    let F = this.forceRvert(this.verts[ui], this.verts[vi]);
                    //console.log("Fnnr", F);
                    //console.log(ui, vi, this.verts[ui], this.verts[vi], "Fnnr", F);
                    if (!isNaN(F[0])) {
                        FX[ui] += F[0];
                        FY[ui] += F[1];
                    }
                }
            }

            // calculate edge attractive force
            for (let vi of this.adjMap[ui]) {
                let F = this.forceAvert(this.verts[ui], this.verts[vi]);

                FX[ui] += F[0];
                FY[ui] += F[1];
            }

            // calculate node-edge repulsive force
            for (let edge of this.surrEdges[ui]) {
                let [ai, bi] = edge;
                if (ui == ai || ui == bi) {
                    continue;
                }
                let ve = this.computeVe(
                    this.verts[ui], this.verts[ai], this.verts[bi]);

                if (this.veOnEdge(ve, this.verts[ai], this.verts[bi])) {
                    let F = this.forceRedge(
                        this.verts[ui], this.verts[ai], this.verts[bi], ve);
                    if (!isNaN(F[0])) {
                        FX[ui] += F[0];
                        FY[ui] += F[1];
                    }
                }
            }

            let MU = M[ui];
            //console.log("Surr:", this.surrEdges);

            for (let edge of this.surrEdges[ui]) {
                let [ai, bi] = edge;
                if (ui == ai || ui == bi) {
                    continue;
                }
                let ve = this.computeVe(
                    this.verts[ui], this.verts[ai], this.verts[bi]);

                let cv;

                if (ui == 0 && ai == 5 && bi == 2) {
                    //console.log(this.verts[ai], this.verts[bi], ve);
                }
                //console.log("v-e", ui, ai, bi);
                if (this.veOnEdge(ve, this.verts[ai], this.verts[bi])) {
                    cv = sub(ve, this.verts[ui]);

                    let i;
                    if (cv[0] >= 0) {
                        if (cv[1] >= 0) {
                            if (cv[0] >= cv[1]) {
                                i = 0;
                            } else {
                                i = 1;
                            }
                        } else {
                            if (cv[0] >= -cv[1]) {
                                i = 7;
                            } else {
                                i = 6;
                            }
                        }
                    } else {
                        if (cv[1] >= 0) {
                            if (-cv[0] >= cv[1]) {
                                i = 3;
                            } else {
                                i = 2;
                            }
                        } else {
                            if (-cv[0] >= -cv[1]) {
                                i = 4;
                            } else {
                                i = 5;
                            }
                        }
                    }

                    let maxR = norm(cv)/2.1;
                    //console.log("???", cv);

                    //console.log(MU, maxR, Math.cos(Math.atan2(cv[1], cv[0])));
                    let ell = (i+4)%8;
                    for (let j = 0; j < MU.length; j++) {
                        if ((i-j+8)%8 == 0) {
                            MU[j] = min(MU[j], maxR);
                        } else if ((i-j+8)%8 == 1 || (i-j+8)%8 == 2) {
                            MU[j] = min(MU[j], maxR /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j+1)*Math.PI/4));
                        } else if ((i-j+8)%8 == 6 || (i-j+8)%8 == 7) {
                            MU[j] = min(MU[j], maxR /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    for (let j = 0; j < MU.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[ai][j] = min(M[ai][j], maxR);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[ai][j] = min(M[ai][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[ai][j] = min(M[ai][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    for (let j = 0; j < MU.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[bi][j] = min(M[bi][j], maxR);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[bi][j] = min(M[bi][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[bi][j] = min(M[bi][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j)*Math.PI/4));
                        }
                    }

                } else {
                    let va = sub(this.verts[ai], this.verts[ui]);
                    let vb = sub(this.verts[bi], this.verts[ui]);
                    if (norm(va) < norm(vb)) {
                        cv = va;
                    } else {
                        cv = vb;
                    }

                    let i;
                    if (cv[0] >= 0) {
                        if (cv[1] >= 0) {
                            if (cv[0] >= cv[1]) {
                                i = 0;
                            } else {
                                i = 1;
                            }
                        } else {
                            if (cv[0] >= -cv[1]) {
                                i = 7;
                            } else {
                                i = 6;
                            }
                        }
                    } else {
                        if (cv[1] >= 0) {
                            if (-cv[0] >= cv[1]) {
                                i = 3;
                            } else {
                                i = 2;
                            }
                        } else {
                            if (-cv[0] >= -cv[1]) {
                                i = 4;
                            } else {
                                i = 5;
                            }
                        }
                    }

                    let maxR = norm(cv)/2.1;
                    //console.log("???", cv);

                    //console.log(MU, maxR, Math.cos(Math.atan2(cv[1], cv[0])));
                    let ell = (i+4)%8;
                    for (let j = 0; j < MU.length; j++) {
                        if ((i-j+8)%8 == 0) {
                            MU[j] = min(MU[j], maxR);
                        } else if ((i-j+8)%8 == 1 || (i-j+8)%8 == 2) {
                            MU[j] = min(MU[j], maxR /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j+1)*Math.PI/4));
                        } else if ((i-j+8)%8 == 6 || (i-j+8)%8 == 7) {
                            MU[j] = min(MU[j], maxR /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    let m = cv[1]/cv[0]; // Slope of cv
                    let n = -1 / m; // Slope of l

                    for (let j = 0; j < MU.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[ai][j] = min(M[ai][j], maxR);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[ai][j] = min(M[ai][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[ai][j] = min(M[ai][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    for (let j = 0; j < MU.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[bi][j] = min(M[bi][j], maxR);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[bi][j] = min(M[bi][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[bi][j] = min(M[bi][j], maxR /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j)*Math.PI/4));
                        }
                    }
                }

            }
            //if (ui == 0) console.log("MU0", MU, FX[ui], FY[ui]);
        }

        //console.log("Fx", FX);
        for (let ui in this.verts) {

            this.move(ui, FX[ui], FY[ui], M[ui]);
        }
    }
}

function sleep(millis)
{
    var date = new Date();
    var curDate = null;
    do { curDate = new Date(); }
    while(curDate-date < millis);
}


var workerFunctions = {
    setLinkDiagram: function(sigma, crossBend) {
        self.shadow = new __WEBPACK_IMPORTED_MODULE_0__lib_shadow_js__["a" /* default */](sigma);
        self.orthshadow = new __WEBPACK_IMPORTED_MODULE_1__lib_orthemb_js__["a" /* default */](self.shadow);

        workerFunctions.embedDiagram();
    },

    embedDiagram: function() {
        let tstart = Date.now();

        let thresh = 5e-10;
    }
}

onmessage = function(e) {
    workerFunctions[e.data.function](...e.data.arguments);
}


/***/ }),
/* 2 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
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

    newArc(idx) {
        this.arcs[idx] = {index: idx, edge:undefined, edgepos:undefined, vert:undefined, vertpos:undefined};
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
            //console.log(i)
            //console.log(this.arcs)
            this.arcs[ais[i]].vert = parseInt(idx);
            this.arcs[ais[i]].vertpos = parseInt(i);
        }
    }
}
/* harmony export (immutable) */ __webpack_exports__["a"] = LinkShadow;



/***/ }),
/* 3 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
/* harmony import */ var __WEBPACK_IMPORTED_MODULE_0__digraph_js__ = __webpack_require__(4);


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

class OrthogonalDiagramEmbedding {
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
/* harmony export (immutable) */ __webpack_exports__["a"] = OrthogonalDiagramEmbedding;



/***/ }),
/* 4 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
class DiGraph {
    /* Lifted from networkx's network_simplex routines */

    constructor() {
        this.nodes = new Map();
        this.edges = new Map();
    }

    copy() {
        let G = new DiGraph();
        G.nodes = new Map(this.nodes);
        for (let [u, sinks] of this.edges) {
            G.edges.set(u, new Map(this.edges.get(u)));
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
    }

    removeNode(name) {
        this.nodes.delete(name);
        this.edges.delete(name);
        for (let [source, sinks] of this.edges) {
            sinks.delete(name);
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

        if (attrs === undefined) {
            return;
        }

        for (let [k, v] of Object.entries(attrs)) {
            this.edges.get(source).get(sink)[k] = v;
        }
    }

    removeEdge(source, sink) {
        if (this.edges.has(source)) {
            this.edges.get(source).delete(sink);
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

    getReverseEdges() {
        let rev_edges = new Map();
        for (let [v, d] of this.nodes) {
            rev_edges.set(v, new Map());
        }
        for (let [u, sinks] of this.edges) {
            for (let [v, d] of sinks) {
                if (!rev_edges.has(v)) {
                    rev_edges.set(v, new Map());
                }
                rev_edges.get(v).set(u, {});
            }
        }
        return rev_edges;
    }

    treePath(start, stop) {
        // Requires that this is a tree to avoid real programming...

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

    *depthFirstSearch(start) {
        let stack = [start];
        let seen = new Set(stack);

        let rev_edges = this.getReverseEdges();

        while (stack.length > 0) {
            let current = stack.pop();
            for (let v in this.edges.get(current).keys()) {
                if (!seen.has(v)) {
                    stack.push(v);
                    seen.add(v);
                }
            }
            for (let v in rev_edges.get(current).keys()) {
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
            for (let [v, d] of sinks) {
                c[[u, v]] = _.get(d, 'weight', 0) + y[u] - y[v];
            }
        }

        let reverse;
        while (true) {
            let newEdge = H._findEnteringEdge(c);
            if (newEdge.length != 2) {
                break;
            }

            let cycleCost = Math.abs(c[newEdge]);

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

            if (_.get(H.getEdge(...newEdge), "flow", 0) == 0) {
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
                    _.set(H.getEdge(u, v), "flow", _.get(H.getEdge(u, v), "flow", 0) + eps);
                    _.set(H.getEdge(v, u), "flow", _.get(H.getEdge(v, u), "flow", 0) + eps);
                } else {
                    for (let j = 0; j < cycle.length-1; j++) {
                        let u = cycle[j];
                        let v = cycle[j+1];
                        if ((u == newEdge[0] && v == newEdge[1]) ||
                            T.hasEdge(u, v)) {
                            _.set(H.getEdge(u, v), "flow", _.get(H.getEdge(u, v), "flow", 0) + eps);
                        } else {
                            _.set(H.getEdge(v, u), "flow", _.get(H.getEdge(v, u), "flow", 0) - eps);
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
                        y[v] += c[newEdge];
                    }
                } else {
                    for (let v of notR.nodes.keys()) {
                        y[v] -= c[newEdge];
                    }
                }

                for (let [u, sinks] of H.edges) {
                    for (let [v, d] of sinks) {
                        if (notR.nodes.has(u) || notR.nodes.has(v)) {
                            c[[u, v]] = _.get(H.getEdge(u, v), "weight", 0) + y[u] - y[v];
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
                if (_.get(d, 'flow', 0) == 0) {
                    if (c[[u, v]] < 0) {
                        newEdge = [u, v];
                        return newEdge;
                    }
                } else {
                    if ("capacity" in d &&
                        _.get(d, 'flow', 0) == d.capacity &&
                        c[[u, v]] > 0) {
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



/***/ })
/******/ ]);
//# sourceMappingURL=orth_worker.js.map