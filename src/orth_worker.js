importScripts('lalolib/lalolib.js');
importScripts('lodash.min.js');

importScripts('planarmap-em/libplanarmap-em.js');

import LinkShadow from "./lib/shadow.js";
import OrthogonalDiagramEmbedding from "./lib/orthemb.js";

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
        let edges = [];
        for (let face of this.faces) {
            if (face.includes(ui)) {
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

        MU[i] = Math.max(0, MU[i]);
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

self.randomDiagram = Module.cwrap('randomDiagram', 'number',
                                  ['number', 'number', 'number', 'number', 'number', 'number']);

function randomDiagram(n_verts, n_comps, max_att, type) {
    let vertPtr = Module._malloc(4);

    let nVerts = self.randomDiagram(n_verts, n_comps, max_att, type,
                                    Math.random()*2**32, vertPtr);
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
            edges.push([source, sink]);
        }
        let faces = rep.faces.map(f => f.evPairs.map(ev => ev[1]));
        //faces.forEach(f => f.reverse());

        self.force_shadow = new ForceLinkDiagram(verts, edges, faces);
        self.faces = self.orthShadow.faces.map(f => {
            let face = f.arcs.map(a => a.vert);
            face.exterior = f.exterior;
            return face;
        });
        console.log(faces);
        self.components = self.orthShadow.components.map(c => c.map(a => a.vert));

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
        let n_steps = 50;

        for (let i = 0; i < n_steps; i++) {
            let procStart = Date.now();

            postMessage({
                function: "setLinkDiagram",
                arguments: [self.force_shadow]
            });

            self.force_shadow.update();
            self.force_shadow.a_exp -= (1 - 0.4)/n_steps;
            self.force_shadow.re_exp += (4 - 2)/n_steps;
            self.force_shadow.dbar -= (3*self.force_shadow.delta)/n_steps;

            //do { curDate = Date.now(); }
            //while( curDate-procStart < 50);
        }

        postMessage({
            function: "finalizeLinkDiagram",
            arguments: [self.force_shadow, self.faces, self.components]
        });
    }
}

onmessage = function(e) {
    workerFunctions[e.data.function](...e.data.arguments);
}
