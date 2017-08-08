export default class ForceLinkDiagram {
    /* Link diagram embedding improved by ImPrEd */
    constructor (verts, edges, faces, components) {
        this.verts = verts;
        this.edges = edges;
        this.edgeMatrix = [];
        for (let edge of this.edges) {
            if (!(edge[0] in this.edgeMatrix)) {
                this.edgeMatrix[edge[0]] = [];
            }
            if (!(edge[1] in this.edgeMatrix)) {
                this.edgeMatrix[edge[1]] = [];
            }
            //console.log(this.edgeMatrix);
            this.edgeMatrix[edge[0]][edge[1]] = edge;
            this.edgeMatrix[edge[1]][edge[0]] = edge;
        }
        this.faces = faces;
        this.components = components;

        this.freeVi = [];

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

        let N = this.verts.length;
        function index(edge) {
            return edge[0] * N + edge[1];
        }
        function deindex(index) {
            return [(index - index%N) / N, index % N];
        }
        this.paths = [];
        let edgeSearch = new Set(this.edges.map(e => index(e)));
        while (edgeSearch.size > 0) {
            let path;
            let edge = deindex(edgeSearch.values().next().value);
            edgeSearch.delete(index(edge));
            //console.log(edge);

            let [prev, next] = edge;
            path = [prev, next];
            //console.log(path);
            while (this.adjMap[prev].length == 2) {
                //console.log(this.adjMap[prev]);
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
                //console.log(this.adjMap[next]);
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
        //console.log(this.paths);

        this.delta = 2;
        this.gamma = 5;

        this.dbar = 3*this.delta;

        this.alpha = 2*this.delta;
        this.beta = 2*this.delta;

        //this.flexEvery = 2;

        this.aExp = 1;
        this.erExp = 2;

        this.calculateSurroundingEdges();

        this.numIter = 0;
    }

    newVertIdx() {
        if (this.freeVi.length > 0) {
            return this.freeVi.pop();
        }
        return this.verts.length;
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
        let edges = new Set();
        for (let face of this.faces) {
            if (face.includes(ui)) {
                for (let i = 0; i < face.length-1; i++) {
                    console.assert(this.edges.filter(
                        e => ((e[0] == face[i] && e[1] == face[i+1]) ||
                              (e[1] == face[i] && e[0] == face[i+1]))).length > 0);
                    edges.add(this.edgeMatrix[face[i]][face[i+1]]);
                    console.assert(!edges.has(undefined));
                }
                edges.add(this.edgeMatrix[face[face.length-1]][face[0]]);
                console.assert(!edges.has(undefined));
            }
        }
        console.assert(!edges.has(undefined));
        return edges;
    }

    calculateSurroundingEdges() {
        this.surrEdges = [];
        for (let i = 0; i < this.verts.length; i++) {
            this.surrEdges[i] = this.surroundingEdges(i);
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

    move (ui, FUx, FUy, MU) {
        console.assert (this.verts[ui] !== undefined, this.verts);
        let i = this.octant(FUx, FUy);

        let FU = [FUx, FUy];

        if (ui == 2) {
            console.log(ui, FUx, FUy, MU[i]);
        }

        MU[i] = Math.max(0, MU[i]);
        let fU = norm(FU);
        let du;
        if (fU <= MU[i]) {
            du = FU;
        } else {
            du = mul(MU[i]/fU, FU);
        }

        if (ui == 2) {
            console.log(ui, FUx, FUy, MU[i], du, norm(du), Number.EPSILON);
        }

        if (norm(du) < 10*Number.EPSILON) { return; }

        this.verts[ui][0] += du[0];
        this.verts[ui][1] += du[1];
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

        if (this.numIter == 4) {
            if (vi == 13) {
                //console.log(a, b, c, p);
                //console.log(u, v, 1-(u+v));
            }
        }

        return (u >= 0) && (v >= 0) && (u+v < 1);
    }

    contract() {
        for (let path of this.paths) {
            for (let i = 0; i < path.length-2; i++) {
                let [ai, bi, ci] = path.slice(i, i+3);
                let [a, b, c] = [this.verts[ai], this.verts[bi], this.verts[ci]];
                if (norm(sub(a, b)) < this.alpha) {
                    if (this.faces.filter(f => f.includes(ai)).some(f => f.length <= 3)) { continue; }

                    if (!this.verts.some(
                        (v, vi) => (vi != ai && vi != bi && vi != ci && this.triangleIncludes(a,b,c,v,vi)))) {
                        // contract

                        //console.log("deleting vertex ", bi);

                        // delete the vertex at bi
                        delete this.verts[bi];
                        delete this.adjMap[bi];

                        // remove bi from its adjacent vertices at ai and ci
                        this.adjMap[ai].splice(this.adjMap[ai].indexOf(bi), 1);
                        this.adjMap[ci].splice(this.adjMap[ci].indexOf(bi), 1);
                        this.adjMap[ai].push(ci);
                        this.adjMap[ci].push(ai);
                        let nEdge = [ai, ci];
                        this.edgeMatrix[ai][ci] = nEdge;
                        this.edgeMatrix[ci][ai] = nEdge;
                        this.edges.push(nEdge);

                        // remove bi from any faces
                        let bFaces = this.faces.filter(f => f.includes(bi));
                        // first, update these faces' surrEdges
                        for (let f of bFaces) {
                            for (let ui of f) {
                                for (let bEdge of this.edgeMatrix[bi]) {
                                    if (bEdge === undefined) { continue; }
                                    //console.log(bEdge);
                                    this.surrEdges[ui].delete(bEdge);
                                }
                                this.surrEdges[ui].add(nEdge);
                            }
                        }

                        // finally, remove bi from the faces
                        //console.log(bFaces);
                        bFaces.forEach(f => f.splice(f.indexOf(bi), 1));
                        //console.log("bFaces has no", bi, bFaces);

                        // remove bi from any components.
                        this.components.filter(cp => cp.includes(bi)).forEach(
                            cp => {while(cp.includes(bi)) { cp.splice(cp.indexOf(bi), 1); }}
                        );

                        // remove edges including bi
                        this.edges.splice(this.edges.indexOf(this.edgeMatrix[ai][bi]), 1);
                        this.edges.splice(this.edges.indexOf(this.edgeMatrix[bi][ci]), 1);
                        delete this.edgeMatrix[ai][bi];
                        delete this.edgeMatrix[ci][bi];
                        delete this.edgeMatrix[bi];

                        // remove bi from this path
                        path.splice(path.indexOf(bi), 1);
                        break; // Only contract one bend of a path at a time???
                    }
                }
            }
            //break;
        }
    }

    expand() {
        for (let path of this.paths) {
            for (let i = 0; i < path.length-1; i++) {
                let [ai, bi] = path.slice(i,i+2);
                let [a, b] = [this.verts[ai], this.verts[bi]];
                let elen = norm(sub(a, b);
                if (elen) > this.beta) {
                    // Insert a new vertex along this path
                    //let nv = 
                }

            }
        }
    }

    update() {
        // if (contract) {
        //console.log(this.surrEdges);

        if (true) {
            this.contract();
        }

        if (true) {
            this.expand();
        }
        //console.log(this.edges);
        //console.log(this.faces);
        //console.log(this.surrEdges);

        //return;
        if (this.numIter == 4) {
           // return;
        }
        this.numIter++;

        //console.log(this.adjMap);
        let FX = zeros(this.verts.length);
        let FY = zeros(this.verts.length);
        let M = [];
        for (let i = 0; i < this.verts.length; i++) {
            M.push([this.dbar, this.dbar, this.dbar, this.dbar,
                    this.dbar, this.dbar, this.dbar, this.dbar]);
        }

        let barycenter = mul(1/this.verts.length, sum(this.verts, 2));

        //for (let ui = 0; ui < this.verts.length; ui++) {
        for (let ui in this.verts) {
            
            // Calculate gravity force
            let db = sub(barycenter, this.verts[ui]);
            let nDb = norm(db);
            FX[ui] += db[0]/nDb;
            FY[ui] += db[1]/nDb;

            // Calculate total node-node repulsive force
            //for (let vi = 0; vi < this.verts.length; vi++) {
            for (let vi in this.verts) {
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

                } else {
                    let va = sub(this.verts[ai], this.verts[ui]);
                    let vb = sub(this.verts[bi], this.verts[ui]);
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

                for (let j = 0; j < MU.length; j++) {
                    if ((ell-j+8)%8 == 0) {
                        M[ai][j] = Math.min(M[ai][j], maxR);
                    } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                        M[ai][j] = Math.min(M[ai][j], maxR /
                                            Math.cos(cw_angle - (j+1)*Math.PI/4));
                    } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                        M[ai][j] = Math.min(M[ai][j], maxR /
                                            Math.cos(cw_angle - (j)*Math.PI/4));
                    }
                }

                for (let j = 0; j < MU.length; j++) {
                    if ((ell-j+8)%8 == 0) {
                        M[bi][j] = Math.min(M[bi][j], maxR);
                    } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                        M[bi][j] = Math.min(M[bi][j], maxR /
                                            Math.cos(cw_angle - (j+1)*Math.PI/4));
                    } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                        M[bi][j] = Math.min(M[bi][j], maxR /
                                            Math.cos(cw_angle - (j)*Math.PI/4));
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
