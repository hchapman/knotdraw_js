export default class ForceLinkDiagram {
    /* Link diagram embedding improved by ImPrEd */
    constructor (verts, edges, faces, components) {
        this.verts = new Map(verts.map((v, vi) => [vi, v]));
        this.edges = edges;
        this.edgeMatrix = [];
        for (let edge of this.edges) {
            if (!(edge[0] in this.edgeMatrix)) {
                this.edgeMatrix[edge[0]] = [];
            }
            if (!(edge[1] in this.edgeMatrix)) {
                this.edgeMatrix[edge[1]] = [];
            }
            this.edgeMatrix[edge[0]][edge[1]] = edge;
            this.edgeMatrix[edge[1]][edge[0]] = edge;
        }
        this.faces = faces;
        this.components = components;

        this.freeVi = [];
        this.maxVi = Math.max(...verts.keys())+1;

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
        while (edgeSearch.size > 0) {
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

        this.delta = 2;
        this.gamma = 5;

        this.dbar = 3*this.delta;

        this.alpha = 2*this.delta;
        this.beta = 3*this.delta;

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
        return this.verts.size;
    }

    distance(u, v) {
        return norm(sub(u, v));
    }

    forceAvert(u, v) {
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
                    console.assert(this.edges.some(
                        e => ((e[0] == face[i] && e[1] == face[i+1]) ||
                              (e[1] == face[i] && e[0] == face[i+1]))), this.edges, face, i, ui);
                    if (face[i] != ui && face[i+1] != ui) {
                        edges.add(this.edgeMatrix[face[i]][face[i+1]]);
                    }
                    console.assert(!edges.has(undefined));
                }
                if (face[face.length-1] != ui && face[0] != ui) {
                    edges.add(this.edgeMatrix[face[face.length-1]][face[0]]);
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

        this.verts.get(ui)[0] += du[0];
        this.verts.get(ui)[1] += du[1];
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
                    let [a, b, c] = [this.verts.get(ai), this.verts.get(bi), this.verts.get(ci)];
                    if (norm(sub(a, b)) < this.alpha) {
                        if (this.faces.filter(f => f.includes(ai)).some(f => f.length <= 3)) { continue; }

                        if (!vArray.some(
                            ([vi, v]) => (vi != ai && vi != bi && vi != ci && this.triangleIncludes(a,b,c,v,vi)))) {
                            // contract

                            // delete the vertex at bi
                            this.verts.delete(bi);
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
                            this.edges.splice(this.edges.indexOf(this.edgeMatrix[ai][bi]), 1);
                            this.edges.splice(this.edges.indexOf(this.edgeMatrix[bi][ci]), 1);
                            delete this.edgeMatrix[ai][bi];
                            delete this.edgeMatrix[ci][bi];
                            delete this.edgeMatrix[bi];

                            // remove bi from this path
                            path.splice(path.indexOf(bi), 1);
                            reContract = true;
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
                let [a, b] = [this.verts.get(ai), this.verts.get(bi)];
                let elen = norm(sub(a, b));
                if (elen > this.beta) {
                    // Insert a new vertex along this path
                    let vi = this.maxVi++;
                    this.verts.set(vi, mul(.5, add(a, b)));

                    // Edge to delete; which was split
                    let splitEdge = this.edgeMatrix[ai][bi];

                    // Add the new edges to the diagram
                    let nEdgeA = [ai, vi];
                    let nEdgeB = [vi, bi];

                    // Update the edgematrix of edges
                    this.edgeMatrix[vi] = [];
                    this.edgeMatrix[ai][vi] = nEdgeA;
                    this.edgeMatrix[vi][ai] = nEdgeA;
                    this.edgeMatrix[bi][vi] = nEdgeB;
                    this.edgeMatrix[vi][bi] = nEdgeB;

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

        // if (contract) {
        if (true) {
            this.contract();
        }

        // if (expand) {
        if (true) {
            this.expand();
        }

        this.numIter++;

        let F = new Map(Array.from(this.verts.keys()).map(k => [k, [0,0]])); //zeros(this.verts.size);
        let M = new Map();
        for (let [vi, v] of this.verts) {
            M.set(vi, [this.dbar, this.dbar, this.dbar, this.dbar,
                       this.dbar, this.dbar, this.dbar, this.dbar]);
        }

        let barycenter = mul(1/this.verts.size, sum(Array.from(this.verts.values()), 2));

        //for (let ui = 0; ui < this.verts.size; ui++) {
        for (let [ui, u] of this.verts) {
            let Fu = F.get(ui);

            // Calculate gravity force
            let db = sub(barycenter, this.verts.get(ui));
            let nDb = norm(db);
            Fu[0] += db[0]/nDb;
            Fu[1] += db[1]/nDb;

            // Calculate total node-node repulsive force
            for (let [vi, v] of this.verts) {
                if (ui != vi) {
                    if (this.paths.some(p => p.includes(ui) && p.includes(vi))) {
                        //console.log("consec");
                        continue;
                    }
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
            //for (let vi of this.paths.find(p => p.includes(ui)).filter(ell => ell != ui)) {
                let dF = this.forceAvert(this.verts.get(ui), this.verts.get(vi));

                Fu[0] += dF[0];
                Fu[1] += dF[1];
            }

            // calculate node-edge repulsive force
            for (let edge of this.surrEdges.get(ui)) {
                let [ai, bi] = edge;
                if (ui == ai || ui == bi) {
                    continue;
                }
                let ve = this.computeVe(
                    u, this.verts.get(ai), this.verts.get(bi));

                if (this.veOnEdge(ve, this.verts.get(ai), this.verts.get(bi))) {
                    let dF = this.forceRedge(
                        u, this.verts.get(ai), this.verts.get(bi), ve);
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
                let ve = this.computeVe(
                    u, this.verts.get(ai), this.verts.get(bi));

                let cv;

                if (this.veOnEdge(ve, this.verts.get(ai), this.verts.get(bi))) {
                    cv = sub(ve, u);

                } else {
                    let va = sub(this.verts.get(ai), u);
                    let vb = sub(this.verts.get(bi), u);
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
