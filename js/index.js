let lab = new Lalolab();

function least_squares(X /* : Matrix */, Y /* : Matrix */) /* : least_squares */ {
    let betaHat = solve(mul(X, transpose(X)), mul(X, Y));

    return betaHat;
}

class MeshEdge {
    constructor(parent, start, stop) {
        this.parent = parent;
        this.start = start;
        this.stop = stop;

        this.svg = this.parent.edgeG.line(start.x(), start.y(), stop.x(), stop.y());
        this.svg.addClass("scaffold");
    }
}

class MeshNode {
    constructor(parent, x, y) {
        this.parent = parent;
        this.svg = this.parent.nodeG.circle(x, y, .3);
        this.svg.addClass('plnode');
    }


    set_r(r) {
        this.svg.attr({'r': r});
    }

    x() {
        return this.svg.attr('cx');
    }

    y() {
        return this.svg.attr('cy');
    }

    onClick(e) {

    }

}

class MeshDraw {
    constructor(div) {
        this.nodes = [];
        this.edges = {};

        this.draw = Snap(div);
        this.edgeG = this.draw.g();
        this.nodeG = this.draw.g();
        this.knotG = this.draw.g();
    }

    clear() {
        this.nodes = [];
        this.edges = {};

        this.edgeG.clear();
        this.nodeG.clear();
        this.knotG.clear();
    }

    set_embedding(g /*metric*/, embedding, map4v /*original map*/) {
        let points = embedding[0];
        let faces = embedding[1];
        let phi = embedding[2];
        /* set the embedding */
        this.g = g;

        // alter the viewbox to fit
        let min_x = min(get(points, range(), 0));
        let min_y = min(get(points, range(), 1));
        let max_x = max(get(points, range(), 0));
        let max_y = max(get(points, range(), 1));

        let wid = max_x-min_x;
        let hgt = max_y-min_y;

        let dx = wid*0.05;
        let dy = hgt*0.05;

        this.draw.attr({viewBox: [min_x-dx, min_y-dy, wid+2*dx, hgt+2*dy].join(",")});

        let i = 0;
        for (let point of points.toArray()) {
            let node = this.add_node(i, point[0], point[1]);
            node.set_r(this.g.gamma[this.nodes.length-1]);
            i += 1;
        }

        for (let edge of this.g.mesh.edges) {
            if (edge[0] in this.nodes && edge[1] in this.nodes) {
                let mesh_edge = this.add_edge(edge[0], edge[1]);

                if (Math.max(...edge) < map4v.nv+map4v.ne) {
                    // This edge has no face connection, and so is not scaffolding
                    mesh_edge.svg.addClass('edge');
                }
            }
        }

        let root_arc = map4v.arcs[0];
        let component = map4v.component(root_arc);

        // A path of the form anchor, control, anchor...
        let path = [];
        let pts = points.toArray();

        let last_ep = undefined;
        for (let arc of component) {
            let vi = arc.vert;
            let ei = map4v.nv+arc.edge;

            let vp = pts[vi];
            let ep = pts[ei];

            if (last_ep != undefined) {
                // We really want to push 'center' points as anchor points and
                // real points as control points
                let dp = sub(last_ep, vp);

                // Find center point by linear interpolation; whole way across is
                // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
                let whole_l = norm(dp);

                // We only want to go this.g.gamma[vi]/|dp| of the way:
                let mid_dp = mul(this.g.gamma[vi]/whole_l, dp);
                let mp = add(vp, mid_dp);

                path.push(mp);
                path.push(vp);
            }

            // We really want to push 'center' points as anchor points and
            // real points as control points
            let dp = sub(ep, vp); // ep - vp; notice vp + dp = ep

            // Find center point by linear interpolation; whole way across is
            // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
            let whole_l = norm(dp);

            // We only want to go this.g.gamma[vi]/|dp| of the way:
            let mid_dp = mul(this.g.gamma[vi]/whole_l, dp);
            let mp = add(vp, mid_dp);

            //path.push(vp);
            path.push(mp);
            path.push(ep);

            last_ep = ep;
        }
        let vi = root_arc.vert;
        let ei = map4v.nv+root_arc.edge;

        let vp = pts[vi];
        let ep = pts[ei];

        // We really want to push 'center' points as anchor points and
        // real points as control points
        let dp = sub(last_ep, vp);

        // Find center point by linear interpolation; whole way across is
        // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
        let whole_l = norm(dp);

        // We only want to go this.g.gamma[vi]/|dp| of the way:
        let mid_dp = mul(this.g.gamma[vi]/whole_l, dp);
        let mp = add(vp, mid_dp);

        path.push(mp);
        path.push(vp);

        path.push(path[0]);

        let pathStr = "M";
        pathStr += path[0].join(",");

        let idx = 0;
        for (let pt of path.slice(1)) {
            if (idx % 2 == 0) {
                pathStr += "Q";
            } else {
                pathStr += " ";
            }
            idx += 1;
            pathStr += pt.join(",");
        }

        //console.log(pathStr);
        let knot = this.knotG.path(pathStr+"Z");
        knot.addClass('knot');
    }

    add_node(i, x, y) {
        let node = new MeshNode(this, x, y);
        this.nodes[i] = node;
        return node;
    }

    add_edge(i, j) {
        let edge = new MeshEdge(this, this.nodes[i], this.nodes[j]);
        this.edges[[i,j]] = edge;
        this.edges[[j,i]] = edge;
        return edge;
    }
}

class Map4v {
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
            this.new_arc(2*ei);
            this.new_arc(2*ei+1);
            //console.log(this.arcs)

            this.set_edge(ei, [2*ei, 2*ei+1]);
        }

        // Set verts by looking through verts
        for (let vi in verts) {
            this.set_vert(vi, verts[vi]);
        }
    }

    triangulate(bdry_face_i=undefined) {
        let faces = this.faces();

        if (bdry_face_i === undefined) {
            let sizes = faces.map((f) => { return f.length; }, this);
            bdry_face_i = sizes.indexOf(max(sizes));
        }

        /* Triangles in the resulting triangulation will always have the form,
           [f, e, v]. We take advantage of this structure.

           This means that every triangle contains precisely one arc, and
           that every arc corresponds to precisely two resultant triangles---
           except in the case of the boundary face, whose triangles we discard */

        /* The edges in the resulting triangulation correspond either to arcs,
           which we want to display, and scaffolding. Every triangle consists
           of precisely one arc and two scaffold edges. */

        let triangles = [];

        /* Vertex indices correspond to components of the combinatorial map as follows:
           + [0:nv] correspond to actual vertices,
           + [nv:nv+ne] correspond to edge vertices,
           + [nv+ne:nv+ne+nf-1] correspond to (nonboundary) faces;

           Vertices are canonically indexed 0:nv+ne+nf-1. */
        let verts = range(0, this.nv+this.ne+faces.length-1);

        // Boundary verts are simply determined by the vertices and edges in the boundary face.
        //console.log(faces);
        let bdy_verts = (
            faces[bdry_face_i].map((a) => {return a.vert;}, this).concat(
                faces[bdry_face_i].map((a) => {return this.nv+a.edge;}, this))
        );
        // Remove the boundary face.
        faces.splice(bdry_face_i, 1);

        // Some edges correspond to arcs [vert->edge pairs]
        //console.log(this.arcs);
        let edges = this.arcs.map(
            (arc) => { return [arc.vert, this.nv+arc.edge]; }, this);

        // Loop through the faces, producing the rest of the edges and faces.
        for (let fi in faces) {
            let face = faces[fi];
            let fvi = this.nv+this.ne+parseInt(fi);
            // This face vertex has already been produced (nv+ne+fi)
            for (let arc of face) {
                // We get a scaffolding edge face->arc.edge
                edges.push([fvi, this.nv+arc.edge]);
                // We get a scaffolding edge face->arc.vert
                edges.push([fvi, arc.vert]);

                // We get two faces; one corresponding to this arc,
                triangles.push([fvi, this.nv+arc.edge, arc.vert]);
                // and one corresponding to its edge opposite
                let o_arc = this.edges[arc.edge][(arc.edgepos+1)%2];
                triangles.push([fvi, this.nv+o_arc.edge, o_arc.vert]);
            }
        }

        return [verts, bdy_verts, edges, triangles];
    }

    faces() {
        let left_arcs = new Set(this.arcs);
        let faces = [];

        while(left_arcs.size > 0) {
            let start_arc = Array.from(left_arcs).pop();

            let face = [];
            let arc = start_arc;
            let _failsafe = 0;
            do {
                left_arcs.delete(arc);
                face.push(arc);
                //console.log(arc);

                let o_arc = this.edges[arc.edge][(arc.edgepos+1)%2];
                arc = this.verts[o_arc.vert][(o_arc.vertpos+1)%4];
                //console.log(arc);
                _failsafe += 1;
                if (_failsafe > 500) {
                    console.log("Failure")
                    return faces;
                }
            } while (arc != start_arc)

            face.reverse();
            faces.push(face);
        }

        return faces;
    }

    component(arc) {
        let start_arc = arc;

        let component = [];
        do {
            component.push(arc);

            let o_arc = this.edges[arc.edge][(arc.edgepos+1)%2];
            arc = this.verts[o_arc.vert][(o_arc.vertpos+2)%4];
        } while (arc != start_arc)

        return component;
    }

    new_arc(idx) {
        this.arcs[idx] = {index: idx, edge:undefined, edgepos:undefined, vert:undefined, vertpos:undefined};
    }

    set_edge(idx, ais) {
        this.edges[idx] = ais.map((ai) => {return this.arcs[ai];}, this);

        for (let i in ais) {
            this.arcs[ais[i]].edge = idx;
            this.arcs[ais[i]].edgepos = parseInt(i);
        }
    }

    set_vert(idx, ais) {
        this.verts[idx] = ais.map((ai) => {return this.arcs[ai];}, this);

        for (let i in ais) {
            //console.log(i)
            //console.log(this.arcs)
            this.arcs[ais[i]].vert = parseInt(idx);
            this.arcs[ais[i]].vertpos = parseInt(i);
        }
    }
}

class TriangleMesh {
    constructor(verts, bdyverts, edges, faces) {
        this.verts = verts;
        this.bdyverts = bdyverts;
        this.edges = edges;
        this.faces = faces;
        this.bg_geom = "euclidean";
    }

    adjacent_edges(vert) {
        let adj = [];
        for (let edge of this.edges) {
            if (edge.indexOf(vert) >= 0) {
                adj.push(edge);
            }
        }
        return adj;
    }

    adjacent_faces(vert) {
        return this.faces.filter(function(f) {
            return f.indexOf(vert) >= 0;
        });
    }

    valence(vert) {
        return this.adjacent_faces(vert).length;
    }

    valences() {
        return this.verts.map(this.valence, this);
    }

    min_valence() {
        return min(this.valences());
    }

    chi() {
        // Euler characteristic
        return this.verts.length - this.edges.length + this.faces.length;
    }
}

function partition_face(face) {
    return [[face[0], [face[1], face[2]]],
            [face[1], [face[2], face[0]]],
            [face[2], [face[0], face[1]]]];
}

class DiscreteRiemannMetric {
    constructor(mesh, radius_map, edge_weights) {
        this.n = mesh.verts.length;
        this.mesh = mesh;

        this.gamma = array2vec(radius_map);
        this.u = this.conf_factor(radius_map);
        this.l = zeros(this.n, this.n);

        this.theta = {};
        this.phi = edge_weights;

        this.update();
    }

    static from_triangle_mesh(mesh) {
        let n = mesh.verts.length;
        let gamma = mesh.verts.map(function() { return 1; }, this);
        let phi = zeros(n, n);

        for (let edge of mesh.edges) {
            phi.val[edge[0]*phi.n + edge[1]] = 0;
            phi.val[edge[1]*phi.n + edge[0]] = 0;
        }

        return new DiscreteRiemannMetric(mesh, gamma, phi);
    }

    conf_factor(gamma) {
        return log(gamma);
    }

    compute_length(edge) {
        // console.log(this.gamma, this.gamma.val);
        let g0 = this.gamma[edge[0]];
        let g1 = this.gamma[edge[1]];

        // Law of cosines
        return sqrt(2*g0*g1*cos(this.phi.val[edge[0]*this.phi.n + edge[1]]) +
                    g0**2 + g1**2);
    }

    length(edge) {
        return this.l.val[edge[0]*this.l.n + edge[1]];
    }

    abc_for_vert(face, vert) {
        let other_vs = face.filter(function(v) { return v != vert; });
        let edge_a = [vert, other_vs[0]];
        let edge_b = [vert, other_vs[1]];
        let edge_c = other_vs;

        //console.log(face, vert, edge_a, edge_b, edge_c)
        return [edge_a, edge_b, edge_c].map(this.length, this);
    }

    compute_angle(face, vert) {
        let abc = this.abc_for_vert(face, vert);

        let ratio = (abc[0]**2 + abc[1]**2 - abc[2]**2)/(2.0*abc[0]*abc[1]);
        return Math.acos(ratio);
    }

    angle(face, vert) {
        let other_vs = face.filter(function(v) { return v != vert; });
        return this.theta[[vert, other_vs[0], other_vs[1]]];
    }

    curvature(vert) {
        let adj_faces = this.mesh.adjacent_faces(vert);
        if (this.mesh.bdyverts.indexOf(vert) >= 0) {
            return   Math.PI - sum(adj_faces.map((face) => {return this.angle(face, vert);}, this));
        } else {
            return 2*Math.PI - sum(adj_faces.map((face) => {return this.angle(face, vert);}, this));
        }
    }

    update() {
        this.gamma = exp(this.u);

        for (let edge of this.mesh.edges) {
            let l = this.compute_length(edge);
            this.l.val[edge[0]*this.l.n + edge[1]] = l;
            this.l.val[edge[1]*this.l.n + edge[0]] = l;
        }

        // Set angles using law of cosines
        for (let face of this.mesh.faces) {
            for (let part of partition_face(face)) {
                let theta = this.compute_angle(face, part[0]);
                this.theta[[part[0], part[1][0], part[1][1]]] = theta;
                this.theta[[part[0], part[1][1], part[1][0]]] = theta;
            }
        }

        // Set curvatures
        this.K = sparse(this.mesh.verts.map(this.curvature, this));
    }

    newton(target_K=null, dt=0.05, thresh=1e-4) {
        if (target_K == null) {

        }

        let g = new DiscreteRiemannMetric(this.mesh, this.gamma, this.phi);

        let K = g.K;
        let DeltaK = sub(target_K, K);

        let _failsafe = 0;
        while (max(abs(DeltaK)) > thresh){
            let H = this.hessian();
            let deltau = least_squares(H, DeltaK);

            g.u = sub(g.u, mul(dt, deltau));

            g.update();

            K = g.K;
            DeltaK = sub(target_K, K);

            //console.log(math.max(math.abs(DeltaK)));

            _failsafe += 1;
            if (_failsafe > 1000) {
                console.log("Took too long to flatten; abort!");
                return g;
            }
        }

        return g;
    }

    tau2(l_jk, g_j, g_k) {
        return .5*(l_jk**2 + g_j**2 - g_k**2);
    }

    face_area(face) {
        let gamma = this.theta[face];
        let a = this.length([face[0], face[1]]);
        let b = this.length([face[0], face[2]]);

        return .5*a*b*sin(gamma);
    }

    hessian() {
        let n = this.mesh.verts.length;
        let H = zeros(n, n);
        let t = this.tau2;

        for (let face of this.mesh.faces) {
            let i = face[0];
            let j = face[1];
            let k = face[2];

            let l_k = this.l.val[i*this.l.n + j];
            let l_i = this.l.val[j*this.l.n + k];
            let l_j = this.l.val[k*this.l.n + i];

            let g_i = this.gamma[i];
            let g_j = this.gamma[j];
            let g_k = this.gamma[k];

            let th_i = this.angle(face, i);
            let th_j = this.angle(face, j);
            let th_k = this.angle(face, k);

            let A = this.face_area(face);

            let L = diag([l_i, l_j, l_k]);
            let D = array2mat([[0,                t(l_i, g_j, g_k), t(l_i, g_k, g_j)],
                               [t(l_j, g_i, g_k), 0,                t(l_j, g_k, g_i)],
                               [t(l_k, g_i, g_j), t(l_k, g_j, g_i), 0               ]]);

            let Theta = cos(array2mat([
                [Math.PI, th_k, th_j],
                [th_k, Math.PI, th_i],
                [th_j, th_i, Math.PI]
            ]));

            let Tijk = mul(-.5/A, mul(mul(mul(L, Theta), inv(L)), D));

            for (let rowi of [0,1,2]) {
                let a = [i,j,k][rowi];
                for (let coli of [0,1,2]) {
                    let b = [i,j,k][coli];
                    H.val[a*H.n + b] += Tijk.val[rowi*Tijk.n + coli];
                }
            }
        }
        return H;
    }
}

function adjacent_faces_and_edge(faces, face) {
    let adj = [];

    for (let aface of faces) {
        let e = aface.filter((v) => { return face.indexOf(v) >= 0; });
        if (e.length == 2) {
            if (face[(face.indexOf(e[0])+1)%3] == e[1]) {
                adj.push([aface, [e[1], e[0]]]);
            } else {
                adj.push([aface, e]);
            }
        }
    }

    return adj;
}

function adj_or_faces_and_edge(faces, face) {
    let adj = [];

    for (let aface of faces) {
        if (aface.indexOf(face[0]) >= 0 && aface.indexOf(face[1]) >= 0) {
            adj.push([aface, [face[0], face[1]]]);
        } else if (aface.indexOf(face[1]) >= 0 && aface.indexOf(face[2]) >= 0) {
            adj.push([aface, [face[1], face[2]]]);
        } else if (aface.indexOf(face[2]) >= 0 && aface.indexOf(face[0]) >= 0) {
            adj.push([aface, [face[2], face[0]]]);
        }
    }

    return adj;
}

function orient_faces(faces) {
    let edges = [];
    let oriented = [];
    let to_orient = faces.slice();
    let adj_queue = new Set();

    let f_0 = to_orient.pop();
    edges.concat([
        [f_0[0], f_0[1]],
        [f_0[1], f_0[2]],
        [f_0[2], f_0[0]]]);
    oriented.push(f_0);

    for (let adj_pair of adjacent_faces_and_edge(to_orient, f_0)) {
        adj_queue.add(adj_pair);
    }

    let _failsafe = 0;
    while (to_orient.length > 0) {
        //console.log(adj_queue, to_orient)
        let adj_pair = Array.from(adj_queue).pop();
        adj_queue.delete(adj_pair);

        //console.log(e);
        let F = adj_pair[0];
        let e = adj_pair[1];

        _failsafe += 1;
        if (_failsafe > 5000) {
            console.log("Too much time spent orienting faces...");
            return oriented;
        }

        if (to_orient.indexOf(F) < 0) {
            continue;
        }

        let v_i, v_j, v_k;
        if (e.indexOf(F[0]) < 0) {
            v_k = F[0];
        } else if (e.indexOf(F[1]) < 0) {
            v_k = F[1];
        } else {
            v_k = F[2];
        }
        v_i = e[0];
        v_j = e[1];

        let i,j,k;
        if (edges.includes([v_i, v_j])) {
            i = v_j; j = v_i; k = v_k;
        } else {
            i = v_i; j = v_j; k = v_k;
        }

        edges.concat([[i,j], [j,k], [k,i]]);
        oriented.push([i,j,k]);

        to_orient.splice(to_orient.indexOf(F), 1);
        for (let adj_pair of adjacent_faces_and_edge(to_orient, [i,j,k])) {
            adj_queue.add(adj_pair);
        }

    }

    return oriented;
}

function u_theta(theta) {
    return [cos(theta), sin(theta)];
}

function embed_faces(g) {
    let mesh = g.mesh;
    let x = mul(-30, ones(g.n, 2));
    let phi = {};
    let pi = Math.PI;

    let faces = orient_faces(mesh.faces);
    let to_embed = faces.slice();
    let embed_queue = new Set();

    let f_0 = to_embed.pop();
    let i = f_0[0], j = f_0[1], k = f_0[2];

    x.val[i*x.n + 0] = 0;
    x.val[i*x.n + 1] = 0;
    x.val[j*x.n + 0] = g.length([i,j]);
    x.val[j*x.n + 1] = 0;

    let phi_ik = g.angle(f_0, i) % (2*pi);
    let g_ik = g.length([i,k]);

    x.val[k*x.n + 0] = g_ik*cos(phi_ik);
    x.val[k*x.n + 1] = g_ik*sin(phi_ik);

    let phi_jk = (pi - g.angle(f_0, j)) % (2*pi);

    phi[[i,j]] = 0;
    phi[[j,i]] = pi;
    phi[[j,k]] = phi_jk;
    phi[[k,j]] = pi+phi_jk;
    phi[[i,k]] = phi_ik;
    phi[[k,i]] = pi+phi_ik;

    for (let adj_pair of adj_or_faces_and_edge(to_embed, f_0)) {
        embed_queue.add(adj_pair);
    }

    let _failsafe = 0;
    while (to_embed.length > 0) {
        let adj_pair = Array.from(embed_queue).pop();
        embed_queue.delete(adj_pair);
        let F = adj_pair[0];
        let e = adj_pair[1];

        _failsafe += 1;
        if (_failsafe > 5000) {
            //console.log(F, e);
            console.log("Too much time spent embedding faces...");
            break;
        }

        if (to_embed.indexOf(F) < 0) {
            continue;
        }

        //if (_failsafe == 11) {
        // F = [1,6,5];
        //}

        //console.log("Iteration number "+_failsafe+":")

        let i = e[0], j = e[1], k;
        if (F[(F.indexOf(i)+1)%3] != j) {
            k = F[(F.indexOf(i)+1)%3];
            j = e[0], i = e[1];
        } else {
            k = F[(F.indexOf(i)+2)%3]
        }

        //console.log(F, e, i, j, k);
        //console.log(x.subset(math.index([i,j,k], [0,1])).toArray());

        // We already know x[i], x[j]. We only have to find x[k].
        phi_ik = (phi[[i,j]] + g.angle(F, i)) % (2*pi);
        phi_jk = (phi[[j,i]] - g.angle(F, j)) % (2*pi);

        //if (_failsafe == 11) {
        // phi_ik = (phi[[i,j]] - g.angle(F, i)) % (2*pi);
        //}
        if (!([j,k] in phi) && !([i,k] in phi)) {
            g_ik = g.length([i,k]);
            x.val[k*x.n + 0] = x.val[i*x.n + 0] + g_ik*cos(phi_ik);
            x.val[k*x.n + 1] = x.val[i*x.n + 1] + g_ik*sin(phi_ik);
        }

        //console.log(g.length([i,k]))
        //console.log(g.angle(F, i));
        //console.log(phi[[i,j]]);
        //console.log(x.subset(math.index(i, [0,1])).toArray()[0])
        //console.log(phi_ik)
        //console.log(u_theta(phi_ik))

        if (!([j,k] in phi)) {
            phi[[j,k]] = phi_jk;
            phi[[k,j]] = pi+phi_jk;
        }

        if (!([i,k] in phi)) {
            phi[[i,k]] = phi_ik;
            phi[[k,i]] = pi+phi_ik;
        }

        to_embed.splice(to_embed.indexOf(F), 1);
        for (let adj_pair of adj_or_faces_and_edge(to_embed, F)) {
            embed_queue.add(adj_pair);
        }

    }

    return [x, faces, phi];
}

function map_to_mesh(m) {
    /* take a combinatorial map and produce a triangular mesh
       with a vertex for each crossing, face, and edge

       also returns a mapping from vertices to c,e,f's and edges to edges

       edges are canonically of the form [2i, 2i+1],
       input is an array of vertices of the form [a, b, c, d]

       returns [mesh, mapping]*/
}

let meshDraw = new MeshDraw("#knot-draw");

function drawMap(sigma, cross_bend=8) {
    var tstart = Date.now();
    meshDraw.clear();

    let m4v = new Map4v(sigma);

    let triangulation = m4v.triangulate();

    //console.log(triangulation);

    let testMesh = new TriangleMesh(
        triangulation[0], triangulation[1], triangulation[2], triangulation[3]
    );

    let cpmetric = DiscreteRiemannMetric.from_triangle_mesh(testMesh);

    let K = zeros(testMesh.verts.length);
    //K.subset(math.index(testMesh.bdyverts, 0),
    //         testMesh.bdyverts.map(function(b) { return 2*Math.PI/testMesh.bdyverts.length; }));

    // Set boundary crossing target curvature
    let bdyCross = testMesh.bdyverts.filter(function(vi) { return vi < m4v.nv; });
    let bdyEdge = testMesh.bdyverts.filter(function(vi) { return (vi >= m4v.nv && vi < m4v.nv+m4v.ne); });
    let fac = 8; // Inverse of how concave crossing vertices should be imbedded
    for (let bci of bdyCross) {
        K[bci] = -Math.PI/fac;
    }
    for (let bei of bdyEdge) {
        K[bei] = (2*Math.PI + bdyCross.length*Math.PI/fac)/bdyEdge.length;
    }

    let flat_poly = cpmetric.newton(K, 1, 5e-2);

    let embedding = embed_faces(flat_poly);

    meshDraw.set_embedding(flat_poly, embedding, m4v);
    console.log("Computation finished in: " + (Date.now() - tstart) + " milliseconds");
}

// Trefoil
//let sigma = [[0, 6, 11, 5], [1, 8, 2, 7], [3, 9, 4, 10]];

// More complicated
//let sigma = [[0, 41, 99, 42], [1, 47, 2, 48], [3, 34, 4, 33], [5, 96, 6, 95], [7, 13, 8, 14], [9, 68, 10, 67], [11, 69, 12, 70], [15, 62, 16, 61], [17, 59, 18, 60], [19, 85, 20, 86], [21, 27, 22, 28], [23, 82, 24, 81], [25, 83, 26, 84], [29, 88, 30, 87], [31, 93, 32, 94], [35, 46, 36, 45], [37, 43, 38, 44], [39, 98, 40, 97], [49, 76, 50, 75], [51, 73, 52, 74], [53, 72, 54, 71], [55, 65, 56, 66], [57, 64, 58, 63], [77, 92, 78, 91], [79, 89, 80, 90]];

// Even moreso
let sigma = [[1, 287, 2, 288], [0, 289, 299, 290], [3, 126, 4, 125], [5, 16, 6, 15], [7, 13, 8, 14], [9, 255, 10, 256], [11, 254, 12, 253], [17, 127, 18, 128], [19, 282, 20, 281], [21, 132, 22, 131], [23, 277, 24, 278], [25, 235, 26, 236], [27, 222, 28, 221], [29, 223, 30, 224], [31, 38, 32, 37], [33, 228, 34, 227], [35, 225, 36, 226], [39, 230, 40, 229], [41, 156, 42, 155], [43, 157, 44, 158], [45, 164, 46, 163], [47, 165, 48, 166], [49, 56, 50, 55], [51, 198, 52, 197], [53, 195, 54, 196], [57, 243, 58, 244], [59, 242, 60, 241], [61, 231, 62, 232], [63, 234, 64, 233], [65, 248, 66, 247], [67, 249, 68, 250], [69, 276, 70, 275], [71, 133, 72, 134], [73, 83, 74, 84], [75, 82, 76, 81], [77, 283, 78, 284], [79, 286, 80, 285], [85, 295, 86, 296], [87, 294, 88, 293], [89, 291, 90, 292], [91, 298, 92, 297], [93, 136, 94, 135], [95, 109, 96, 110], [97, 108, 98, 107], [99, 122, 100, 121], [101, 259, 102, 260], [103, 262, 104, 261], [105, 119, 106, 120], [111, 273, 112, 274], [113, 272, 114, 271], [115, 265, 116, 266], [117, 264, 118, 263], [123, 137, 124, 138], [129, 280, 130, 279], [139, 257, 140, 258], [141, 268, 142, 267], [143, 269, 144, 270], [145, 252, 146, 251], [147, 237, 148, 238], [149, 216, 150, 215], [151, 217, 152, 218], [153, 220, 154, 219], [159, 210, 160, 209], [161, 211, 162, 212], [167, 190, 168, 189], [169, 187, 170, 188], [171, 205, 172, 206], [173, 200, 174, 199], [175, 201, 176, 202], [177, 204, 178, 203], [179, 186, 180, 185], [181, 191, 182, 192], [183, 194, 184, 193], [207, 214, 208, 213], [239, 245, 240, 246]];

drawMap(sigma);

document.getElementById("map_submit").onclick = function(ev) {
    document.querySelector(ev.target.dataset.errbox).textContent = "";
    try {
        drawMap(JSON.parse(document.querySelector(ev.target.dataset.target).value));
    } catch(err) {
        document.querySelector(ev.target.dataset.errbox).textContent = "Error";
        console.log("Error:", err);
    }
};
