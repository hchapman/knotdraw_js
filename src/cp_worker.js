importScripts('lalolib/lalolib.js');
importScripts('planarmap-em/libplanarmap-em.js');

function least_squares(X /* : Matrix */, Y /* : Matrix */) /* : least_squares */ {
    //console.log("X", X, inv(X), det(X), qr(X, true), Y);
    let betaHat = solve(mul(X, transpose(X)), mul(X, Y));

    return betaHat;
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
            this.new_arc(2*ei);
            this.new_arc(2*ei+1);
            //console.log(this.arcs)

            this.set_edge(ei, [2*ei, 2*ei+1]);
        }

        // Set verts by looking through verts
        for (let vi in verts) {
            this.set_vert(vi, verts[vi]);
        }

        this.generate_faces();
        this.generate_components();
    }

    triangulate(bdry_face_i=undefined) {
        if (bdry_face_i === undefined) {
            let sizes = this.faces.map((f) => f.length);
            bdry_face_i = sizes.indexOf(max(sizes));
        }

        /* Allowing for triangulations of maps with isthmi require that some
         * triangles are not necessarily of the form f,e,v and for some edges
         * (i.e. monogon edges) to have multiple triangulation vertices in order
         * to avoid singular flat embeddings.
         *
         * This also complicates the old idea that arcs mapped injectively into
         * triangulation edges---so it might be easiest to keep a new data
         * structure with component mappings instead.
         */

        let triangles = [];

        // We will now build the array of triangulation verts dynamically;
        // furthermore we will make them data objects which hold references to
        // appropriate objects
        let verts = [];

        // We will build edges of the triangulation as we process, too.
        let edges = [];

        // We will now calculate boundary vertices later; external faces with,
        // say, isthmi will require scaffolding which makes some vertices which
        // used to be boundary vertices internal.
        let bdy_face = this.faces[bdry_face_i];
        let bdy_verts;

        /* Since we want to be mindful of the paths---the components---when we
         * draw the diagram, we will create the triangulation as follows: We
         * will first run through each component, creating verts and edges as we
         * go. We will then consider faces (and the boundary face). */
        let tri_map = {
            verts: [],
            edges: [],
            faces: [],
            regions: [],
            comps: [],
            arcs : []
        };

        //console.log(this.faces);

        for (let component of this.components) {
            let tri_comp = [];
            for (let arc of component) {
                // Add a new vert for the out_vert of this arc, unless we've already
                // hit this vertex before.
                if (tri_map.verts[this.out_vert_i(arc)] === undefined) {
                    tri_map.verts[this.out_vert_i(arc)] = (verts.length);
                    verts.push(verts.length);
                }
                // Regardless, add this vert to the component
                tri_comp.push(tri_map.verts[this.out_vert_i(arc)]);

                // Add a new vert for the edge of this arc
                //console.log(this.faces[arc.face].length);
                //if (this.faces[arc.face].length <= 2 ||
                //    this.faces[(this.edges[arc.edge][(arc.edgepos+1)%2]).face].length <= 2) {
                    tri_map.edges[arc.edge] = [verts.length];
                    tri_comp.push(verts.length);
                    verts.push(verts.length);
                //} else {
                //    tri_map.edges[arc.edge] = [];
                //}

                if (this.out_vert_i(arc) == this.in_vert_i(arc)) {
                    // This arc corresponds to a monogon, and so we must add two
                    // edge vertices rather than one
                    tri_map.edges[arc.edge].push(verts.length);
                    tri_comp.push(verts.length);
                    verts.push(verts.length);
                }

                // Add this edge to tri_map.arcs, which contains direction information
                tri_map.arcs[arc.index] = tri_map.edges[arc.edge];
                tri_map.arcs[this.edges[arc.edge][(arc.edgepos+1)%2].index
                            ] = tri_map.edges[arc.edge].slice().reverse();
            }

            // Push this component on
            tri_map.comps.push(tri_comp);

            // Add triangulation edges for this component
            for (let cvi = 0; cvi < tri_comp.length-1; cvi++) {
                edges.push([tri_comp[cvi], tri_comp[cvi+1]]);
            }
            edges.push([tri_comp[tri_comp.length-1], tri_comp[0]]);
        }

        // Every vertex now has a corresponding tri vertex
        // Every edge now has a corresponding tri vertex (or two)
        // We must now go through and process the faces.
        for (let fi in this.faces) {
            // Regardless, we must identify what isthmi this face has
            // isthmus_arcs will contain one arc for each vertex---the other arc
            // will necessarily be vertex-opposite.
            let isthmus_arcs = [];
            let isthmus_vert = [];

            let _seen_vi = [];
            for (let arc of this.faces[fi]) {
                if (_seen_vi.includes(arc.vert)) {
                    isthmus_arcs.push(arc);
                    isthmus_vert.push(arc.vert);
                } else {
                    _seen_vi.push(arc.vert);
                }
            }

            // Independent of whether this face is a boundary face, we must
            // add scaffolding for any isthmi, if there are any.
            for (let arc of isthmus_arcs) {
                // Add scaffolding for this arc
                let pre_arc = this.verts[arc.vert][(arc.vertpos+3)%4];
                //console.log("!!", arc, tri_map.edges[arc.edge], tri_map.arcs[arc.index]);
                //console.log(pre_arc);
                edges.push([tri_map.arcs[arc.index][0],
                            tri_map.arcs[pre_arc.index][0]]);

                triangles.push([tri_map.arcs[arc.index][0],
                                tri_map.arcs[pre_arc.index][0],
                                tri_map.verts[arc.vert]]);

                // Add scaffolding for opposite arc
                let o_arc = this.verts[arc.vert][(arc.vertpos+2)%4];
                //console.log("~~", o_arc);
                pre_arc = this.verts[o_arc.vert][(o_arc.vertpos+3)%4];
                edges.push([tri_map.arcs[o_arc.index][0],
                            tri_map.arcs[pre_arc.index][0]]);

                triangles.push([tri_map.arcs[o_arc.index][0],
                                tri_map.arcs[pre_arc.index][0],
                                tri_map.verts[o_arc.vert]]);
            }

            // Finally, get a list of tri verts around this face.
            let tri_face = [];
            let region = [];
            for (let arc of this.faces[fi]) {
                if (!isthmus_vert.includes(arc.vert)) {
                    // If this vertex is not an isthmus, then it'll be in the
                    // face
                    tri_face.push(tri_map.verts[arc.vert]);
                }

                // Region wants isthmus to be included regardless
                region.push(tri_map.verts[arc.vert]);
                //console.log("{", region);

                // This arcs triangulation vertices need to be added, too
                tri_face.splice(tri_face.length, 0, ...tri_map.arcs[arc.index]);
                region.splice(region.length, 0, ...tri_map.arcs[arc.index]);
                //console.log(region, "}")
            }

            tri_map.regions[fi] = region;
            if (region.includes(63)) {
                //console.log(region, "&&&&");
            }

            //console.log(tri_face);

            if (fi == bdry_face_i) {
                // This face is the boundary face and must be processed
                tri_map.faces[fi] = [];
                bdy_verts = tri_face;

            } else {
                // This face is an internal face. As we triangulate, each face
                // has precisely one vertex, although through R-moves it is
                // possible that faces will ultimately consist of several
                // vertices.
                tri_map.faces[fi] = [verts.length];

                // For each vertex in this face we add an edge, and a face
                edges.push([verts.length, tri_face[0]]);
                for (let vi = 1; vi < tri_face.length; vi++) {
                    edges.push([verts.length, tri_face[vi]]);
                    triangles.push([verts.length, tri_face[vi-1], tri_face[vi]]);
                }
                triangles.push([verts.length, tri_face[tri_face.length-1], tri_face[0]]);

                // Finally, add this vertex to verts
                verts.push(verts.length);
            }
        }

        return [verts, bdy_verts, edges, triangles, tri_map];
    }

    old_triangulate(bdry_face_i=undefined) {
        let faces = this.faces;

        if (bdry_face_i === undefined) {
            let sizes = faces.map((f) => f.length, this);
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

    generate_faces() {
        let left_arcs = new Set(this.arcs);
        this.faces = [];

        while(left_arcs.size > 0) {
            let start_arc = Array.from(left_arcs).pop();

            let face = [];
            let arc = start_arc;
            let _failsafe = 0;
            do {
                left_arcs.delete(arc);
                face.push(arc);
                arc.face = this.faces.length;
                //console.log(arc);

                let o_arc = this.edges[arc.edge][(arc.edgepos+1)%2];
                arc = this.verts[o_arc.vert][(o_arc.vertpos+1)%4];
                //console.log(arc);
                _failsafe += 1;
                if (_failsafe > 500) {
                    console.log("Failure");
                    return this.faces;
                }
            } while (arc != start_arc)

            //face.reverse();
            this.faces.push(face);
        }
        return this.faces;
    }

    generate_components(one_orient=true) {
        let left_arcs = new Set(this.arcs);

        this.components = [];
        while (left_arcs.size > 0) {
            let start_arc = Array.from(left_arcs).pop();

            let component = this.component(start_arc);
            for (let arc of component) {
                left_arcs.delete(arc);
                if (one_orient) {
                    // Delete the other arc edge-opposite this one
                    left_arcs.delete(this.edges[arc.edge][(arc.edgepos+1)%2]);
                }
            }

            this.components.push(component);
        }
        return this.components;
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

    out_vert_i(arc) {
        return arc.vert;
    }

    in_vert_i(arc) {
        return this.edges[arc.edge][(arc.edgepos+1)%2].vert;
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

        //this.lab = new Lalolab();

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

        let K = this.K;
        let DeltaK = sub(target_K, K);

        let _failsafe = 0;
        while (this.loss(target_K) > thresh){
            let H = this.hessian();
            let deltau = least_squares(H, DeltaK);

            this.u = sub(this.u, mul(dt, deltau));

            this.update();

            K = this.K;
            DeltaK = sub(target_K, K);

            //console.log(math.max(math.abs(DeltaK)));

            _failsafe += 1;
            if (_failsafe > 1000) {
                console.log("Took too long to flatten; abort!");
            }
        }

    }

    newton_step(target_K=null, dt=0.05) {
        let DeltaK = sub(target_K, this.K);

        let H = this.hessian();
        let deltau = least_squares(H, DeltaK);

        this.u = sub(this.u, mul(dt, deltau));

        this.update();
    }


    loss(target_K) {
        return max(abs(sub(target_K, this.K)));
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
        //console.log(det(H));
        return H;
    }
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

        this.adj_map = {};
        for (let edge of edges) {
            let [a, b] = edge;
            if (a in this.adj_map) {
                this.adj_map[a].push(b);
            } else {
                this.adj_map[a] = [b];
            }

            if (b in this.adj_map) {
                this.adj_map[b].push(a);
            } else {
                this.adj_map[b] = [a];
            }
        }

        this.delta = 2;
        this.gamma = 5;

        this.dbar = 3*this.delta;

        this.a_exp = 1;
        this.er_exp = 2;

        this.calculate_surrounding_edges();
    }

    distance(u, v) {
        return norm(sub(u, v));
    }

    force_avert(u, v) {
        //console.log(this.distance(u,v));
        return mul(Math.pow(this.distance(u, v)/this.delta, this.a_exp),
                   sub(v, u));
    }

    force_rvert(u, v) {
        let d = this.distance(u, v);
        return mul(Math.pow(this.delta/d, this.er_exp),
                   sub(u, v));
    }

    compute_ve(v, a, b) {
        let m = (a[1] - b[1])/(a[0] - b[0]);
        let n = -1 / m;
        let c = a[1] - m*a[0];
        let d = v[1] - n*v[0];

        let x = (d - c) / (m - n);
        return [x, m*x + c];
    }

    ve_on_edge(ve, a, b) {
        return (((ve[0] <= a[0] && ve[0] >= b[0]) ||
                 (ve[0] <= b[0] && ve[0] >= a[0])) &&
                ((ve[1] <= a[1] && ve[1] >= b[1]) ||
                 (ve[1] <= b[1] && ve[1] >= a[1])));
    }

    force_redge(u, a, b, ve) {
        let d = this.distance(u, ve);
        if (d >= this.gamma) {
            // node and "virtual edge" too far
            return [0, 0];
        }

        return mul(-Math.pow(this.gamma - d, this.er_exp)/d,
                   sub(ve, u));
    }

    surrounding_edges(ui) {
        // calculate the surrounding edges S_ui
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

    calculate_surrounding_edges() {
        this.surr_edges = [];
        for (let i = 0; i < this.verts.length; i++) {
            this.surr_edges[i] = this.surrounding_edges(i);
        }
    }

    move (ui, F_ux, F_uy, M_u) {
        let i;
        if (F_ux >= 0) {
            if (F_uy >= 0) {
                if (F_ux >= F_uy) {
                    i = 0;
                } else {
                    i = 1;
                }
            } else {
                if (F_ux >= -F_uy) {
                    i = 7;
                } else {
                    i = 6;
                }
            }
        } else {
            if (F_uy >= 0) {
                if (-F_ux >= F_uy) {
                    i = 3;
                } else {
                    i = 2;
                }
            } else {
                if (-F_ux >= -F_uy) {
                    i = 4;
                } else {
                    i = 5;
                }
            }
        }

        let F_u = [F_ux, F_uy];

        let f_u = norm(F_u);
        let du;
        if (f_u <= M_u[i]) {
            du = F_u;
        } else {
            du = mul(M_u[i]/f_u, F_u);
        }

        //if (ui == 9) {
        //    console.log(ui, i, du, F_u, M_u[i]);
        //}

        //console.log(this.verts[ui], du, F_u, M_u);
        this.verts[ui][0] += du[0];
        this.verts[ui][1] += du[1];
        //console.log(F_u, M_u[i], du, this.verts[ui]);
    }

    update() {
        //console.log(this.adj_map);
        let F_x = zeros(this.verts.length);
        let F_y = zeros(this.verts.length);
        let M = [];
        for (let i = 0; i < this.verts.length; i++) {
            M.push([this.dbar, this.dbar, this.dbar, this.dbar,
                    this.dbar, this.dbar, this.dbar, this.dbar]);
        }

        let barycenter = mul(1/this.verts.length, sum(this.verts, 2));

        for (let ui = 0; ui < this.verts.length; ui++) {
            // Calculate gravity force
            let db = sub(barycenter, this.verts[ui]);
            let n_db = norm(db);
            F_x[ui] += db[0]/n_db;
            F_y[ui] += db[1]/n_db;

            // Calculate total node-node repulsive force
            for (let vi = 0; vi < this.verts.length; vi++) {
                if (ui != vi) {
                    if (this.distance(ui, vi) >= 3*this.delta) {
                        continue;
                    }

                    if (this.adj_map[ui].length == 2) {
                        if (this.adj_map[ui].includes(vi)) {
                            continue;
                        }
                    }

                    let F = this.force_rvert(this.verts[ui], this.verts[vi]);
                    //console.log("Fnnr", F);
                    //console.log(ui, vi, this.verts[ui], this.verts[vi], "Fnnr", F);
                    if (!isNaN(F[0])) {
                        F_x[ui] += F[0];
                        F_y[ui] += F[1];
                    }
                }
            }

            // calculate edge attractive force
            for (let vi of this.adj_map[ui]) {
                let F = this.force_avert(this.verts[ui], this.verts[vi]);

                F_x[ui] += F[0];
                F_y[ui] += F[1];
            }

            // calculate node-edge repulsive force
            for (let edge of this.surr_edges[ui]) {
                let [ai, bi] = edge;
                if (ui == ai || ui == bi) {
                    continue;
                }
                let ve = this.compute_ve(
                    this.verts[ui], this.verts[ai], this.verts[bi]);

                if (this.ve_on_edge(ve, this.verts[ai], this.verts[bi])) {
                    let F = this.force_redge(
                        this.verts[ui], this.verts[ai], this.verts[bi], ve);
                    if (!isNaN(F[0])) {
                        F_x[ui] += F[0];
                        F_y[ui] += F[1];
                    }
                }
            }

            let M_u = M[ui];
            //console.log("Surr:", this.surr_edges);

            for (let edge of this.surr_edges[ui]) {
                let [ai, bi] = edge;
                if (ui == ai || ui == bi) {
                    continue;
                }
                let ve = this.compute_ve(
                    this.verts[ui], this.verts[ai], this.verts[bi]);

                let cv;

                if (ui == 0 && ai == 5 && bi == 2) {
                    //console.log(this.verts[ai], this.verts[bi], ve);
                }
                //console.log("v-e", ui, ai, bi);
                if (this.ve_on_edge(ve, this.verts[ai], this.verts[bi])) {
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

                    let max_r = norm(cv)/2.1;
                    //console.log("???", cv);

                    //console.log(M_u, max_r, Math.cos(Math.atan2(cv[1], cv[0])));
                    let ell = (i+4)%8;
                    for (let j = 0; j < M_u.length; j++) {
                        if ((i-j+8)%8 == 0) {
                            M_u[j] = min(M_u[j], max_r);
                        } else if ((i-j+8)%8 == 1 || (i-j+8)%8 == 2) {
                            M_u[j] = min(M_u[j], max_r /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j+1)*Math.PI/4));
                        } else if ((i-j+8)%8 == 6 || (i-j+8)%8 == 7) {
                            M_u[j] = min(M_u[j], max_r /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    for (let j = 0; j < M_u.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[ai][j] = min(M[ai][j], max_r);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[ai][j] = min(M[ai][j], max_r /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[ai][j] = min(M[ai][j], max_r /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    for (let j = 0; j < M_u.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[bi][j] = min(M[bi][j], max_r);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[bi][j] = min(M[bi][j], max_r /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[bi][j] = min(M[bi][j], max_r /
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

                    let max_r = norm(cv)/2.1;
                    //console.log("???", cv);

                    //console.log(M_u, max_r, Math.cos(Math.atan2(cv[1], cv[0])));
                    let ell = (i+4)%8;
                    for (let j = 0; j < M_u.length; j++) {
                        if ((i-j+8)%8 == 0) {
                            M_u[j] = min(M_u[j], max_r);
                        } else if ((i-j+8)%8 == 1 || (i-j+8)%8 == 2) {
                            M_u[j] = min(M_u[j], max_r /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j+1)*Math.PI/4));
                        } else if ((i-j+8)%8 == 6 || (i-j+8)%8 == 7) {
                            M_u[j] = min(M_u[j], max_r /
                                         Math.cos(Math.atan2(cv[1], cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    let m = cv[1]/cv[0]; // Slope of cv
                    let n = -1 / m; // Slope of l

                    for (let j = 0; j < M_u.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[ai][j] = min(M[ai][j], max_r);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[ai][j] = min(M[ai][j], max_r /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[ai][j] = min(M[ai][j], max_r /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j)*Math.PI/4));
                        }
                    }

                    for (let j = 0; j < M_u.length; j++) {
                        if ((ell-j+8)%8 == 0) {
                            M[bi][j] = min(M[bi][j], max_r);
                        } else if ((ell-j+8)%8 == 1 || (ell-j+8)%8 == 2) {
                            M[bi][j] = min(M[bi][j], max_r /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j+1)*Math.PI/4));
                        } else if ((ell-j+8)%8 == 6 || (ell-j+8)%8 == 7) {
                            M[bi][j] = min(M[bi][j], max_r /
                                           Math.cos(Math.atan2(-cv[1], -cv[0]) - (j)*Math.PI/4));
                        }
                    }
                }

            }
            //if (ui == 0) console.log("M_u0", M_u, F_x[ui], F_y[ui]);
        }

        //console.log("Fx", F_x);
        for (let ui in this.verts) {

            this.move(ui, F_x[ui], F_y[ui], M[ui]);
        }
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
    while (to_orient.length > 0 && adj_queue.size > 0) {
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

    //let faces = orient_faces(mesh.faces);
    let faces = mesh.faces;
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

    //console.log("toEmbed:", to_embed);
    //console.log(mesh.faces);

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
            k = F[(F.indexOf(i)+2)%3];
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

    //console.log("toEmbed:", to_embed);
    return [x, faces, phi];
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

    let nVerts = self.randomDiagram(10, 1, 50, 0, Math.random()*(2**32), vertPtr);
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
    randomDiagram: function(n_verts, n_comps, max_att, type) {
        return randomDiagram(n_verts, n_comps, max_att, type);
    },

    setRandomLinkDiagram: function(n_verts, n_comps, max_att, type) {
        let sigma = randomDiagram(n_verts, n_comps, max_att, type);
        workerFunctions.setLinkDiagram(sigma, 8);
    },

    setLinkDiagram: function(sigma, cross_bend) {
        self.shadow = new LinkShadow(sigma);

        self.trign = self.shadow.triangulate();

        let ofaces = orient_faces(self.trign[3]);
        let testMesh = new TriangleMesh(
            self.trign[0], self.trign[1], self.trign[2], ofaces
        );

        let cpmetric = DiscreteRiemannMetric.from_triangle_mesh(testMesh);

        self.tgt_K = zeros(testMesh.verts.length, 1);

        // Set boundary crossing target curvature
        let bdyCross = testMesh.bdyverts.filter(vi => self.trign[4].verts.includes(vi));
        let bdyEdge = testMesh.bdyverts.filter(vi => !self.trign[4].verts.includes(vi));

        let fac = 8; // Inverse of how concave crossing vertices should be imbedded
        for (let bci of bdyCross) {
            self.tgt_K[bci] = -Math.PI/fac;
        }
        for (let bei of bdyEdge) {
            self.tgt_K[bei] = (2*Math.PI + bdyCross.length*Math.PI/fac)/bdyEdge.length;
        }

        self.flat_poly = cpmetric;

        workerFunctions.embedDiagram();
    },

    embedDiagram: function() {
        let tstart = Date.now();

        let thresh = 5e-10;
        let embedding = embed_faces(self.flat_poly);

        /*postMessage({
            function: "setEmbedding",
            arguments: [
                self.flat_poly,
                embedding,
                self.shadow,
                self.trign[4]]
        });*/

        //let curDate;
        //do { curDate = Date.now(); }
        //while( curDate-tstart < 100);

        while(self.flat_poly.loss(self.tgt_K) > thresh) {
            //let procStart = Date.now();

            self.flat_poly.newton_step(self.tgt_K, 1);

            embedding = embed_faces(self.flat_poly);

            /*postMessage({
                function: "updateEmbedding",
                arguments: [
                    self.flat_poly,
                    embedding,
                    self.shadow,
                    self.trign[4]]
            });*/

            //do { curDate = Date.now(); }
            //while( curDate-procStart < 1000);
        }

        console.log("Triangulation flattened in: " + (Date.now() - tstart) + " milliseconds");

        let pts = embedding[0];
        //console.log(pts);
        let min_x = min(get(pts, range(), 0));
        let min_y = min(get(pts, range(), 1));

        let max_x = max(get(pts, range(), 0));
        let max_y = max(get(pts, range(), 1));

        let wid = max_x-min_x;
        let hgt = max_y-min_y;

        let sqw = Math.min(wid, hgt);

        let mind = Infinity;
        for (let i = 0; i < embedding[0].m; i++) {
            for (let j = 0; i < embedding[0].m; i++) {
                if (i == j) { continue; }
                mind = Math.min(norm(sub(get(embedding[0], i, range()),
                                         get(embedding[0], j, range()))));
            }
        }

        embedding[0] = mul(10/mind , embedding[0]);
        console.log("minD: ", mind);

        pts = embedding[0];
        //console.log(pts);
        min_x = min(get(pts, range(), 0));
        min_y = min(get(pts, range(), 1));

        max_x = max(get(pts, range(), 0));
        max_y = max(get(pts, range(), 1));

        wid = max_x-min_x;
        hgt = max_y-min_y;

        sqw = Math.min(wid, hgt);

        mind = Infinity;
        for (let i = 0; i < embedding[0].m; i++) {
            for (let j = 0; i < embedding[0].m; i++) {
                if (i == j) { continue; }
                mind = Math.min(norm(sub(get(embedding[0], i, range()),
                                         get(embedding[0], j, range()))));
            }
        }

        console.log("minD: ", mind);

        // Create an embedded graph without scaffolding
        let l_verts = [];
        let l_edges = [];

        // Verts is a list of all triangulation vertices which are "real" graph
        // verts, as opposed to scaffolding (i.e. faces)
        let verts = self.trign[4].comps.reduce((res, a) => res.concat(a), []);
        verts = Array.from(new Set(verts));
        let vert_map = [];
        for (let i = 0; i < verts.length; i++) {
            // Create a backref in vert_map, for edges
            if (i == 63) { console.log("63 is", verts[i]); }
            vert_map[verts[i]] = i;

            // Push the point onto l_verts
            let x = embedding[0].val[verts[i]*embedding[0].n + 0];
            let y = embedding[0].val[verts[i]*embedding[0].n + 1];

            l_verts.push([x, y]);
            //console.log(l_verts);
        }

        console.log(vert_map);

        // Push the edges; edges are determined by components
        for (let comp of self.trign[4].comps) {
            for (let pi = 0; pi < comp.length; pi++) {
                l_edges.push([vert_map[comp[pi]], vert_map[comp[(pi+1)%comp.length]]]);
            }
        }

        self.force_shadow = new ForceLinkDiagram(
            l_verts, l_edges,
            self.trign[4].regions.map(r => r.map(vi => vert_map[vi])));

        //while(true) {
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

        console.log("Surr 56");
        console.log(self.force_shadow.surr_edges[56]);

        postMessage({
            function: "setLinkDiagram",
            arguments: [self.force_shadow]
        });
    }
}

onmessage = function(e) {
    workerFunctions[e.data.function](...e.data.arguments);
}
