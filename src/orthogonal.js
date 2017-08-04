//let lab = new Lalolab();

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

    set_nodes(start, end, anim_ms=0) {
        if (anim_ms == 0) {
            this.svg.attr({'x1': start.x(), 'y1': start.y(), 'x2': end.x(), 'y2': end.y()});
        } else {
            this.svg.animate({'x1': start.x(), 'y1': start.y(), 'x2': end.x(), 'y2': end.y()},
                             anim_ms);
        }
    }
}

class MeshNode {
    constructor(parent, x, y) {
        this.parent = parent;
        this.svg = this.parent.nodeG.circle(x, y, .3);
        this._x = x;
        this._y = y;
        this.svg.addClass('plnode');

        this.svg.node.addEventListener('click', this.onClick.bind(this));
    }

    set_obj(i, obj) {
        this.i = i;
        this.obj = obj;
    }

    move(x, y, new_r=undefined, anim_ms=0) {
        this._x = x;
        this._y = y;
        if (anim_ms == 0) {
            // Do not animate
            this.svg.attr({'cx': x, 'cy': y, 'r': new_r});
        } else {
            // Animate the motion
            this.svg.animate({'cx': x, 'cy': y, 'r': new_r},
                             anim_ms);
        }
    }

    set_r(r, anim_ms=0) {
        if (anim_ms == 0) {
            this.svg.attr({'r': r});
        } else {
            this.svg.animate({'r': r},
                             anim_ms);
        }
    }

    cur_x() {
        return this.svg.attr('cx');
    }
    x() {
        return this._x;
    }

    cur_y() {
        return this.svg.attr('cy');
    }
    y() {
        return this._y;
    }

    onClick(e) {
        //console.log(this);
        //console.log(this.obj);
        if (this.obj.length == 2) {
            this.parent.delete_face(this.i);
        }
    }

}

class CompEdgeNode {
    constructor(parent, x, y, r) {
        this.parent = parent;
        this.svg = this.parent.knotG.circle(x, y, r);
        this._x = x;
        this._y = y;
        this._r = r;

        this.svg.node.addEventListener('click', this.onClick.bind(this));

        this.dragging = false;
    }

    set_obj(i, obj) {
        this.i = i;
        this.obj = obj;
    }

    move(x, y, r, anim_ms=0) {
        this._x = x;
        this._y = y;
        this._r = r;
        if (anim_ms == 0) {
            // Do not animate
            this.svg.attr({'cx': x, 'cy': y, 'r': r});
        } else {
            // Animate the motion
            this.svg.animate({'cx': x, 'cy': y, 'r': r},
                             anim_ms);
        }
    }

    set_r(r, anim_ms=0) {
        if (anim_ms == 0) {
            this.svg.attr({'r': r});
        } else {
            this.svg.animate({'r': r},
                             anim_ms);
        }
    }

    cur_x() {
        return this.svg.attr('cx');
    }
    x() {
        return this._x;
    }

    cur_y() {
        return this.svg.attr('cy');
    }
    y() {
        return this._y;
    }

    onClick(e) {
        //console.log(this);
        //console.log(this.obj);
        if (this.obj.length == 2) {
            this.parent.delete_face(this.i);
        }
    }
}

class MeshDraw {
    constructor(div) {
        this.nodes = [];
        this.edges = {};
        this.comp_edgenodes = [];
        this.comps = [];

        this.anim_ms = 1000;
        /*this.pan = svgPanZoom(div, {
            minZoom: 0.1,
            maxZoom: 50,
            contain: true,
            controlIconsEnabled: true,
            zoomScaleSensitivity: 1
        });

        this.draw = Snap(document.querySelector(div).children[0]);*/
        this.draw = Snap(div);
        this.edgeG = this.draw.g();
        this.nodeG = this.draw.g();
        this.knotG = this.draw.g();
        this.knotG.addClass("Set1");
    }

    clear() {
        this.nodes = [];
        this.edges = {};
        this.comps = [];

        this.edgeG.clear();
        this.nodeG.clear();
        this.knotG.clear();
    }

    set_link_diagram(ld) {
        let pts = array2mat(ld.verts);
        //console.log(pts);
        let min_x = min(get(pts, range(), 0));
        let min_y = min(get(pts, range(), 1));

        let max_x = max(get(pts, range(), 0));
        let max_y = max(get(pts, range(), 1));

        let wid = max_x-min_x;
        let hgt = max_y-min_y;

        let dx = wid*0.05;
        let dy = hgt*0.05;

        //console.log(min_x, max_x);

        //this.nodeG.clear();
        //this.edgeG.clear();

        this.draw.attr({viewBox: [(min_x-dx),
                                  (min_y-dy),
                                  (wid+2*dx),
                                  (hgt+2*dy)].join(",")});

        let i = 0;
        //console.log(ld.verts, "!!");
        for (let vert of ld.verts) {
            //this.nodeG.circle(vert[0], vert[1], .25);
            if (this.nodes[i] === undefined) {
                let t = this.nodeG.text(vert[0], vert[1], i.toString());
                t.attr({"style": "font-size: .5px;"});
                this.nodes[i] = t;
            } else {
                this.nodes[i].attr({'x': vert[0], 'y': vert[1]});
            }
            i++;
        }

        for (let edge of ld.edges) {
            let [a, b] = [ld.verts[edge[0]], ld.verts[edge[1]]];
            if (this.edges[edge] === undefined) {
                this.edges[edge] = this.edgeG.line(a[0], a[1], b[0], b[1]).addClass("edge");
            } else {
                this.edges[edge].attr({'x1': a[0], 'y1': a[1], 'x2': b[0], 'y2': b[1]});
            }
        }
    }

    set_embedding(g /*metric*/, embedding, map4v /*original map*/, conv) {
        let points = embedding[0];
        let faces = embedding[1];
        let phi = embedding[2];
        /* set the embedding */
        this.g = g;
        this.map4v = map4v;

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

        //console.log(points);
        for (let i = 0; i < points.m; i++) {
            let node = this.add_node(i, points.val[i*points.n], points.val[i*points.n+1]);
            node.set_r(this.g.gamma[this.nodes.length-1]);
        }

        for (let comp of conv.comps) {
            for (let pi = 0; pi < comp.length; pi++) {
                let edge = [comp[pi], comp[(pi+1)%comp.length]];
                //console.log(edge);
                if (edge[0] in this.nodes && edge[1] in this.nodes) {
                    let mesh_edge = this.add_edge(edge[0], edge[1]);
                    mesh_edge.svg.addClass('edge');
                }
            }
        }

        //console.log(this.map4v.faces);
        for (let idx in conv.faces) {
            //console.log(fi);
            if(conv.faces[idx].length > 0) {
                let mesh_face = this.add_face(conv.faces[idx][0], parseInt(idx));
            }
        }

        let ci = 0;
        for (let component of conv.comps) {
            this.comps[ci] = this.add_component(component, map4v, points, conv);
            this.comps[ci].addClass('q'+ci+"-9");
            ci += 1;
        }
    }

    update_embedding(g /*metric*/, embedding, map4v /*original map*/, conv, anim_ms=0) {
        let points = embedding[0];
        let faces = embedding[1];
        let phi = embedding[2];
        /* set the embedding */
        this.g = g;
        this.map4v = map4v;

        // alter the viewbox to fit
        let min_x = min(get(points, range(), 0));
        let min_y = min(get(points, range(), 1));
        let max_x = max(get(points, range(), 0));
        let max_y = max(get(points, range(), 1));

        let wid = max_x-min_x;
        let hgt = max_y-min_y;

        let dx = wid*0.05;
        let dy = hgt*0.05;

        Snap.animate(this.draw.attr("viewBox").vb.split(" "),
                     [min_x-dx, min_y-dy, wid+2*dx, hgt+2*dy],
                     (v) => { this.draw.attr("viewBox", v.join(" ")); },
                     this.anim_ms);

        //console.log(points);
        for (let i = 0; i < points.m; i++) {
            let node = this.nodes[i];
            node.move(points.val[i*points.n], points.val[i*points.n+1],
                      this.g.gamma[i], this.anim_ms);
        }

        for (let comp of conv.comps) {
            for (let pi = 0; pi < comp.length; pi++) {
                let edge = [comp[pi], comp[(pi+1)%comp.length]];
                //console.log(edge);
                if (edge[0] in this.nodes && edge[1] in this.nodes) {
                    this.update_edge(edge[0], edge[1]);
                }
            }
        }

        let ci = 0;
        for (let component of conv.comps) {
            this.update_component(ci, component, map4v, points, conv);
            ci += 1;
        }
    }

    delete_face(i) {
        /* Deleting a face (valid curve-preserving faces to delete are a bigon *
         * or a monogon) is something of an ordeal: We want to delete only the
         * local triangles which the face affects, modifying the underlying
         * triangulation only locally, so that we can keep the old discrete
         * riemannian metric, hopefully so that the resulting embedding looks
         * very similar. Is there a proof that this will happen? */
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

    update_edge(i, j) {
        this.edges[[i,j]].set_nodes(this.nodes[i], this.nodes[j], this.anim_ms);
    }


    add_face(i, fi) {
        let face_node = this.nodes[i];
        face_node.svg.addClass("face");
        face_node.set_obj(i, this.map4v.faces[fi]);
    }

    component_path_gen(component, map4v, points) {
        let path = [];
        let pts = points;

        let root_arc = component[0];

        let last_ep = undefined;

        for (let i = 0; i < component.length; i++) {
            let j = (i+1) % component.length;

            let p_i = component[i];
            let p_j = component[j];
            let p_a = [pts.val[p_i*pts.n], pts.val[p_i*pts.n+1]];
            let p_b = [pts.val[p_j*pts.n], pts.val[p_j*pts.n+1]];

            let dp = sub(p_b, p_a);
            let whole_l = norm(dp);
            let mid_dp = mul(this.g.gamma[p_i]/whole_l, dp);
            let p_mid = add(p_a, mid_dp);

            path.push(p_mid); // midpoint is anchor
            path.push(p_b);   // p_b is quadratic bezier control
        }

        // Close the path
        path.push(path[0]);

        //console.log(path);

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

        return [path, pathStr];
    }

    quadratic_segments(path) {
        let slices = [];
        for (let i = 0; i < path.length-4; i += 2) {
            slices.push(path.slice(i, i+3));
        }
        slices.unshift(path.slice(path.length-3, path.length));

        return slices;
    }

    quadratic_segments_to_cubic(segs) {
        let csegs = [];
        for (let seg of segs) {
            let [q0, q1, q2] = seg;
            csegs.push([q0,
                        add(q0,mul(2/3 ,sub(q1,q0))),
                        add(q2,mul(2/3, sub(q1,q2))),
                        q2]);
        }
        return csegs;
    }

    add_component(component, map4v, points, conv) {
        // A path of the form anchor, control, anchor...
        let [path, pathStr] = this.component_path_gen(component, map4v, points);

        let segs = this.quadratic_segments_to_cubic(this.quadratic_segments(path));

        for (let si = 0; si < segs.length; si++) {
            let tri_i = component[si];
            if (!conv.edges.some(e => e.includes(tri_i))) {
                continue;
            }

            let seg = segs[si];

            let p = Snap.path.findDotsAtSegment(
                seg[0][0], seg[0][1],
                seg[1][0], seg[1][1],
                seg[2][0], seg[2][1],
                seg[3][0], seg[3][1], .5);
            this.comp_edgenodes[tri_i] = new CompEdgeNode(this, p.x, p.y, this.g.gamma[tri_i]/2);
        }


        //console.log(pathStr);
        let knot = this.knotG.path(pathStr+"Z");
        knot.addClass('knot');

        return knot;
    }

    update_component(ci, component, map4v, points, conv) {
        // A path of the form anchor, control, anchor...
        let [path, pathStr] = this.component_path_gen(component, map4v, points);

        let segs = this.quadratic_segments_to_cubic(this.quadratic_segments(path));

        for (let si = 0; si < segs.length; si++) {
            let tri_i = component[si];
            if (!conv.edges.some(e => e.includes(tri_i))) {
                continue;
            }

            let seg = segs[si];

            let p = Snap.path.findDotsAtSegment(
                seg[0][0], seg[0][1],
                seg[1][0], seg[1][1],
                seg[2][0], seg[2][1],
                seg[3][0], seg[3][1], .5);

            this.comp_edgenodes[tri_i].move(p.x, p.y, this.g.gamma[tri_i]/2, this.anim_ms);
        }


        //console.log(pathStr);
        this.comps[ci].stop();
        this.comps[ci].animate({"d": pathStr+"Z"}, this.anim_ms);
    }
}


let meshDraw = new MeshDraw("#knot-draw");

let cpWorker = new Worker("js/orth_worker.js");

function drawMapAsync(sigma, cross_bend=8) {
    cpWorker.postMessage({
        function: "setLinkDiagram",
        arguments: [sigma, cross_bend]
    });
}

var cpWorkerFunctions = {
    setEmbedding: function(flat_poly, embedding, m4v, conv) {
        meshDraw.clear();
        meshDraw.set_embedding(flat_poly, embedding, m4v, conv);
    },

    updateEmbedding: function(flat_poly, embedding, m4v, conv) {
        meshDraw.update_embedding(flat_poly, embedding, m4v, conv);
    },

    setLinkDiagram: function(link_diagram) {
        meshDraw.set_link_diagram(link_diagram);
    }
}

cpWorker.onmessage = function(ev) {
    //console.log("hello", ev);
    cpWorkerFunctions[ev.data.function].apply(this, ev.data.arguments);
}

// Trefoil
//let sigma = [[0, 6, 11, 5], [1, 8, 2, 7], [3, 9, 4, 10]];

// Twist
//let sigma = [[0, 1, 3, 2]];

// 2 Twist
//let sigma = [[0, 1, 7, 2], [3, 6, 4, 5]];

// Monogon in internal face
//let sigma = [[1, 8, 2, 7], [0, 14, 15, 13], [3, 9, 4, 10], [5, 12, 6, 11]];

// Small 2-link
// Fails orientation?
//let sigma = [[1, 4, 2, 3], [0, 8, 7, 15], [5, 14, 6, 13], [9, 12, 10, 11]];

// (2,5) torus
//let sigma = [[0, 10, 19, 9], [1, 12, 2, 11], [3, 13, 4, 14], [5, 16, 6, 15], [7, 17, 8, 18]];

// Random composite
//let sigma = [[0, 57, 79, 58], [1, 59, 2, 60], [3, 78, 4, 77], [5, 11, 6, 12], [7, 9, 8, 10], [13, 16, 14, 15], [17, 56, 18, 55], [19, 54, 20, 53], [21, 36, 22, 35], [23, 25, 24, 26], [27, 29, 28, 30], [31, 46, 32, 45], [33, 39, 34, 40], [37, 52, 38, 51], [41, 44, 42, 43], [47, 50, 48, 49], [61, 72, 62, 71], [63, 73, 64, 74], [65, 68, 66, 67], [69, 76, 70, 75]];

//let sigma = [[1, 36, 2, 35], [0, 38, 75, 37], [3, 18, 4, 17], [5, 28, 6, 27], [7, 21, 8, 22], [9, 31, 10, 32], [11, 34, 12, 33], [13, 24, 14, 23], [15, 25, 16, 26], [19, 30, 20, 29], [39, 73, 40, 74], [41, 72, 42, 71], [43, 66, 44, 65], [45, 51, 46, 52], [47, 61, 48, 62], [49, 60, 50, 59], [53, 63, 54, 64], [55, 69, 56, 70], [57, 68, 58, 67]];

//let sigma = [[1, 3, 2, 4], [0, 60, 59, 119], [5, 15, 6, 16], [7, 84, 8, 83], [9, 101, 10, 102], [11, 92, 12, 91], [13, 89, 14, 90], [17, 62, 18, 61], [19, 79, 20, 80], [21, 70, 22, 69], [23, 42, 24, 41], [25, 28, 26, 27], [29, 32, 30, 31], [33, 52, 34, 51], [35, 37, 36, 38], [39, 49, 40, 50], [43, 54, 44, 53], [45, 48, 46, 47], [55, 67, 56, 68], [57, 82, 58, 81], [63, 77, 64, 78], [65, 72, 66, 71], [73, 75, 74, 76], [85, 88, 86, 87], [93, 99, 94, 100], [95, 98, 96, 97], [103, 118, 104, 117], [105, 112, 106, 111], [107, 109, 108, 110], [113, 115, 114, 116]];

// Complicated 2-link
// TODO: Fails orthogonal -- getting negative flow
let sigma = [[1, 4, 2, 3], [0, 66, 295, 65], [5, 64, 6, 63], [7, 49, 8, 50], [9, 12, 10, 11], [13, 43, 14, 44], [15, 46, 16, 45], [17, 19, 18, 20], [21, 24, 22, 23], [25, 32, 26, 31], [27, 30, 28, 29], [33, 47, 34, 48], [35, 38, 36, 37], [39, 42, 40, 41], [51, 54, 52, 53], [55, 294, 56, 293], [57, 188, 58, 187], [59, 289, 60, 290], [61, 175, 62, 176], [67, 174, 68, 173], [69, 140, 70, 139], [71, 278, 72, 277], [73, 131, 74, 132], [75, 78, 76, 77], [79, 117, 80, 118], [81, 87, 82, 88], [83, 86, 84, 85], [89, 104, 90, 103], [91, 94, 92, 93], [95, 98, 96, 97], [99, 101, 100, 102], [105, 108, 106, 107], [109, 111, 110, 112], [113, 256, 114, 255], [115, 257, 116, 258], [119, 122, 120, 121], [123, 134, 124, 133], [125, 275, 126, 276], [127, 137, 128, 138], [129, 244, 130, 243], [135, 269, 136, 270], [141, 280, 142, 279], [143, 145, 144, 146], [147, 241, 148, 242], [149, 151, 150, 152], [153, 155, 154, 156], [157, 192, 158, 191], [159, 169, 160, 170], [161, 163, 162, 164], [165, 168, 166, 167], [171, 189, 172, 190], [177, 184, 178, 183], [179, 181, 180, 182], [185, 292, 186, 291], [193, 227, 194, 228], [195, 198, 196, 197], [199, 285, 200, 286], [201, 208, 202, 207], [203, 222, 204, 221], [205, 211, 206, 212], [209, 224, 210, 223], [213, 219, 214, 220], [215, 217, 216, 218], [225, 284, 226, 283], [229, 287, 230, 288], [231, 234, 232, 233], [235, 298, 236, 297], [237, 299, 238, 296], [239, 282, 240, 281], [245, 268, 246, 267], [247, 249, 248, 250], [251, 265, 252, 266], [253, 260, 254, 259], [261, 264, 262, 263], [271, 274, 272, 273]];

// More complicated
//let sigma = [[0, 41, 99, 42], [1, 47, 2, 48], [3, 34, 4, 33], [5, 96, 6, 95], [7, 13, 8, 14], [9, 68, 10, 67], [11, 69, 12, 70], [15, 62, 16, 61], [17, 59, 18, 60], [19, 85, 20, 86], [21, 27, 22, 28], [23, 82, 24, 81], [25, 83, 26, 84], [29, 88, 30, 87], [31, 93, 32, 94], [35, 46, 36, 45], [37, 43, 38, 44], [39, 98, 40, 97], [49, 76, 50, 75], [51, 73, 52, 74], [53, 72, 54, 71], [55, 65, 56, 66], [57, 64, 58, 63], [77, 92, 78, 91], [79, 89, 80, 90]];

//let sigma = [[1, 44, 2, 59], [0, 58, 43, 57], [3, 41, 4, 42], [5, 16, 6, 15], [7, 25, 8, 26], [9, 32, 10, 31], [11, 54, 12, 53], [13, 47, 14, 48], [17, 36, 18, 35], [19, 37, 20, 38], [21, 40, 22, 39], [23, 33, 24, 34], [27, 50, 28, 49], [29, 51, 30, 52], [45, 56, 46, 55]];

//let sigma = [[1, 3, 2, 4], [0, 48, 47, 59], [5, 11, 6, 12], [7, 9, 8, 10], [13, 43, 14, 44], [15, 25, 16, 26], [17, 24, 18, 23], [19, 22, 20, 21], [27, 42, 28, 41], [29, 31, 30, 32], [33, 35, 34, 36], [37, 40, 38, 39], [45, 58, 46, 57], [49, 51, 50, 52], [53, 55, 54, 56]];

// Even moreso
//let sigma = [[1, 287, 2, 288], [0, 289, 299, 290], [3, 126, 4, 125], [5, 16, 6, 15], [7, 13, 8, 14], [9, 255, 10, 256], [11, 254, 12, 253], [17, 127, 18, 128], [19, 282, 20, 281], [21, 132, 22, 131], [23, 277, 24, 278], [25, 235, 26, 236], [27, 222, 28, 221], [29, 223, 30, 224], [31, 38, 32, 37], [33, 228, 34, 227], [35, 225, 36, 226], [39, 230, 40, 229], [41, 156, 42, 155], [43, 157, 44, 158], [45, 164, 46, 163], [47, 165, 48, 166], [49, 56, 50, 55], [51, 198, 52, 197], [53, 195, 54, 196], [57, 243, 58, 244], [59, 242, 60, 241], [61, 231, 62, 232], [63, 234, 64, 233], [65, 248, 66, 247], [67, 249, 68, 250], [69, 276, 70, 275], [71, 133, 72, 134], [73, 83, 74, 84], [75, 82, 76, 81], [77, 283, 78, 284], [79, 286, 80, 285], [85, 295, 86, 296], [87, 294, 88, 293], [89, 291, 90, 292], [91, 298, 92, 297], [93, 136, 94, 135], [95, 109, 96, 110], [97, 108, 98, 107], [99, 122, 100, 121], [101, 259, 102, 260], [103, 262, 104, 261], [105, 119, 106, 120], [111, 273, 112, 274], [113, 272, 114, 271], [115, 265, 116, 266], [117, 264, 118, 263], [123, 137, 124, 138], [129, 280, 130, 279], [139, 257, 140, 258], [141, 268, 142, 267], [143, 269, 144, 270], [145, 252, 146, 251], [147, 237, 148, 238], [149, 216, 150, 215], [151, 217, 152, 218], [153, 220, 154, 219], [159, 210, 160, 209], [161, 211, 162, 212], [167, 190, 168, 189], [169, 187, 170, 188], [171, 205, 172, 206], [173, 200, 174, 199], [175, 201, 176, 202], [177, 204, 178, 203], [179, 186, 180, 185], [181, 191, 182, 192], [183, 194, 184, 193], [207, 214, 208, 213], [239, 245, 240, 246]];

drawMapAsync(sigma);

document.getElementById("map_submit").onclick = function(ev) {
    document.querySelector(ev.target.dataset.errbox).textContent = "";
    try {
        meshDraw.clear();
        drawMapAsync(JSON.parse(document.querySelector(ev.target.dataset.target).value));
    } catch(err) {
        document.querySelector(ev.target.dataset.errbox).textContent = "Error";
        console.log("Error:", err);
    }
};
