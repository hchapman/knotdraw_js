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
}

class MeshNode {
    constructor(parent, x, y) {
        this.parent = parent;
        this.svg = this.parent.nodeG.circle(x, y, .3);
        this.svg.addClass('plnode');

        this.svg.node.addEventListener('click', this.onClick.bind(this));
    }

    set_obj(i, obj) {
        this.i = i;
        this.obj = obj;
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
        console.log(this);
        console.log(this.obj);
        if (this.obj.length == 2) {
            this.parent.delete_face(this.i);
        }
    }

}

class MeshDraw {
    constructor(div) {
        this.nodes = [];
        this.edges = {};
        /*this.pan = svgPanZoom(div, {
            minZoom: 0.1,
            maxZoom: 50,
            fit: false,
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

        this.edgeG.clear();
        this.nodeG.clear();
        this.knotG.clear();
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
                console.log(edge);
                if (edge[0] in this.nodes && edge[1] in this.nodes) {
                    let mesh_edge = this.add_edge(edge[0], edge[1]);
                    mesh_edge.svg.addClass('edge');
                }
            }
        }

        //console.log(this.map4v.faces);
        for (let fi in this.map4v.faces) {
            //console.log(fi);
            //let mesh_face = this.add_face(parseInt(fi));
        }

        let ci = 0;
        for (let component of map4v.components) {
            //let knot = this.add_component(component, map4v, points);
            //knot.addClass('q'+ci+"-9");
            //ci += 1;
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

    add_face(i) {
        let face_node = this.nodes[this.map4v.nv+this.map4v.ne+i];
        face_node.svg.addClass("face");
        face_node.set_obj(i, this.map4v.faces[i]);
    }

    add_component(component, map4v, points) {
        // A path of the form anchor, control, anchor...
        let path = [];
        let pts = points;

        let root_arc = component[0];

        let last_ep = undefined;
        for (let arc of component) {
            let vi = arc.vert;
            let ei = map4v.nv+arc.edge;

            let vp = [pts.val[vi*pts.n], pts.val[vi*pts.n+1]];
            let ep = [pts.val[ei*pts.n], pts.val[ei*pts.n+1]];

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

        let vp = [pts.val[vi*pts.n], pts.val[vi*pts.n+1]];
        let ep = [pts.val[ei*pts.n], pts.val[vi*pts.n+1]];

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

        return knot;
    }
}

let meshDraw = new MeshDraw("#knot-draw");

let cpWorker = new Worker("js/cp_worker.js");

function drawMapAsync(sigma, cross_bend=8) {
    cpWorker.postMessage([sigma, cross_bend]);
}

cpWorker.onmessage = function(ev) {
    //console.log(ev.data.m4v);
    meshDraw.clear();
    meshDraw.set_embedding(ev.data.flat_poly, ev.data.embedding, ev.data.m4v, ev.data.conv);
}

// Trefoil
//let sigma = [[0, 6, 11, 5], [1, 8, 2, 7], [3, 9, 4, 10]];

// Twist
//let sigma = [[0, 1, 3, 2]];

// Monogon in internal face
let sigma = [[1, 8, 2, 7], [0, 14, 15, 13], [3, 9, 4, 10], [5, 12, 6, 11]];

// More complicated
//let sigma = [[0, 41, 99, 42], [1, 47, 2, 48], [3, 34, 4, 33], [5, 96, 6, 95], [7, 13, 8, 14], [9, 68, 10, 67], [11, 69, 12, 70], [15, 62, 16, 61], [17, 59, 18, 60], [19, 85, 20, 86], [21, 27, 22, 28], [23, 82, 24, 81], [25, 83, 26, 84], [29, 88, 30, 87], [31, 93, 32, 94], [35, 46, 36, 45], [37, 43, 38, 44], [39, 98, 40, 97], [49, 76, 50, 75], [51, 73, 52, 74], [53, 72, 54, 71], [55, 65, 56, 66], [57, 64, 58, 63], [77, 92, 78, 91], [79, 89, 80, 90]];

// Even moreso
//let sigma = [[1, 287, 2, 288], [0, 289, 299, 290], [3, 126, 4, 125], [5, 16, 6, 15], [7, 13, 8, 14], [9, 255, 10, 256], [11, 254, 12, 253], [17, 127, 18, 128], [19, 282, 20, 281], [21, 132, 22, 131], [23, 277, 24, 278], [25, 235, 26, 236], [27, 222, 28, 221], [29, 223, 30, 224], [31, 38, 32, 37], [33, 228, 34, 227], [35, 225, 36, 226], [39, 230, 40, 229], [41, 156, 42, 155], [43, 157, 44, 158], [45, 164, 46, 163], [47, 165, 48, 166], [49, 56, 50, 55], [51, 198, 52, 197], [53, 195, 54, 196], [57, 243, 58, 244], [59, 242, 60, 241], [61, 231, 62, 232], [63, 234, 64, 233], [65, 248, 66, 247], [67, 249, 68, 250], [69, 276, 70, 275], [71, 133, 72, 134], [73, 83, 74, 84], [75, 82, 76, 81], [77, 283, 78, 284], [79, 286, 80, 285], [85, 295, 86, 296], [87, 294, 88, 293], [89, 291, 90, 292], [91, 298, 92, 297], [93, 136, 94, 135], [95, 109, 96, 110], [97, 108, 98, 107], [99, 122, 100, 121], [101, 259, 102, 260], [103, 262, 104, 261], [105, 119, 106, 120], [111, 273, 112, 274], [113, 272, 114, 271], [115, 265, 116, 266], [117, 264, 118, 263], [123, 137, 124, 138], [129, 280, 130, 279], [139, 257, 140, 258], [141, 268, 142, 267], [143, 269, 144, 270], [145, 252, 146, 251], [147, 237, 148, 238], [149, 216, 150, 215], [151, 217, 152, 218], [153, 220, 154, 219], [159, 210, 160, 209], [161, 211, 162, 212], [167, 190, 168, 189], [169, 187, 170, 188], [171, 205, 172, 206], [173, 200, 174, 199], [175, 201, 176, 202], [177, 204, 178, 203], [179, 186, 180, 185], [181, 191, 182, 192], [183, 194, 184, 193], [207, 214, 208, 213], [239, 245, 240, 246]];

drawMapAsync(sigma);

document.getElementById("map_submit").onclick = function(ev) {
    document.querySelector(ev.target.dataset.errbox).textContent = "";
    try {
        drawMapAsync(JSON.parse(document.querySelector(ev.target.dataset.target).value));
    } catch(err) {
        document.querySelector(ev.target.dataset.errbox).textContent = "Error";
        console.log("Error:", err);
    }
};
