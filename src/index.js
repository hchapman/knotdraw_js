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
                //console.log(edge);
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
        for (let component of conv.comps) {
            let knot = this.add_component(component, map4v, points);
            knot.addClass('q'+ci+"-9");
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

        console.log(path);

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

// 2 Twist
//let sigma = [[0, 1, 7, 2], [3, 6, 4, 5]];

// Monogon in internal face
//let sigma = [[1, 8, 2, 7], [0, 14, 15, 13], [3, 9, 4, 10], [5, 12, 6, 11]];

// Small 2-link
//let sigma = [[1, 4, 2, 3], [0, 8, 7, 15], [5, 14, 6, 13], [9, 12, 10, 11]];

// Random composite
//let sigma = [[0, 57, 79, 58], [1, 59, 2, 60], [3, 78, 4, 77], [5, 11, 6, 12], [7, 9, 8, 10], [13, 16, 14, 15], [17, 56, 18, 55], [19, 54, 20, 53], [21, 36, 22, 35], [23, 25, 24, 26], [27, 29, 28, 30], [31, 46, 32, 45], [33, 39, 34, 40], [37, 52, 38, 51], [41, 44, 42, 43], [47, 50, 48, 49], [61, 72, 62, 71], [63, 73, 64, 74], [65, 68, 66, 67], [69, 76, 70, 75]];

// Complicated 2-link
//let sigma = [[1, 4, 2, 3], [0, 66, 295, 65], [5, 64, 6, 63], [7, 49, 8, 50], [9, 12, 10, 11], [13, 43, 14, 44], [15, 46, 16, 45], [17, 19, 18, 20], [21, 24, 22, 23], [25, 32, 26, 31], [27, 30, 28, 29], [33, 47, 34, 48], [35, 38, 36, 37], [39, 42, 40, 41], [51, 54, 52, 53], [55, 294, 56, 293], [57, 188, 58, 187], [59, 289, 60, 290], [61, 175, 62, 176], [67, 174, 68, 173], [69, 140, 70, 139], [71, 278, 72, 277], [73, 131, 74, 132], [75, 78, 76, 77], [79, 117, 80, 118], [81, 87, 82, 88], [83, 86, 84, 85], [89, 104, 90, 103], [91, 94, 92, 93], [95, 98, 96, 97], [99, 101, 100, 102], [105, 108, 106, 107], [109, 111, 110, 112], [113, 256, 114, 255], [115, 257, 116, 258], [119, 122, 120, 121], [123, 134, 124, 133], [125, 275, 126, 276], [127, 137, 128, 138], [129, 244, 130, 243], [135, 269, 136, 270], [141, 280, 142, 279], [143, 145, 144, 146], [147, 241, 148, 242], [149, 151, 150, 152], [153, 155, 154, 156], [157, 192, 158, 191], [159, 169, 160, 170], [161, 163, 162, 164], [165, 168, 166, 167], [171, 189, 172, 190], [177, 184, 178, 183], [179, 181, 180, 182], [185, 292, 186, 291], [193, 227, 194, 228], [195, 198, 196, 197], [199, 285, 200, 286], [201, 208, 202, 207], [203, 222, 204, 221], [205, 211, 206, 212], [209, 224, 210, 223], [213, 219, 214, 220], [215, 217, 216, 218], [225, 284, 226, 283], [229, 287, 230, 288], [231, 234, 232, 233], [235, 298, 236, 297], [237, 299, 238, 296], [239, 282, 240, 281], [245, 268, 246, 267], [247, 249, 248, 250], [251, 265, 252, 266], [253, 260, 254, 259], [261, 264, 262, 263], [271, 274, 272, 273]];

// More complicated
//let sigma = [[0, 41, 99, 42], [1, 47, 2, 48], [3, 34, 4, 33], [5, 96, 6, 95], [7, 13, 8, 14], [9, 68, 10, 67], [11, 69, 12, 70], [15, 62, 16, 61], [17, 59, 18, 60], [19, 85, 20, 86], [21, 27, 22, 28], [23, 82, 24, 81], [25, 83, 26, 84], [29, 88, 30, 87], [31, 93, 32, 94], [35, 46, 36, 45], [37, 43, 38, 44], [39, 98, 40, 97], [49, 76, 50, 75], [51, 73, 52, 74], [53, 72, 54, 71], [55, 65, 56, 66], [57, 64, 58, 63], [77, 92, 78, 91], [79, 89, 80, 90]];

// Even moreso
let sigma = [[1, 287, 2, 288], [0, 289, 299, 290], [3, 126, 4, 125], [5, 16, 6, 15], [7, 13, 8, 14], [9, 255, 10, 256], [11, 254, 12, 253], [17, 127, 18, 128], [19, 282, 20, 281], [21, 132, 22, 131], [23, 277, 24, 278], [25, 235, 26, 236], [27, 222, 28, 221], [29, 223, 30, 224], [31, 38, 32, 37], [33, 228, 34, 227], [35, 225, 36, 226], [39, 230, 40, 229], [41, 156, 42, 155], [43, 157, 44, 158], [45, 164, 46, 163], [47, 165, 48, 166], [49, 56, 50, 55], [51, 198, 52, 197], [53, 195, 54, 196], [57, 243, 58, 244], [59, 242, 60, 241], [61, 231, 62, 232], [63, 234, 64, 233], [65, 248, 66, 247], [67, 249, 68, 250], [69, 276, 70, 275], [71, 133, 72, 134], [73, 83, 74, 84], [75, 82, 76, 81], [77, 283, 78, 284], [79, 286, 80, 285], [85, 295, 86, 296], [87, 294, 88, 293], [89, 291, 90, 292], [91, 298, 92, 297], [93, 136, 94, 135], [95, 109, 96, 110], [97, 108, 98, 107], [99, 122, 100, 121], [101, 259, 102, 260], [103, 262, 104, 261], [105, 119, 106, 120], [111, 273, 112, 274], [113, 272, 114, 271], [115, 265, 116, 266], [117, 264, 118, 263], [123, 137, 124, 138], [129, 280, 130, 279], [139, 257, 140, 258], [141, 268, 142, 267], [143, 269, 144, 270], [145, 252, 146, 251], [147, 237, 148, 238], [149, 216, 150, 215], [151, 217, 152, 218], [153, 220, 154, 219], [159, 210, 160, 209], [161, 211, 162, 212], [167, 190, 168, 189], [169, 187, 170, 188], [171, 205, 172, 206], [173, 200, 174, 199], [175, 201, 176, 202], [177, 204, 178, 203], [179, 186, 180, 185], [181, 191, 182, 192], [183, 194, 184, 193], [207, 214, 208, 213], [239, 245, 240, 246]];

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
