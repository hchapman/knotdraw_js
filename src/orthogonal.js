
class MeshDraw {
    constructor(div) {
        this.nodes = [];
        this.edges = {};
        this.comp_edgenodes = [];
        this.comps = [];

        this.draw = document.querySelector(div);
        let self = this;
        paper.setup(this.draw);
        paper.project.view.on('resize', function(e) {
            self.resizeView();
            return;
        });

        this.edgeG = new paper.Layer();
        this.nodeG = new paper.Layer();
        this.faceG = new paper.Layer();
        this.knotG = new paper.Layer();
    }

    clear() {
        this.nodes = [];
        this.edges = {};
        this.comps = [];
        this.comp_edgenodes = [];

        this.edgeG.clear();
        this.nodeG.clear();
        this.faceG.clear();
        this.knotG.clear();
    }

    resizeView() {
        let viewSize = paper.project.view.getViewSize();

        if (viewSize.height/this.viewHeight < viewSize.width/this.viewWidth) {
            this.zoom = viewSize.height / this.viewHeight;
        } else {
            this.zoom = viewSize.width / this.viewWidth;
        }
        paper.project.view.setZoom(this.zoom);

        paper.project.view.setCenter([this.min_x+this.wid/2, this.min_y+this.hgt/2]);
    }

    set_link_diagram(ld) {
        let pts = array2mat(ld.verts);
        //console.log(pts);
        this.min_x = Array.from(ld.verts.values()).reduce((m, [x,y]) => Math.min(m, x), Infinity);
        this.min_y = Array.from(ld.verts.values()).reduce((m, [x,y]) => Math.min(m, y), Infinity);

        this.max_x = Array.from(ld.verts.values()).reduce((m, [x,y]) => Math.max(m, x), -Infinity);
        this.max_y = Array.from(ld.verts.values()).reduce((m, [x,y]) => Math.max(m, y), -Infinity);

        this.wid = this.max_x-this.min_x;
        this.hgt = this.max_y-this.min_y;

        this.dx = this.wid*0.05;
        this.dy = this.hgt*0.05;

        this.viewHeight = this.hgt + this.dy*2;
        this.viewWidth = this.wid + this.dx*2;

        //console.log(min_x, max_x);

        //this.nodeG.clear();
        //this.edgeG.clear();

        paper.project.clear();

        /*let bds = paper.project.view.bounds;
        if (bds.height < bds.width) {
            paper.project.view.scale(bds.height / (this.hgt+2*this.dy));
        } else {
            paper.project.view.scale(bds.width / (this.wid+2*this.dy));
            }*/

        this.resizeView();

        let c = new paper.Path.Rectangle(this.min_x, this.min_y, this.wid, this.hgt);

        c.strokeColor = 'blue';
        c.strokeWidth = 1;
        c.strokeScaling = false;

        let i = 0;
        //console.log(ld.verts, "!!");
        for (let [vi, vert] of ld.verts) {
            //break;
            if (vert === undefined) { continue; }
            //this.nodeG.circle(vert[0], vert[1], .25);
            if (this.nodes[i] === undefined) {
                let t = new paper.PointText({
                    point: [vert[0], vert[1]],
                    content: ""+vi,
                    fontSize: 12,
                });
                //console.log(vi, vert);
                t.scale(1/this.zoom, [vert[0], vert[1]]);
                //console.log(t.scaling, t.applyMatrix);
                //let t = this.nodeG.text(vert[0], vert[1], i.toString());
                //t.attr({"style": "font-size: .5px;"});
                //this.nodes[i] = t;
            } else {
                //this.nodes[i].attr({'x': vert[0], 'y': vert[1]});
            }
            i++;
        }

        for (let edge of ld.edges) {
            let [a, b] = [ld.verts.get(edge[0]), ld.verts.get(edge[1])];
            if (this.edges[edge] === undefined) {
                let edgePath = new paper.Path();
                edgePath.strokeColor = 'green';
                edgePath.strokeWidth = .25;
                edgePath.strokeScaling = false;
                edgePath.add(new paper.Point(a[0], a[1]));
                edgePath.add(new paper.Point(b[0], b[1]));
                //this.edges[edge] = this.edgeG.line(a[0], a[1], b[0], b[1]).addClass("edge");
            } else {
                //this.edges[edge].attr({'x1': a[0], 'y1': a[1], 'x2': b[0], 'y2': b[1]});
            }
        }
    }

    finalize_link_diagram(ld, faces, components) {
        let min_radii = [];
        let vEArray = Array.from(ld.verts.entries());
        for (let [ai, a] of ld.verts) {
            //console.log(ld.adjMap[ai]);
            if (ld.adjMap[ai].length != 2) {
                // ai is a crossing
                min_radii[ai] = vEArray.reduce((r, [pi, p]) => ai==pi ? r : Math.min(r, norm(sub(a, p))), Infinity)/2;
                console.assert(min_radii[ai] != 0, min_radii, ai);
            } else {
                // ai is an edge
                let [bi, ci] = ld.adjMap[ai];
                let [b, c] = [ld.verts.get(bi), ld.verts.get(ci)];
                let u = [b[0]-a[0], b[1]-a[1]];
                let v = [c[0]-a[0], c[1]-a[1]];

                let phi_u = Math.atan2(u[1], u[0]);
                let phi_v = Math.atan2(v[1], v[0]);

                // we only want to consider points which are in the
                // angle between u and v.
                // we check phi_u - phi_v.
                if ((phi_u - phi_v + (Math.PI*2)) % (Math.PI * 2) > Math.PI) {
                    let tmp = phi_u;
                    phi_v = phi_u;
                    phi_u = tmp; // (phi_u - phi_v) % Math.PI*2 <= Math.PI
                }

                let phi = (phi_u - phi_v + (Math.PI*2)) % (Math.PI * 2);
                console.assert(phi <= Math.PI);

                min_radii[ai] = Math.min(...vEArray.filter(
                    ([xi, x]) => {
                        if (xi == ai) { return false; }
                        if (xi == bi || xi == ci) { return true; }
                        let w = [x[0]-a[0], x[1]-a[1]];
                        let phi_w = Math.atan2(w[1], w[0]);
                        return (phi_u - phi_w + (Math.PI*2)) % (Math.PI*2) < phi;
                    }
                ).map(([xi, x]) => norm(sub(a, x))))/2;
                //console.log(min_radii[ai]);
            }

            //console.log(ai, a);
            let nodeCircle = new paper.Path.Circle(a[0], a[1], min_radii[ai]);
            nodeCircle.strokeColor = 'red';
            nodeCircle.strokeWidth = .1;
            nodeCircle.strokeScaling = false;
            //let t = this.nodeG.circle(a[0], a[1], min_radii[ai]);
            //t.addClass("plnode");
        }
        paper.project.view.update();

        let vertPaths = [];
        let edgePaths = [];
        let crssPaths = [];
        let p;

        const N = Math.max(...ld.verts.keys())+1;
        //console.log(N);
        function tup(a, b, c) {
            if (b === undefined) {
                return a;
            } else if (c === undefined) {
                return a*N+b;
            } else {
                return a*N*N + b*N + c;
            }
        }

        for (let ci = 0; ci < components.length; ci++) {
            let component = components[ci];
            let cmp = [component[component.length-1]].concat(component).concat([component[0]]);
            let laststop, firststart;

            for (let i = 0; i < cmp.length-2; i++) {
                let bend = cmp.slice(i, i+3);

                let dstart = sub(ld.verts.get(bend[0]), ld.verts.get(bend[1]));
                let start = add(ld.verts.get(bend[1]),
                                mul(min_radii[bend[1]]/norm(dstart), dstart));
                if (firststart == undefined) { firststart = start; }
                let ctrl = ld.verts.get(bend[1]);
                let dstop = sub(ld.verts.get(bend[2]), ld.verts.get(bend[1]));
                let stop = add(ld.verts.get(bend[1]),
                               mul(min_radii[bend[1]]/norm(dstop), dstop));
                if (laststop !== undefined) {
                    //pStr = "M"+laststop+"L"+start+"Q"+ctrl+" "+stop;
                    p = new paper.Path("M"+laststop+"L"+start+"Q"+ctrl+" "+stop);
                } else {
                    //pStr = "M"+start+"C"+ctrl+ctrl+" "+stop;
                    p = new paper.Path("M"+start+"Q"+ctrl+" "+stop);
                }
                p.strokeColor = 'black';
                p.strokeWidth = 1;
                p.strokeScaling = false;

                if (vertPaths[bend[1]] === undefined) {
                    vertPaths[bend[1]] = {type: "edge",
                                          idx: bend};

                    edgePaths[tup(...bend)] = p.clone({insert: false});
                    if (bend[1] == 5) {
                        //console.log(bend, edgePaths[tup(...bend)]);
                    }
                } else {
                    // This 'edge' is actually a crossing, so let's break it up
                    // pathA goes bend0 -> bend1 -> bend2
                    // pathB goes in component order, too
                    let obend = vertPaths[bend[1]].idx;
                    let pathA = p.clone({insert: false});
                    let pathB = edgePaths[tup(...obend)].clone({insert: false});
                    let outIdxs = [bend[0], bend[2], obend[0], obend[2]];

                    // Delete edgePaths for this, since it's not really an edge
                    delete edgePaths[tup(...vertPaths[bend[1]].idx)];

                    let cross = pathA.getCrossings(pathB);
                    // Create a crossing effect (TODO: signage, robustness?)
                    let xing = new paper.Path.Circle(cross[0].point, .1);
                    //let q = p.divide(xing);
                    xing.fillColor = 'violet';
                    let ch = p.subtract(xing, {trace: false});
                    //ch.selected = true;
                    //let ch = xing.subtract(p);
                    ch.strokeColor = 'black';
                    ch.strokeWidth = 1;
                    ch.strokeScaling = false;
                    p.remove();
                    xing.remove();

                    xing = new paper.Path.Circle(cross[0].point, .01);
                    xing.remove();

                    //laststop = stop;
                    //continue;

                    // Create paths for each of the crossing directions,
                    // originating from the crossing (for consistency);
                    let cutPaths = pathA.subtract(xing, {trace: false});
                    if (cutPaths.children === undefined) {
                        //laststop = stop;
                        //continue;
                        //console.log("pathA", pathA);
                        //console.log("pathB", pathB);
                        //console.log(cross);
                        //console.log(cutPaths);
                        pathA.addTo(paper.project);
                        pathB.addTo(paper.project);
                        pathA.strokeColor = "blue";
                        pathB.strokeColor = "green";
                        cutPaths.addTo(paper.project);
                        cutPaths.fullySelected = true;
                        //console.log(cutPaths);
                        return;
                    }
                    let arc1 = cutPaths.children[0];
                    arc1.reverse();
                    let arc2 = cutPaths.children[1];

                    cutPaths = pathB.subtract(xing, {trace: false});

                    if (cutPaths.children === undefined) {
                        //laststop = stop;
                        //continue;
                        //console.log("pathA", pathA);
                        //console.log("pathB", pathB);
                        //console.log(cross);
                        //console.log(cutPaths);
                        pathA.addTo(paper.project);
                        pathB.addTo(paper.project);
                        pathA.strokeColor = "violet";
                        pathB.strokeColor = "orange";
                        cutPaths.addTo(paper.project);
                        cutPaths.fullySelected = true;
                        //console.log(cutPaths);
                        return;
                    }
                    let arc3 = cutPaths.children[0];
                    arc3.reverse();
                    let arc4 = cutPaths.children[1];

                    vertPaths[bend[1]] = {type: "cross",
                                          outs: outIdxs};

                    crssPaths[tup(bend[1], outIdxs[0])] = arc1;
                    crssPaths[tup(bend[1], outIdxs[1])] = arc2;
                    crssPaths[tup(bend[1], outIdxs[2])] = arc3;
                    crssPaths[tup(bend[1], outIdxs[3])] = arc4;
                }

                laststop = stop;
            }
            p = new paper.Path("M"+laststop+"L"+firststart);
            p.strokeColor = 'black';
            p.strokeWidth = 1;
            p.strokeScaling = false;

            //vertPaths[cmp[cmp.length-2]].join(p);
        }

        paper.project.insertLayer(0, new paper.Layer()).activate();
        //console.log(paper.project.activeLayer);
        //return;

        for (let fi = 0; fi < faces.length; fi++) {
            let face = faces[fi];
            let cface = [face[face.length-1]].concat(face).concat([face[0]]);

            if (face.exterior) {
                continue;
            }

            //console.log(vertPaths);
            let laststop, firststart;
            let pathStr;
            let segs;
            segs = [];
            for (let i = 0; i < cface.length-2; i++) {
                let bend = cface.slice(i, i+3);

                if (bend[1] in vertPaths) {
                    let vp = vertPaths[bend[1]];
                    if (vp.type == "edge") {
                        // Push the edge
                        if (bend[0] == vp.idx[0]) {
                            segs.splice(segs.length,0,...edgePaths[tup(...vp.idx)].segments);
                        } else {
                            let tempP = edgePaths[tup(...vp.idx)].clone({insert: false});
                            tempP.reverse();
                            segs.splice(segs.length,0,
                                        ...tempP.segments);
                        }
                    } else if (vp.type == "cross") {
                        // Push the two appropriate segments
                        // ---->crs (in)
                        let tempP = crssPaths[tup(bend[1], bend[0])].clone({insert: false});
                        tempP.reverse();
                        segs.splice(segs.length,0,
                                    ...tempP.segments);

                        // crs----> (out)
                        segs.splice(segs.length,0,
                                    ...crssPaths[tup(bend[1], bend[2])].segments);
                    }
                }

                for (let j = 0; j < vertPaths[cface[i]].length; j++) {
                    segs.splice(segs.length,0,...vertPaths[cface[i]][j].segments.slice(0));
                    p = new paper.Path(vertPaths[cface[i]][j].segments);
                    p.strokeColor = "blue";
                    p.strokeWidth = 1;
                    p.strokeScaling = false;
                }
                //break;
            }
            p = new paper.Path(segs);
            //console.log(pathStr);
            //console.log(segs);
            p.fillColor = "red";
            p.fillColor.alpha = .2;
            p.strokeScaling = false;
            p.onMouseEnter = function(e) {
                this.fillColor.alpha = .5;
            };
            p.onMouseLeave = function(e) {
                this.fillColor.alpha = .2;
            };
            p.onClick = function(e) {
                deleteMonogon(fi);
            }
            //p.fullySelected = true;
            //break;
        }
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

function drawRandomMapAsync(n_verts=20, n_comps=0, max_att=50, type=0) {
    n_verts = Number.parseInt(document.getElementById("random-ncross").value);
    n_comps = Number.parseInt(document.getElementById("random-ncomps").value);
    type = Number.parseInt(document.getElementById("random-type").value);
    cpWorker.postMessage({
        function: "setRandomLinkDiagram",
        arguments: [n_verts, n_comps, max_att, type]
    });
}

function stepMapUpdate() {
    cpWorker.postMessage({
        function: "stepUpdate",
        arguments: []
    });
}

function deleteMonogon(i) {
    cpWorker.postMessage({
        function: "deleteMonogon",
        arguments: [i]
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

    updatePDCode: function(pdcode) {
        document.getElementById("map_input").value = pdcode;
    },

    setLinkDiagram: function(link_diagram) {
        meshDraw.set_link_diagram(link_diagram);
    },

    finalizeLinkDiagram: function(link_diagram, faces, components) {
        meshDraw.set_link_diagram(link_diagram);
        meshDraw.finalize_link_diagram(link_diagram, faces, components);
    },

    raiseError: function(msg) {
        document.getElementById("map_errbox").textContent = msg;
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
let sigma = [[0, 10, 19, 9], [1, 12, 2, 11], [3, 13, 4, 14], [5, 16, 6, 15], [7, 17, 8, 18]];

//let sigma = [[3,5,4,1],[0,2,7,9],[8,6,11,10]];
//let sigma = [[7,1,3,5],[4,11,0,9],[12,8,6,13],[10,2,15,14]];

// Random composite
//let sigma = [[0, 57, 79, 58], [1, 59, 2, 60], [3, 78, 4, 77], [5, 11, 6, 12], [7, 9, 8, 10], [13, 16, 14, 15], [17, 56, 18, 55], [19, 54, 20, 53], [21, 36, 22, 35], [23, 25, 24, 26], [27, 29, 28, 30], [31, 46, 32, 45], [33, 39, 34, 40], [37, 52, 38, 51], [41, 44, 42, 43], [47, 50, 48, 49], [61, 72, 62, 71], [63, 73, 64, 74], [65, 68, 66, 67], [69, 76, 70, 75]];

//let sigma = [[1, 36, 2, 35], [0, 38, 75, 37], [3, 18, 4, 17], [5, 28, 6, 27], [7, 21, 8, 22], [9, 31, 10, 32], [11, 34, 12, 33], [13, 24, 14, 23], [15, 25, 16, 26], [19, 30, 20, 29], [39, 73, 40, 74], [41, 72, 42, 71], [43, 66, 44, 65], [45, 51, 46, 52], [47, 61, 48, 62], [49, 60, 50, 59], [53, 63, 54, 64], [55, 69, 56, 70], [57, 68, 58, 67]];

//let sigma = [[1, 3, 2, 4], [0, 60, 59, 119], [5, 15, 6, 16], [7, 84, 8, 83], [9, 101, 10, 102], [11, 92, 12, 91], [13, 89, 14, 90], [17, 62, 18, 61], [19, 79, 20, 80], [21, 70, 22, 69], [23, 42, 24, 41], [25, 28, 26, 27], [29, 32, 30, 31], [33, 52, 34, 51], [35, 37, 36, 38], [39, 49, 40, 50], [43, 54, 44, 53], [45, 48, 46, 47], [55, 67, 56, 68], [57, 82, 58, 81], [63, 77, 64, 78], [65, 72, 66, 71], [73, 75, 74, 76], [85, 88, 86, 87], [93, 99, 94, 100], [95, 98, 96, 97], [103, 118, 104, 117], [105, 112, 106, 111], [107, 109, 108, 110], [113, 115, 114, 116]];

// Complicated 2-link
// TODO: Fails orthogonal -- getting negative flow
//let sigma = [[1, 4, 2, 3], [0, 66, 295, 65], [5, 64, 6, 63], [7, 49, 8, 50], [9, 12, 10, 11], [13, 43, 14, 44], [15, 46, 16, 45], [17, 19, 18, 20], [21, 24, 22, 23], [25, 32, 26, 31], [27, 30, 28, 29], [33, 47, 34, 48], [35, 38, 36, 37], [39, 42, 40, 41], [51, 54, 52, 53], [55, 294, 56, 293], [57, 188, 58, 187], [59, 289, 60, 290], [61, 175, 62, 176], [67, 174, 68, 173], [69, 140, 70, 139], [71, 278, 72, 277], [73, 131, 74, 132], [75, 78, 76, 77], [79, 117, 80, 118], [81, 87, 82, 88], [83, 86, 84, 85], [89, 104, 90, 103], [91, 94, 92, 93], [95, 98, 96, 97], [99, 101, 100, 102], [105, 108, 106, 107], [109, 111, 110, 112], [113, 256, 114, 255], [115, 257, 116, 258], [119, 122, 120, 121], [123, 134, 124, 133], [125, 275, 126, 276], [127, 137, 128, 138], [129, 244, 130, 243], [135, 269, 136, 270], [141, 280, 142, 279], [143, 145, 144, 146], [147, 241, 148, 242], [149, 151, 150, 152], [153, 155, 154, 156], [157, 192, 158, 191], [159, 169, 160, 170], [161, 163, 162, 164], [165, 168, 166, 167], [171, 189, 172, 190], [177, 184, 178, 183], [179, 181, 180, 182], [185, 292, 186, 291], [193, 227, 194, 228], [195, 198, 196, 197], [199, 285, 200, 286], [201, 208, 202, 207], [203, 222, 204, 221], [205, 211, 206, 212], [209, 224, 210, 223], [213, 219, 214, 220], [215, 217, 216, 218], [225, 284, 226, 283], [229, 287, 230, 288], [231, 234, 232, 233], [235, 298, 236, 297], [237, 299, 238, 296], [239, 282, 240, 281], [245, 268, 246, 267], [247, 249, 248, 250], [251, 265, 252, 266], [253, 260, 254, 259], [261, 264, 262, 263], [271, 274, 272, 273]];

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

document.getElementById("map_random").onclick = function(ev) {
    try {
        meshDraw.clear();
        drawRandomMapAsync(20, 3, 50, 0);
    } catch(err) {
        document.querySelector(ev.target.dataset.errbox).textContent = "Error";
        console.log("Error:", err);
    }
};

document.getElementById("map_step").onclick = function(ev) {
    try {
        meshDraw.clear();
        stepMapUpdate();
    } catch(err) {
        
    }
}
