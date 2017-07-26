'use strict';

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

//let lab = new Lalolab();

function least_squares(X /* : Matrix */, Y /* : Matrix */) /* : least_squares */{
    var betaHat = solve(mul(X, transpose(X)), mul(X, Y));

    return betaHat;
}

var MeshEdge = function MeshEdge(parent, start, stop) {
    _classCallCheck(this, MeshEdge);

    this.parent = parent;
    this.start = start;
    this.stop = stop;

    this.svg = this.parent.edgeG.line(start.x(), start.y(), stop.x(), stop.y());
    this.svg.addClass("scaffold");
};

var MeshNode = function () {
    function MeshNode(parent, x, y) {
        _classCallCheck(this, MeshNode);

        this.parent = parent;
        this.svg = this.parent.nodeG.circle(x, y, .3);
        this.svg.addClass('plnode');

        this.svg.node.addEventListener('click', this.onClick.bind(this));
    }

    _createClass(MeshNode, [{
        key: 'set_obj',
        value: function set_obj(i, obj) {
            this.i = i;
            this.obj = obj;
        }
    }, {
        key: 'set_r',
        value: function set_r(r) {
            this.svg.attr({ 'r': r });
        }
    }, {
        key: 'x',
        value: function x() {
            return this.svg.attr('cx');
        }
    }, {
        key: 'y',
        value: function y() {
            return this.svg.attr('cy');
        }
    }, {
        key: 'onClick',
        value: function onClick(e) {
            console.log(this);
            console.log(this.obj);
            if (this.obj.length == 2) {
                this.parent.delete_face(this.i);
            }
        }
    }]);

    return MeshNode;
}();

var MeshDraw = function () {
    function MeshDraw(div) {
        _classCallCheck(this, MeshDraw);

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

    _createClass(MeshDraw, [{
        key: 'clear',
        value: function clear() {
            this.nodes = [];
            this.edges = {};

            this.edgeG.clear();
            this.nodeG.clear();
            this.knotG.clear();
        }
    }, {
        key: 'set_embedding',
        value: function set_embedding(g /*metric*/, embedding, map4v /*original map*/, conv) {
            var points = embedding[0];
            var faces = embedding[1];
            var phi = embedding[2];
            /* set the embedding */
            this.g = g;
            this.map4v = map4v;

            // alter the viewbox to fit
            var min_x = min(get(points, range(), 0));
            var min_y = min(get(points, range(), 1));
            var max_x = max(get(points, range(), 0));
            var max_y = max(get(points, range(), 1));

            var wid = max_x - min_x;
            var hgt = max_y - min_y;

            var dx = wid * 0.05;
            var dy = hgt * 0.05;

            this.draw.attr({ viewBox: [min_x - dx, min_y - dy, wid + 2 * dx, hgt + 2 * dy].join(",") });

            //console.log(points);
            for (var i = 0; i < points.m; i++) {
                var node = this.add_node(i, points.val[i * points.n], points.val[i * points.n + 1]);
                node.set_r(this.g.gamma[this.nodes.length - 1]);
            }

            var _iteratorNormalCompletion = true;
            var _didIteratorError = false;
            var _iteratorError = undefined;

            try {
                for (var _iterator = conv.comps[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
                    var comp = _step.value;

                    for (var pi = 0; pi < comp.length; pi++) {
                        var edge = [comp[pi], comp[(pi + 1) % comp.length]];
                        console.log(edge);
                        if (edge[0] in this.nodes && edge[1] in this.nodes) {
                            var mesh_edge = this.add_edge(edge[0], edge[1]);
                            mesh_edge.svg.addClass('edge');
                        }
                    }
                }

                //console.log(this.map4v.faces);
            } catch (err) {
                _didIteratorError = true;
                _iteratorError = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion && _iterator.return) {
                        _iterator.return();
                    }
                } finally {
                    if (_didIteratorError) {
                        throw _iteratorError;
                    }
                }
            }

            for (var fi in this.map4v.faces) {
                //console.log(fi);
                //let mesh_face = this.add_face(parseInt(fi));
            }

            var ci = 0;
            var _iteratorNormalCompletion2 = true;
            var _didIteratorError2 = false;
            var _iteratorError2 = undefined;

            try {
                for (var _iterator2 = map4v.components[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                    //let knot = this.add_component(component, map4v, points);
                    //knot.addClass('q'+ci+"-9");
                    //ci += 1;

                    var component = _step2.value;
                }
            } catch (err) {
                _didIteratorError2 = true;
                _iteratorError2 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion2 && _iterator2.return) {
                        _iterator2.return();
                    }
                } finally {
                    if (_didIteratorError2) {
                        throw _iteratorError2;
                    }
                }
            }
        }
    }, {
        key: 'delete_face',
        value: function delete_face(i) {
            /* Deleting a face (valid curve-preserving faces to delete are a bigon *
             * or a monogon) is something of an ordeal: We want to delete only the
             * local triangles which the face affects, modifying the underlying
             * triangulation only locally, so that we can keep the old discrete
             * riemannian metric, hopefully so that the resulting embedding looks
             * very similar. Is there a proof that this will happen? */
        }
    }, {
        key: 'add_node',
        value: function add_node(i, x, y) {
            var node = new MeshNode(this, x, y);
            this.nodes[i] = node;
            return node;
        }
    }, {
        key: 'add_edge',
        value: function add_edge(i, j) {
            var edge = new MeshEdge(this, this.nodes[i], this.nodes[j]);
            this.edges[[i, j]] = edge;
            this.edges[[j, i]] = edge;
            return edge;
        }
    }, {
        key: 'add_face',
        value: function add_face(i) {
            var face_node = this.nodes[this.map4v.nv + this.map4v.ne + i];
            face_node.svg.addClass("face");
            face_node.set_obj(i, this.map4v.faces[i]);
        }
    }, {
        key: 'add_component',
        value: function add_component(component, map4v, points) {
            // A path of the form anchor, control, anchor...
            var path = [];
            var pts = points;

            var root_arc = component[0];

            var last_ep = undefined;
            var _iteratorNormalCompletion3 = true;
            var _didIteratorError3 = false;
            var _iteratorError3 = undefined;

            try {
                for (var _iterator3 = component[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
                    var arc = _step3.value;

                    var _vi = arc.vert;
                    var _ei = map4v.nv + arc.edge;

                    var _vp = [pts.val[_vi * pts.n], pts.val[_vi * pts.n + 1]];
                    var _ep = [pts.val[_ei * pts.n], pts.val[_ei * pts.n + 1]];

                    if (last_ep != undefined) {
                        // We really want to push 'center' points as anchor points and
                        // real points as control points
                        var _dp2 = sub(last_ep, _vp);

                        // Find center point by linear interpolation; whole way across is
                        // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
                        var _whole_l2 = norm(_dp2);

                        // We only want to go this.g.gamma[vi]/|dp| of the way:
                        var _mid_dp2 = mul(this.g.gamma[_vi] / _whole_l2, _dp2);
                        var _mp2 = add(_vp, _mid_dp2);

                        path.push(_mp2);
                        path.push(_vp);
                    }

                    // We really want to push 'center' points as anchor points and
                    // real points as control points
                    var _dp = sub(_ep, _vp); // ep - vp; notice vp + dp = ep

                    // Find center point by linear interpolation; whole way across is
                    // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
                    var _whole_l = norm(_dp);

                    // We only want to go this.g.gamma[vi]/|dp| of the way:
                    var _mid_dp = mul(this.g.gamma[_vi] / _whole_l, _dp);
                    var _mp = add(_vp, _mid_dp);

                    //path.push(vp);
                    path.push(_mp);
                    path.push(_ep);

                    last_ep = _ep;
                }
            } catch (err) {
                _didIteratorError3 = true;
                _iteratorError3 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion3 && _iterator3.return) {
                        _iterator3.return();
                    }
                } finally {
                    if (_didIteratorError3) {
                        throw _iteratorError3;
                    }
                }
            }

            var vi = root_arc.vert;
            var ei = map4v.nv + root_arc.edge;

            var vp = [pts.val[vi * pts.n], pts.val[vi * pts.n + 1]];
            var ep = [pts.val[ei * pts.n], pts.val[vi * pts.n + 1]];

            // We really want to push 'center' points as anchor points and
            // real points as control points
            var dp = sub(last_ep, vp);

            // Find center point by linear interpolation; whole way across is
            // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
            var whole_l = norm(dp);

            // We only want to go this.g.gamma[vi]/|dp| of the way:
            var mid_dp = mul(this.g.gamma[vi] / whole_l, dp);
            var mp = add(vp, mid_dp);

            path.push(mp);
            path.push(vp);

            path.push(path[0]);

            var pathStr = "M";
            pathStr += path[0].join(",");

            var idx = 0;
            var _iteratorNormalCompletion4 = true;
            var _didIteratorError4 = false;
            var _iteratorError4 = undefined;

            try {
                for (var _iterator4 = path.slice(1)[Symbol.iterator](), _step4; !(_iteratorNormalCompletion4 = (_step4 = _iterator4.next()).done); _iteratorNormalCompletion4 = true) {
                    var pt = _step4.value;

                    if (idx % 2 == 0) {
                        pathStr += "Q";
                    } else {
                        pathStr += " ";
                    }
                    idx += 1;
                    pathStr += pt.join(",");
                }

                //console.log(pathStr);
            } catch (err) {
                _didIteratorError4 = true;
                _iteratorError4 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion4 && _iterator4.return) {
                        _iterator4.return();
                    }
                } finally {
                    if (_didIteratorError4) {
                        throw _iteratorError4;
                    }
                }
            }

            var knot = this.knotG.path(pathStr + "Z");
            knot.addClass('knot');

            return knot;
        }
    }]);

    return MeshDraw;
}();

var meshDraw = new MeshDraw("#knot-draw");

var cpWorker = new Worker("js/cp_worker.js");

function drawMapAsync(sigma) {
    var cross_bend = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 8;

    cpWorker.postMessage([sigma, cross_bend]);
}

cpWorker.onmessage = function (ev) {
    //console.log(ev.data.m4v);
    meshDraw.clear();
    meshDraw.set_embedding(ev.data.flat_poly, ev.data.embedding, ev.data.m4v, ev.data.conv);
};

// Trefoil
//let sigma = [[0, 6, 11, 5], [1, 8, 2, 7], [3, 9, 4, 10]];

// Twist
//let sigma = [[0, 1, 3, 2]];

// Monogon in internal face
var sigma = [[1, 8, 2, 7], [0, 14, 15, 13], [3, 9, 4, 10], [5, 12, 6, 11]];

// More complicated
//let sigma = [[0, 41, 99, 42], [1, 47, 2, 48], [3, 34, 4, 33], [5, 96, 6, 95], [7, 13, 8, 14], [9, 68, 10, 67], [11, 69, 12, 70], [15, 62, 16, 61], [17, 59, 18, 60], [19, 85, 20, 86], [21, 27, 22, 28], [23, 82, 24, 81], [25, 83, 26, 84], [29, 88, 30, 87], [31, 93, 32, 94], [35, 46, 36, 45], [37, 43, 38, 44], [39, 98, 40, 97], [49, 76, 50, 75], [51, 73, 52, 74], [53, 72, 54, 71], [55, 65, 56, 66], [57, 64, 58, 63], [77, 92, 78, 91], [79, 89, 80, 90]];

// Even moreso
//let sigma = [[1, 287, 2, 288], [0, 289, 299, 290], [3, 126, 4, 125], [5, 16, 6, 15], [7, 13, 8, 14], [9, 255, 10, 256], [11, 254, 12, 253], [17, 127, 18, 128], [19, 282, 20, 281], [21, 132, 22, 131], [23, 277, 24, 278], [25, 235, 26, 236], [27, 222, 28, 221], [29, 223, 30, 224], [31, 38, 32, 37], [33, 228, 34, 227], [35, 225, 36, 226], [39, 230, 40, 229], [41, 156, 42, 155], [43, 157, 44, 158], [45, 164, 46, 163], [47, 165, 48, 166], [49, 56, 50, 55], [51, 198, 52, 197], [53, 195, 54, 196], [57, 243, 58, 244], [59, 242, 60, 241], [61, 231, 62, 232], [63, 234, 64, 233], [65, 248, 66, 247], [67, 249, 68, 250], [69, 276, 70, 275], [71, 133, 72, 134], [73, 83, 74, 84], [75, 82, 76, 81], [77, 283, 78, 284], [79, 286, 80, 285], [85, 295, 86, 296], [87, 294, 88, 293], [89, 291, 90, 292], [91, 298, 92, 297], [93, 136, 94, 135], [95, 109, 96, 110], [97, 108, 98, 107], [99, 122, 100, 121], [101, 259, 102, 260], [103, 262, 104, 261], [105, 119, 106, 120], [111, 273, 112, 274], [113, 272, 114, 271], [115, 265, 116, 266], [117, 264, 118, 263], [123, 137, 124, 138], [129, 280, 130, 279], [139, 257, 140, 258], [141, 268, 142, 267], [143, 269, 144, 270], [145, 252, 146, 251], [147, 237, 148, 238], [149, 216, 150, 215], [151, 217, 152, 218], [153, 220, 154, 219], [159, 210, 160, 209], [161, 211, 162, 212], [167, 190, 168, 189], [169, 187, 170, 188], [171, 205, 172, 206], [173, 200, 174, 199], [175, 201, 176, 202], [177, 204, 178, 203], [179, 186, 180, 185], [181, 191, 182, 192], [183, 194, 184, 193], [207, 214, 208, 213], [239, 245, 240, 246]];

drawMapAsync(sigma);

document.getElementById("map_submit").onclick = function (ev) {
    document.querySelector(ev.target.dataset.errbox).textContent = "";
    try {
        drawMapAsync(JSON.parse(document.querySelector(ev.target.dataset.target).value));
    } catch (err) {
        document.querySelector(ev.target.dataset.errbox).textContent = "Error";
        console.log("Error:", err);
    }
};
//# sourceMappingURL=index.js.map