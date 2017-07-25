"use strict";

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

importScripts('lalolib/lalolib.js');

function least_squares(X /* : Matrix */, Y /* : Matrix */) /* : least_squares */{
    var betaHat = solve(mul(X, transpose(X)), mul(X, Y));

    return betaHat;
}

var Map4v = function () {
    function Map4v(verts) {
        _classCallCheck(this, Map4v);

        this.nv = verts.length;
        this.ne = this.nv * 2;
        this.na = this.nv * 4;

        this.arcs = [];
        this.edges = [];
        this.verts = [];

        // Edge i is canonically of the form [2i, 2i+1]
        for (var ei = 0; ei < this.ne; ei++) {
            //console.log(ei);
            this.new_arc(2 * ei);
            this.new_arc(2 * ei + 1);
            //console.log(this.arcs)

            this.set_edge(ei, [2 * ei, 2 * ei + 1]);
        }

        // Set verts by looking through verts
        for (var vi in verts) {
            this.set_vert(vi, verts[vi]);
        }

        this.generate_faces();
        this.generate_components();
    }

    _createClass(Map4v, [{
        key: "triangulate",
        value: function triangulate() {
            var _this = this;

            var bdry_face_i = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : undefined;

            var faces = this.faces;

            if (bdry_face_i === undefined) {
                var sizes = faces.map(function (f) {
                    return f.length;
                }, this);
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

            var triangles = [];

            /* Vertex indices correspond to components of the combinatorial map as follows:
               + [0:nv] correspond to actual vertices,
               + [nv:nv+ne] correspond to edge vertices,
               + [nv+ne:nv+ne+nf-1] correspond to (nonboundary) faces;
                Vertices are canonically indexed 0:nv+ne+nf-1. */
            var verts = range(0, this.nv + this.ne + faces.length - 1);

            // Boundary verts are simply determined by the vertices and edges in the boundary face.
            //console.log(faces);
            var bdy_verts = faces[bdry_face_i].map(function (a) {
                return a.vert;
            }, this).concat(faces[bdry_face_i].map(function (a) {
                return _this.nv + a.edge;
            }, this));
            // Remove the boundary face.
            faces.splice(bdry_face_i, 1);

            // Some edges correspond to arcs [vert->edge pairs]
            //console.log(this.arcs);
            var edges = this.arcs.map(function (arc) {
                return [arc.vert, _this.nv + arc.edge];
            }, this);

            // Loop through the faces, producing the rest of the edges and faces.
            for (var fi in faces) {
                var face = faces[fi];
                var fvi = this.nv + this.ne + parseInt(fi);
                // This face vertex has already been produced (nv+ne+fi)
                var _iteratorNormalCompletion = true;
                var _didIteratorError = false;
                var _iteratorError = undefined;

                try {
                    for (var _iterator = face[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
                        var arc = _step.value;

                        // We get a scaffolding edge face->arc.edge
                        edges.push([fvi, this.nv + arc.edge]);
                        // We get a scaffolding edge face->arc.vert
                        edges.push([fvi, arc.vert]);

                        // We get two faces; one corresponding to this arc,
                        triangles.push([fvi, this.nv + arc.edge, arc.vert]);
                        // and one corresponding to its edge opposite
                        var o_arc = this.edges[arc.edge][(arc.edgepos + 1) % 2];
                        triangles.push([fvi, this.nv + o_arc.edge, o_arc.vert]);
                    }
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
            }

            return [verts, bdy_verts, edges, triangles];
        }
    }, {
        key: "generate_faces",
        value: function generate_faces() {
            var left_arcs = new Set(this.arcs);
            var faces = [];

            while (left_arcs.size > 0) {
                var start_arc = Array.from(left_arcs).pop();

                var face = [];
                var arc = start_arc;
                var _failsafe = 0;
                do {
                    left_arcs.delete(arc);
                    face.push(arc);
                    //console.log(arc);

                    var o_arc = this.edges[arc.edge][(arc.edgepos + 1) % 2];
                    arc = this.verts[o_arc.vert][(o_arc.vertpos + 1) % 4];
                    //console.log(arc);
                    _failsafe += 1;
                    if (_failsafe > 500) {
                        console.log("Failure");
                        return faces;
                    }
                } while (arc != start_arc);

                face.reverse();
                faces.push(face);
            }

            this.faces = faces;
        }
    }, {
        key: "generate_components",
        value: function generate_components() {
            var one_orient = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : true;

            var left_arcs = new Set(this.arcs);

            var components = [];
            while (left_arcs.size > 0) {
                var start_arc = Array.from(left_arcs).pop();

                var component = this.component(start_arc);
                var _iteratorNormalCompletion2 = true;
                var _didIteratorError2 = false;
                var _iteratorError2 = undefined;

                try {
                    for (var _iterator2 = component[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                        var arc = _step2.value;

                        left_arcs.delete(arc);
                        if (one_orient) {
                            // Delete the other arc edge-opposite this one
                            left_arcs.delete(this.edges[arc.edge][(arc.edgepos + 1) % 2]);
                        }
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

                components.push(component);
            }
            this.components = components;
        }
    }, {
        key: "component",
        value: function component(arc) {
            var start_arc = arc;

            var component = [];
            do {
                component.push(arc);

                var o_arc = this.edges[arc.edge][(arc.edgepos + 1) % 2];
                arc = this.verts[o_arc.vert][(o_arc.vertpos + 2) % 4];
            } while (arc != start_arc);

            return component;
        }
    }, {
        key: "new_arc",
        value: function new_arc(idx) {
            this.arcs[idx] = { index: idx, edge: undefined, edgepos: undefined, vert: undefined, vertpos: undefined };
        }
    }, {
        key: "set_edge",
        value: function set_edge(idx, ais) {
            var _this2 = this;

            this.edges[idx] = ais.map(function (ai) {
                return _this2.arcs[ai];
            }, this);

            for (var i in ais) {
                this.arcs[ais[i]].edge = idx;
                this.arcs[ais[i]].edgepos = parseInt(i);
            }
        }
    }, {
        key: "set_vert",
        value: function set_vert(idx, ais) {
            var _this3 = this;

            this.verts[idx] = ais.map(function (ai) {
                return _this3.arcs[ai];
            }, this);

            for (var i in ais) {
                //console.log(i)
                //console.log(this.arcs)
                this.arcs[ais[i]].vert = parseInt(idx);
                this.arcs[ais[i]].vertpos = parseInt(i);
            }
        }
    }]);

    return Map4v;
}();

var TriangleMesh = function () {
    function TriangleMesh(verts, bdyverts, edges, faces) {
        _classCallCheck(this, TriangleMesh);

        this.verts = verts;
        this.bdyverts = bdyverts;
        this.edges = edges;
        this.faces = faces;
        this.bg_geom = "euclidean";
    }

    _createClass(TriangleMesh, [{
        key: "adjacent_edges",
        value: function adjacent_edges(vert) {
            var adj = [];
            var _iteratorNormalCompletion3 = true;
            var _didIteratorError3 = false;
            var _iteratorError3 = undefined;

            try {
                for (var _iterator3 = this.edges[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
                    var edge = _step3.value;

                    if (edge.indexOf(vert) >= 0) {
                        adj.push(edge);
                    }
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

            return adj;
        }
    }, {
        key: "adjacent_faces",
        value: function adjacent_faces(vert) {
            return this.faces.filter(function (f) {
                return f.indexOf(vert) >= 0;
            });
        }
    }, {
        key: "valence",
        value: function valence(vert) {
            return this.adjacent_faces(vert).length;
        }
    }, {
        key: "valences",
        value: function valences() {
            return this.verts.map(this.valence, this);
        }
    }, {
        key: "min_valence",
        value: function min_valence() {
            return min(this.valences());
        }
    }, {
        key: "chi",
        value: function chi() {
            // Euler characteristic
            return this.verts.length - this.edges.length + this.faces.length;
        }
    }]);

    return TriangleMesh;
}();

function partition_face(face) {
    return [[face[0], [face[1], face[2]]], [face[1], [face[2], face[0]]], [face[2], [face[0], face[1]]]];
}

var DiscreteRiemannMetric = function () {
    function DiscreteRiemannMetric(mesh, radius_map, edge_weights) {
        _classCallCheck(this, DiscreteRiemannMetric);

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

    _createClass(DiscreteRiemannMetric, [{
        key: "conf_factor",
        value: function conf_factor(gamma) {
            return log(gamma);
        }
    }, {
        key: "compute_length",
        value: function compute_length(edge) {
            // console.log(this.gamma, this.gamma.val);
            var g0 = this.gamma[edge[0]];
            var g1 = this.gamma[edge[1]];

            // Law of cosines
            return sqrt(2 * g0 * g1 * cos(this.phi.val[edge[0] * this.phi.n + edge[1]]) + Math.pow(g0, 2) + Math.pow(g1, 2));
        }
    }, {
        key: "length",
        value: function length(edge) {
            return this.l.val[edge[0] * this.l.n + edge[1]];
        }
    }, {
        key: "abc_for_vert",
        value: function abc_for_vert(face, vert) {
            var other_vs = face.filter(function (v) {
                return v != vert;
            });
            var edge_a = [vert, other_vs[0]];
            var edge_b = [vert, other_vs[1]];
            var edge_c = other_vs;

            //console.log(face, vert, edge_a, edge_b, edge_c)
            return [edge_a, edge_b, edge_c].map(this.length, this);
        }
    }, {
        key: "compute_angle",
        value: function compute_angle(face, vert) {
            var abc = this.abc_for_vert(face, vert);

            var ratio = (Math.pow(abc[0], 2) + Math.pow(abc[1], 2) - Math.pow(abc[2], 2)) / (2.0 * abc[0] * abc[1]);
            return Math.acos(ratio);
        }
    }, {
        key: "angle",
        value: function angle(face, vert) {
            var other_vs = face.filter(function (v) {
                return v != vert;
            });
            return this.theta[[vert, other_vs[0], other_vs[1]]];
        }
    }, {
        key: "curvature",
        value: function curvature(vert) {
            var _this4 = this;

            var adj_faces = this.mesh.adjacent_faces(vert);
            if (this.mesh.bdyverts.indexOf(vert) >= 0) {
                return Math.PI - sum(adj_faces.map(function (face) {
                    return _this4.angle(face, vert);
                }, this));
            } else {
                return 2 * Math.PI - sum(adj_faces.map(function (face) {
                    return _this4.angle(face, vert);
                }, this));
            }
        }
    }, {
        key: "update",
        value: function update() {
            this.gamma = exp(this.u);

            var _iteratorNormalCompletion4 = true;
            var _didIteratorError4 = false;
            var _iteratorError4 = undefined;

            try {
                for (var _iterator4 = this.mesh.edges[Symbol.iterator](), _step4; !(_iteratorNormalCompletion4 = (_step4 = _iterator4.next()).done); _iteratorNormalCompletion4 = true) {
                    var edge = _step4.value;

                    var l = this.compute_length(edge);
                    this.l.val[edge[0] * this.l.n + edge[1]] = l;
                    this.l.val[edge[1] * this.l.n + edge[0]] = l;
                }

                // Set angles using law of cosines
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

            var _iteratorNormalCompletion5 = true;
            var _didIteratorError5 = false;
            var _iteratorError5 = undefined;

            try {
                for (var _iterator5 = this.mesh.faces[Symbol.iterator](), _step5; !(_iteratorNormalCompletion5 = (_step5 = _iterator5.next()).done); _iteratorNormalCompletion5 = true) {
                    var face = _step5.value;
                    var _iteratorNormalCompletion6 = true;
                    var _didIteratorError6 = false;
                    var _iteratorError6 = undefined;

                    try {
                        for (var _iterator6 = partition_face(face)[Symbol.iterator](), _step6; !(_iteratorNormalCompletion6 = (_step6 = _iterator6.next()).done); _iteratorNormalCompletion6 = true) {
                            var part = _step6.value;

                            var theta = this.compute_angle(face, part[0]);
                            this.theta[[part[0], part[1][0], part[1][1]]] = theta;
                            this.theta[[part[0], part[1][1], part[1][0]]] = theta;
                        }
                    } catch (err) {
                        _didIteratorError6 = true;
                        _iteratorError6 = err;
                    } finally {
                        try {
                            if (!_iteratorNormalCompletion6 && _iterator6.return) {
                                _iterator6.return();
                            }
                        } finally {
                            if (_didIteratorError6) {
                                throw _iteratorError6;
                            }
                        }
                    }
                }

                // Set curvatures
            } catch (err) {
                _didIteratorError5 = true;
                _iteratorError5 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion5 && _iterator5.return) {
                        _iterator5.return();
                    }
                } finally {
                    if (_didIteratorError5) {
                        throw _iteratorError5;
                    }
                }
            }

            this.K = sparse(this.mesh.verts.map(this.curvature, this));
        }
    }, {
        key: "newton",
        value: function newton() {
            var target_K = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : null;
            var dt = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0.05;
            var thresh = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 1e-4;

            if (target_K == null) {}

            var g = new DiscreteRiemannMetric(this.mesh, this.gamma, this.phi);

            var K = g.K;
            var DeltaK = sub(target_K, K);

            var _failsafe = 0;
            while (max(abs(DeltaK)) > thresh) {
                var H = this.hessian();
                var deltau = least_squares(H, DeltaK);

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
    }, {
        key: "newton_async",
        value: function newton_async(f) {
            var target_K = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : null;
            var dt = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0.05;
            var thresh = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 1e-4;

            /* Asynchronous; passes resultant g to function f */
            if (target_K == null) {}

            var g = new DiscreteRiemannMetric(this.mesh, this.gamma, this.phi);

            var K = g.K;
            var DeltaK = sub(target_K, K);

            var _failsafe = 0;
            while (max(abs(DeltaK)) > thresh) {
                var H = this.hessian();
                var deltau = least_squares(H, DeltaK);

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
    }, {
        key: "tau2",
        value: function tau2(l_jk, g_j, g_k) {
            return .5 * (Math.pow(l_jk, 2) + Math.pow(g_j, 2) - Math.pow(g_k, 2));
        }
    }, {
        key: "face_area",
        value: function face_area(face) {
            var gamma = this.theta[face];
            var a = this.length([face[0], face[1]]);
            var b = this.length([face[0], face[2]]);

            return .5 * a * b * sin(gamma);
        }
    }, {
        key: "hessian",
        value: function hessian() {
            var n = this.mesh.verts.length;
            var H = zeros(n, n);
            var t = this.tau2;

            var _iteratorNormalCompletion7 = true;
            var _didIteratorError7 = false;
            var _iteratorError7 = undefined;

            try {
                for (var _iterator7 = this.mesh.faces[Symbol.iterator](), _step7; !(_iteratorNormalCompletion7 = (_step7 = _iterator7.next()).done); _iteratorNormalCompletion7 = true) {
                    var face = _step7.value;

                    var i = face[0];
                    var j = face[1];
                    var k = face[2];

                    var l_k = this.l.val[i * this.l.n + j];
                    var l_i = this.l.val[j * this.l.n + k];
                    var l_j = this.l.val[k * this.l.n + i];

                    var g_i = this.gamma[i];
                    var g_j = this.gamma[j];
                    var g_k = this.gamma[k];

                    var th_i = this.angle(face, i);
                    var th_j = this.angle(face, j);
                    var th_k = this.angle(face, k);

                    var A = this.face_area(face);

                    var L = diag([l_i, l_j, l_k]);
                    var D = array2mat([[0, t(l_i, g_j, g_k), t(l_i, g_k, g_j)], [t(l_j, g_i, g_k), 0, t(l_j, g_k, g_i)], [t(l_k, g_i, g_j), t(l_k, g_j, g_i), 0]]);

                    var Theta = cos(array2mat([[Math.PI, th_k, th_j], [th_k, Math.PI, th_i], [th_j, th_i, Math.PI]]));

                    var Tijk = mul(-.5 / A, mul(mul(mul(L, Theta), inv(L)), D));

                    var _arr = [0, 1, 2];
                    for (var _i = 0; _i < _arr.length; _i++) {
                        var rowi = _arr[_i];
                        var a = [i, j, k][rowi];
                        var _arr2 = [0, 1, 2];
                        for (var _i2 = 0; _i2 < _arr2.length; _i2++) {
                            var coli = _arr2[_i2];
                            var b = [i, j, k][coli];
                            H.val[a * H.n + b] += Tijk.val[rowi * Tijk.n + coli];
                        }
                    }
                }
            } catch (err) {
                _didIteratorError7 = true;
                _iteratorError7 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion7 && _iterator7.return) {
                        _iterator7.return();
                    }
                } finally {
                    if (_didIteratorError7) {
                        throw _iteratorError7;
                    }
                }
            }

            return H;
        }
    }], [{
        key: "from_triangle_mesh",
        value: function from_triangle_mesh(mesh) {
            var n = mesh.verts.length;
            var gamma = mesh.verts.map(function () {
                return 1;
            }, this);
            var phi = zeros(n, n);

            var _iteratorNormalCompletion8 = true;
            var _didIteratorError8 = false;
            var _iteratorError8 = undefined;

            try {
                for (var _iterator8 = mesh.edges[Symbol.iterator](), _step8; !(_iteratorNormalCompletion8 = (_step8 = _iterator8.next()).done); _iteratorNormalCompletion8 = true) {
                    var edge = _step8.value;

                    phi.val[edge[0] * phi.n + edge[1]] = 0;
                    phi.val[edge[1] * phi.n + edge[0]] = 0;
                }
            } catch (err) {
                _didIteratorError8 = true;
                _iteratorError8 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion8 && _iterator8.return) {
                        _iterator8.return();
                    }
                } finally {
                    if (_didIteratorError8) {
                        throw _iteratorError8;
                    }
                }
            }

            return new DiscreteRiemannMetric(mesh, gamma, phi);
        }
    }]);

    return DiscreteRiemannMetric;
}();

function adjacent_faces_and_edge(faces, face) {
    var adj = [];

    var _iteratorNormalCompletion9 = true;
    var _didIteratorError9 = false;
    var _iteratorError9 = undefined;

    try {
        for (var _iterator9 = faces[Symbol.iterator](), _step9; !(_iteratorNormalCompletion9 = (_step9 = _iterator9.next()).done); _iteratorNormalCompletion9 = true) {
            var aface = _step9.value;

            var e = aface.filter(function (v) {
                return face.indexOf(v) >= 0;
            });
            if (e.length == 2) {
                if (face[(face.indexOf(e[0]) + 1) % 3] == e[1]) {
                    adj.push([aface, [e[1], e[0]]]);
                } else {
                    adj.push([aface, e]);
                }
            }
        }
    } catch (err) {
        _didIteratorError9 = true;
        _iteratorError9 = err;
    } finally {
        try {
            if (!_iteratorNormalCompletion9 && _iterator9.return) {
                _iterator9.return();
            }
        } finally {
            if (_didIteratorError9) {
                throw _iteratorError9;
            }
        }
    }

    return adj;
}

function adj_or_faces_and_edge(faces, face) {
    var adj = [];

    var _iteratorNormalCompletion10 = true;
    var _didIteratorError10 = false;
    var _iteratorError10 = undefined;

    try {
        for (var _iterator10 = faces[Symbol.iterator](), _step10; !(_iteratorNormalCompletion10 = (_step10 = _iterator10.next()).done); _iteratorNormalCompletion10 = true) {
            var aface = _step10.value;

            if (aface.indexOf(face[0]) >= 0 && aface.indexOf(face[1]) >= 0) {
                adj.push([aface, [face[0], face[1]]]);
            } else if (aface.indexOf(face[1]) >= 0 && aface.indexOf(face[2]) >= 0) {
                adj.push([aface, [face[1], face[2]]]);
            } else if (aface.indexOf(face[2]) >= 0 && aface.indexOf(face[0]) >= 0) {
                adj.push([aface, [face[2], face[0]]]);
            }
        }
    } catch (err) {
        _didIteratorError10 = true;
        _iteratorError10 = err;
    } finally {
        try {
            if (!_iteratorNormalCompletion10 && _iterator10.return) {
                _iterator10.return();
            }
        } finally {
            if (_didIteratorError10) {
                throw _iteratorError10;
            }
        }
    }

    return adj;
}

function orient_faces(faces) {
    var edges = [];
    var oriented = [];
    var to_orient = faces.slice();
    var adj_queue = new Set();

    var f_0 = to_orient.pop();
    edges.concat([[f_0[0], f_0[1]], [f_0[1], f_0[2]], [f_0[2], f_0[0]]]);
    oriented.push(f_0);

    var _iteratorNormalCompletion11 = true;
    var _didIteratorError11 = false;
    var _iteratorError11 = undefined;

    try {
        for (var _iterator11 = adjacent_faces_and_edge(to_orient, f_0)[Symbol.iterator](), _step11; !(_iteratorNormalCompletion11 = (_step11 = _iterator11.next()).done); _iteratorNormalCompletion11 = true) {
            var adj_pair = _step11.value;

            adj_queue.add(adj_pair);
        }
    } catch (err) {
        _didIteratorError11 = true;
        _iteratorError11 = err;
    } finally {
        try {
            if (!_iteratorNormalCompletion11 && _iterator11.return) {
                _iterator11.return();
            }
        } finally {
            if (_didIteratorError11) {
                throw _iteratorError11;
            }
        }
    }

    var _failsafe = 0;
    while (to_orient.length > 0 && adj_queue.size > 0) {
        //console.log(adj_queue, to_orient)
        var _adj_pair = Array.from(adj_queue).pop();
        adj_queue.delete(_adj_pair);

        //console.log(e);
        var F = _adj_pair[0];
        var e = _adj_pair[1];

        _failsafe += 1;
        if (_failsafe > 5000) {
            console.log("Too much time spent orienting faces...");
            return oriented;
        }

        if (to_orient.indexOf(F) < 0) {
            continue;
        }

        var v_i = void 0,
            v_j = void 0,
            v_k = void 0;
        if (e.indexOf(F[0]) < 0) {
            v_k = F[0];
        } else if (e.indexOf(F[1]) < 0) {
            v_k = F[1];
        } else {
            v_k = F[2];
        }
        v_i = e[0];
        v_j = e[1];

        var i = void 0,
            j = void 0,
            k = void 0;
        if (edges.includes([v_i, v_j])) {
            i = v_j;j = v_i;k = v_k;
        } else {
            i = v_i;j = v_j;k = v_k;
        }

        edges.concat([[i, j], [j, k], [k, i]]);
        oriented.push([i, j, k]);

        to_orient.splice(to_orient.indexOf(F), 1);
        var _iteratorNormalCompletion12 = true;
        var _didIteratorError12 = false;
        var _iteratorError12 = undefined;

        try {
            for (var _iterator12 = adjacent_faces_and_edge(to_orient, [i, j, k])[Symbol.iterator](), _step12; !(_iteratorNormalCompletion12 = (_step12 = _iterator12.next()).done); _iteratorNormalCompletion12 = true) {
                var _adj_pair2 = _step12.value;

                adj_queue.add(_adj_pair2);
            }
        } catch (err) {
            _didIteratorError12 = true;
            _iteratorError12 = err;
        } finally {
            try {
                if (!_iteratorNormalCompletion12 && _iterator12.return) {
                    _iterator12.return();
                }
            } finally {
                if (_didIteratorError12) {
                    throw _iteratorError12;
                }
            }
        }
    }

    return oriented;
}

function u_theta(theta) {
    return [cos(theta), sin(theta)];
}

function embed_faces(g) {
    var mesh = g.mesh;
    var x = mul(-30, ones(g.n, 2));
    var phi = {};
    var pi = Math.PI;

    var faces = orient_faces(mesh.faces);
    var to_embed = faces.slice();
    var embed_queue = new Set();

    var f_0 = to_embed.pop();
    var i = f_0[0],
        j = f_0[1],
        k = f_0[2];

    x.val[i * x.n + 0] = 0;
    x.val[i * x.n + 1] = 0;
    x.val[j * x.n + 0] = g.length([i, j]);
    x.val[j * x.n + 1] = 0;

    var phi_ik = g.angle(f_0, i) % (2 * pi);
    var g_ik = g.length([i, k]);

    x.val[k * x.n + 0] = g_ik * cos(phi_ik);
    x.val[k * x.n + 1] = g_ik * sin(phi_ik);

    var phi_jk = (pi - g.angle(f_0, j)) % (2 * pi);

    phi[[i, j]] = 0;
    phi[[j, i]] = pi;
    phi[[j, k]] = phi_jk;
    phi[[k, j]] = pi + phi_jk;
    phi[[i, k]] = phi_ik;
    phi[[k, i]] = pi + phi_ik;

    var _iteratorNormalCompletion13 = true;
    var _didIteratorError13 = false;
    var _iteratorError13 = undefined;

    try {
        for (var _iterator13 = adj_or_faces_and_edge(to_embed, f_0)[Symbol.iterator](), _step13; !(_iteratorNormalCompletion13 = (_step13 = _iterator13.next()).done); _iteratorNormalCompletion13 = true) {
            var adj_pair = _step13.value;

            embed_queue.add(adj_pair);
        }
    } catch (err) {
        _didIteratorError13 = true;
        _iteratorError13 = err;
    } finally {
        try {
            if (!_iteratorNormalCompletion13 && _iterator13.return) {
                _iterator13.return();
            }
        } finally {
            if (_didIteratorError13) {
                throw _iteratorError13;
            }
        }
    }

    var _failsafe = 0;
    while (to_embed.length > 0) {
        var _adj_pair3 = Array.from(embed_queue).pop();
        embed_queue.delete(_adj_pair3);
        var F = _adj_pair3[0];
        var e = _adj_pair3[1];

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

        var _i3 = e[0],
            _j = e[1],
            _k = void 0;
        if (F[(F.indexOf(_i3) + 1) % 3] != _j) {
            _k = F[(F.indexOf(_i3) + 1) % 3];
            _j = e[0], _i3 = e[1];
        } else {
            _k = F[(F.indexOf(_i3) + 2) % 3];
        }

        //console.log(F, e, i, j, k);
        //console.log(x.subset(math.index([i,j,k], [0,1])).toArray());

        // We already know x[i], x[j]. We only have to find x[k].
        phi_ik = (phi[[_i3, _j]] + g.angle(F, _i3)) % (2 * pi);
        phi_jk = (phi[[_j, _i3]] - g.angle(F, _j)) % (2 * pi);

        //if (_failsafe == 11) {
        // phi_ik = (phi[[i,j]] - g.angle(F, i)) % (2*pi);
        //}
        if (!([_j, _k] in phi) && !([_i3, _k] in phi)) {
            g_ik = g.length([_i3, _k]);
            x.val[_k * x.n + 0] = x.val[_i3 * x.n + 0] + g_ik * cos(phi_ik);
            x.val[_k * x.n + 1] = x.val[_i3 * x.n + 1] + g_ik * sin(phi_ik);
        }

        //console.log(g.length([i,k]))
        //console.log(g.angle(F, i));
        //console.log(phi[[i,j]]);
        //console.log(x.subset(math.index(i, [0,1])).toArray()[0])
        //console.log(phi_ik)
        //console.log(u_theta(phi_ik))

        if (!([_j, _k] in phi)) {
            phi[[_j, _k]] = phi_jk;
            phi[[_k, _j]] = pi + phi_jk;
        }

        if (!([_i3, _k] in phi)) {
            phi[[_i3, _k]] = phi_ik;
            phi[[_k, _i3]] = pi + phi_ik;
        }

        to_embed.splice(to_embed.indexOf(F), 1);
        var _iteratorNormalCompletion14 = true;
        var _didIteratorError14 = false;
        var _iteratorError14 = undefined;

        try {
            for (var _iterator14 = adj_or_faces_and_edge(to_embed, F)[Symbol.iterator](), _step14; !(_iteratorNormalCompletion14 = (_step14 = _iterator14.next()).done); _iteratorNormalCompletion14 = true) {
                var _adj_pair4 = _step14.value;

                embed_queue.add(_adj_pair4);
            }
        } catch (err) {
            _didIteratorError14 = true;
            _iteratorError14 = err;
        } finally {
            try {
                if (!_iteratorNormalCompletion14 && _iterator14.return) {
                    _iterator14.return();
                }
            } finally {
                if (_didIteratorError14) {
                    throw _iteratorError14;
                }
            }
        }
    }

    return [x, faces, phi];
}

onmessage = function onmessage(e) {
    var sigma = e.data[0];
    var cross_bend = e.data[1];
    var tstart = Date.now();

    var m4v = new Map4v(sigma);

    var triangulation = m4v.triangulate();

    //console.log(triangulation);

    var testMesh = new TriangleMesh(triangulation[0], triangulation[1], triangulation[2], triangulation[3]);

    var cpmetric = DiscreteRiemannMetric.from_triangle_mesh(testMesh);

    var K = zeros(testMesh.verts.length);
    //K.subset(math.index(testMesh.bdyverts, 0),
    //         testMesh.bdyverts.map(function(b) { return 2*Math.PI/testMesh.bdyverts.length; }));

    // Set boundary crossing target curvature
    var bdyCross = testMesh.bdyverts.filter(function (vi) {
        return vi < m4v.nv;
    });
    var bdyEdge = testMesh.bdyverts.filter(function (vi) {
        return vi >= m4v.nv && vi < m4v.nv + m4v.ne;
    });
    var fac = 8; // Inverse of how concave crossing vertices should be imbedded
    var _iteratorNormalCompletion15 = true;
    var _didIteratorError15 = false;
    var _iteratorError15 = undefined;

    try {
        for (var _iterator15 = bdyCross[Symbol.iterator](), _step15; !(_iteratorNormalCompletion15 = (_step15 = _iterator15.next()).done); _iteratorNormalCompletion15 = true) {
            var bci = _step15.value;

            K[bci] = -Math.PI / fac;
        }
    } catch (err) {
        _didIteratorError15 = true;
        _iteratorError15 = err;
    } finally {
        try {
            if (!_iteratorNormalCompletion15 && _iterator15.return) {
                _iterator15.return();
            }
        } finally {
            if (_didIteratorError15) {
                throw _iteratorError15;
            }
        }
    }

    var _iteratorNormalCompletion16 = true;
    var _didIteratorError16 = false;
    var _iteratorError16 = undefined;

    try {
        for (var _iterator16 = bdyEdge[Symbol.iterator](), _step16; !(_iteratorNormalCompletion16 = (_step16 = _iterator16.next()).done); _iteratorNormalCompletion16 = true) {
            var bei = _step16.value;

            K[bei] = (2 * Math.PI + bdyCross.length * Math.PI / fac) / bdyEdge.length;
        }
    } catch (err) {
        _didIteratorError16 = true;
        _iteratorError16 = err;
    } finally {
        try {
            if (!_iteratorNormalCompletion16 && _iterator16.return) {
                _iterator16.return();
            }
        } finally {
            if (_didIteratorError16) {
                throw _iteratorError16;
            }
        }
    }

    var flat_poly = cpmetric.newton(K, 1, 5e-2);

    var embedding = embed_faces(flat_poly);

    postMessage({
        flat_poly: flat_poly,
        embedding: embedding,
        m4v: m4v });
    console.log("Computation finished in: " + (Date.now() - tstart) + " milliseconds");
};