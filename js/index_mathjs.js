'use strict';

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function least_squares(X /* : Matrix */, Y /* : Matrix */) /* : least_squares */{
    var XX = math.transpose(X);

    var YY = Y;
    var betaHat = math.lusolve(math.multiply(math.transpose(XX), XX), math.multiply(math.transpose(XX), YY));

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
    }

    _createClass(MeshNode, [{
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
        value: function onClick(e) {}
    }]);

    return MeshNode;
}();

var MeshDraw = function () {
    function MeshDraw(div) {
        _classCallCheck(this, MeshDraw);

        this.nodes = [];
        this.edges = {};

        this.draw = Snap(div);
        this.edgeG = this.draw.g();
        this.nodeG = this.draw.g();
        this.knotG = this.draw.g();
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
        value: function set_embedding(g /*metric*/, embedding, map4v /*original map*/) {
            var points = embedding[0];
            var faces = embedding[1];
            var phi = embedding[2];
            /* set the embedding */
            this.g = g;

            // alter the viewbox to fit
            var min_x = math.min(math.subset(points, math.index(math.range(0, math.size(points).toArray()[0]), 0)));
            var max_x = math.max(math.subset(points, math.index(math.range(0, math.size(points).toArray()[0]), 0)));
            var min_y = math.min(math.subset(points, math.index(math.range(0, math.size(points).toArray()[0]), 1)));
            var max_y = math.max(math.subset(points, math.index(math.range(0, math.size(points).toArray()[0]), 1)));

            var wid = max_x - min_x;
            var hgt = max_y - min_y;

            var dx = wid * 0.05;
            var dy = hgt * 0.05;

            this.draw.attr({ viewBox: [min_x - dx, min_y - dy, wid + 2 * dx, hgt + 2 * dy].join(",") });

            var i = 0;
            var _iteratorNormalCompletion = true;
            var _didIteratorError = false;
            var _iteratorError = undefined;

            try {
                for (var _iterator = points.toArray()[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
                    var point = _step.value;

                    var node = this.add_node(i, point[0], point[1]);
                    node.set_r(this.g.gamma[this.nodes.length - 1]);
                    i += 1;
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

            var _iteratorNormalCompletion2 = true;
            var _didIteratorError2 = false;
            var _iteratorError2 = undefined;

            try {
                for (var _iterator2 = this.g.mesh.edges[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                    var edge = _step2.value;

                    if (edge[0] in this.nodes && edge[1] in this.nodes) {
                        var mesh_edge = this.add_edge(edge[0], edge[1]);

                        if (Math.max.apply(Math, _toConsumableArray(edge)) < map4v.nv + map4v.ne) {
                            // This edge has no face connection, and so is not scaffolding
                            mesh_edge.svg.addClass('edge');
                        }
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

            var root_arc = map4v.arcs[0];
            var component = map4v.component(root_arc);

            // A path of the form anchor, control, anchor...
            var path = [];
            var pts = points.toArray();

            var last_ep = undefined;
            var _iteratorNormalCompletion3 = true;
            var _didIteratorError3 = false;
            var _iteratorError3 = undefined;

            try {
                for (var _iterator3 = component[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
                    var arc = _step3.value;

                    var _vi = arc.vert;
                    var _ei = map4v.nv + arc.edge;

                    var _vp = pts[_vi];
                    var _ep = pts[_ei];

                    if (last_ep != undefined) {
                        // We really want to push 'center' points as anchor points and
                        // real points as control points
                        var _dp2 = math.subtract(last_ep, _vp);

                        // Find center point by linear interpolation; whole way across is
                        // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
                        var _whole_l2 = math.norm(_dp2);

                        // We only want to go this.g.gamma[vi]/|dp| of the way:
                        var _mid_dp2 = math.multiply(this.g.gamma[_vi] / _whole_l2, _dp2);
                        var _mp2 = math.add(_vp, _mid_dp2);

                        path.push(_mp2);
                        path.push(_vp);
                    }

                    // We really want to push 'center' points as anchor points and
                    // real points as control points
                    var _dp = math.subtract(_ep, _vp); // ep - vp; notice vp + dp = ep

                    // Find center point by linear interpolation; whole way across is
                    // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
                    var _whole_l = math.norm(_dp);

                    // We only want to go this.g.gamma[vi]/|dp| of the way:
                    var _mid_dp = math.multiply(this.g.gamma[_vi] / _whole_l, _dp);
                    var _mp = math.add(_vp, _mid_dp);

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

            var vp = pts[vi];
            var ep = pts[ei];

            // We really want to push 'center' points as anchor points and
            // real points as control points
            var dp = math.subtract(last_ep, vp);

            // Find center point by linear interpolation; whole way across is
            // |dp| \approx this.g.gamma[vi]+this.g.gamma[ei]
            var whole_l = math.norm(dp);

            // We only want to go this.g.gamma[vi]/|dp| of the way:
            var mid_dp = math.multiply(this.g.gamma[vi] / whole_l, dp);
            var mp = math.add(vp, mid_dp);

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
    }]);

    return MeshDraw;
}();

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
    }

    _createClass(Map4v, [{
        key: 'triangulate',
        value: function triangulate() {
            var _this = this;

            var bdry_face_i = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : undefined;

            var faces = this.faces();

            if (bdry_face_i === undefined) {
                var sizes = faces.map(function (f) {
                    return f.length;
                }, this);
                bdry_face_i = sizes.indexOf(math.max(sizes));
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
            var verts = math.range(0, this.nv + this.ne + faces.length - 1).toArray();

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
                var _iteratorNormalCompletion5 = true;
                var _didIteratorError5 = false;
                var _iteratorError5 = undefined;

                try {
                    for (var _iterator5 = face[Symbol.iterator](), _step5; !(_iteratorNormalCompletion5 = (_step5 = _iterator5.next()).done); _iteratorNormalCompletion5 = true) {
                        var arc = _step5.value;

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
            }

            return [verts, bdy_verts, edges, triangles];
        }
    }, {
        key: 'faces',
        value: function faces() {
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

            return faces;
        }
    }, {
        key: 'component',
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
        key: 'new_arc',
        value: function new_arc(idx) {
            this.arcs[idx] = { index: idx, edge: undefined, edgepos: undefined, vert: undefined, vertpos: undefined };
        }
    }, {
        key: 'set_edge',
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
        key: 'set_vert',
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
        key: 'adjacent_edges',
        value: function adjacent_edges(vert) {
            var adj = [];
            var _iteratorNormalCompletion6 = true;
            var _didIteratorError6 = false;
            var _iteratorError6 = undefined;

            try {
                for (var _iterator6 = this.edges[Symbol.iterator](), _step6; !(_iteratorNormalCompletion6 = (_step6 = _iterator6.next()).done); _iteratorNormalCompletion6 = true) {
                    var edge = _step6.value;

                    if (edge.indexOf(vert) >= 0) {
                        adj.push(edge);
                    }
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

            return adj;
        }
    }, {
        key: 'adjacent_faces',
        value: function adjacent_faces(vert) {
            return this.faces.filter(function (f) {
                return f.indexOf(vert) >= 0;
            });
        }
    }, {
        key: 'valence',
        value: function valence(vert) {
            return this.adjacent_faces(vert).length;
        }
    }, {
        key: 'valences',
        value: function valences() {
            return this.verts.map(this.valence, this);
        }
    }, {
        key: 'min_valence',
        value: function min_valence() {
            return math.min(this.valences());
        }
    }, {
        key: 'chi',
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

        this.gamma = radius_map;
        this.u = this.conf_factor(radius_map);
        this.l = math.zeros(this.n, this.n, 'sparse');

        this.theta = {};
        this.phi = edge_weights;

        this.update();
    }

    _createClass(DiscreteRiemannMetric, [{
        key: 'conf_factor',
        value: function conf_factor(gamma) {
            return math.log(gamma);
        }
    }, {
        key: 'compute_length',
        value: function compute_length(edge) {
            var _math;

            var g = math.subset(this.gamma, math.index(edge));
            return math.sqrt(2 * g[0] * g[1] * math.cos(this.phi.subset((_math = math).index.apply(_math, _toConsumableArray(edge)))) + Math.pow(g[0], 2) + Math.pow(g[1], 2));
        }
    }, {
        key: 'length',
        value: function length(edge) {
            var _math2;

            return this.l.subset((_math2 = math).index.apply(_math2, _toConsumableArray(edge)));
        }
    }, {
        key: 'abc_for_vert',
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
        key: 'compute_angle',
        value: function compute_angle(face, vert) {
            var abc = this.abc_for_vert(face, vert);

            var ratio = (Math.pow(abc[0], 2) + Math.pow(abc[1], 2) - Math.pow(abc[2], 2)) / (2.0 * abc[0] * abc[1]);
            return Math.acos(ratio);
        }
    }, {
        key: 'angle',
        value: function angle(face, vert) {
            var other_vs = face.filter(function (v) {
                return v != vert;
            });
            return this.theta[[vert, other_vs[0], other_vs[1]]];
        }
    }, {
        key: 'curvature',
        value: function curvature(vert) {
            var _this4 = this;

            var adj_faces = this.mesh.adjacent_faces(vert);
            if (this.mesh.bdyverts.indexOf(vert) >= 0) {
                return Math.PI - math.sum(adj_faces.map(function (face) {
                    return _this4.angle(face, vert);
                }, this));
            } else {
                return 2 * Math.PI - math.sum(adj_faces.map(function (face) {
                    return _this4.angle(face, vert);
                }, this));
            }
        }
    }, {
        key: 'update',
        value: function update() {
            this.gamma = math.exp(this.u);

            var _iteratorNormalCompletion7 = true;
            var _didIteratorError7 = false;
            var _iteratorError7 = undefined;

            try {
                for (var _iterator7 = this.mesh.edges[Symbol.iterator](), _step7; !(_iteratorNormalCompletion7 = (_step7 = _iterator7.next()).done); _iteratorNormalCompletion7 = true) {
                    var edge = _step7.value;

                    var l = this.compute_length(edge);
                    this.l.subset(math.index(edge[0], edge[1]), l);
                    this.l.subset(math.index(edge[1], edge[0]), l);
                }

                // Set angles using law of cosines
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

            var _iteratorNormalCompletion8 = true;
            var _didIteratorError8 = false;
            var _iteratorError8 = undefined;

            try {
                for (var _iterator8 = this.mesh.faces[Symbol.iterator](), _step8; !(_iteratorNormalCompletion8 = (_step8 = _iterator8.next()).done); _iteratorNormalCompletion8 = true) {
                    var face = _step8.value;
                    var _iteratorNormalCompletion9 = true;
                    var _didIteratorError9 = false;
                    var _iteratorError9 = undefined;

                    try {
                        for (var _iterator9 = partition_face(face)[Symbol.iterator](), _step9; !(_iteratorNormalCompletion9 = (_step9 = _iterator9.next()).done); _iteratorNormalCompletion9 = true) {
                            var part = _step9.value;

                            var theta = this.compute_angle(face, part[0]);
                            this.theta[[part[0], part[1][0], part[1][1]]] = theta;
                            this.theta[[part[0], part[1][1], part[1][0]]] = theta;
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
                }

                // Set curvatures
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

            this.K = math.sparse(this.mesh.verts.map(this.curvature, this));
        }
    }, {
        key: 'newton',
        value: function newton() {
            var target_K = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : null;
            var dt = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0.05;
            var thresh = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 1e-4;

            if (target_K == null) {}

            var g = new DiscreteRiemannMetric(this.mesh, this.gamma, this.phi);

            var K = g.K;
            var DeltaK = math.subtract(target_K, K);

            var _failsafe = 0;
            while (math.max(math.abs(DeltaK)) > thresh) {
                var H = this.hessian();
                var deltau = least_squares(H, DeltaK);

                g.u = math.subtract(g.u, math.transpose(math.multiply(dt, deltau)).toArray()[0]);

                g.update();

                K = g.K;
                DeltaK = math.subtract(target_K, K);

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
        key: 'tau2',
        value: function tau2(l_jk, g_j, g_k) {
            return .5 * (Math.pow(l_jk, 2) + Math.pow(g_j, 2) - Math.pow(g_k, 2));
        }
    }, {
        key: 'face_area',
        value: function face_area(face) {
            var gamma = this.theta[face];
            var a = this.length([face[0], face[1]]);
            var b = this.length([face[0], face[2]]);

            return .5 * a * b * math.sin(gamma);
        }
    }, {
        key: 'hessian',
        value: function hessian() {
            var n = this.mesh.verts.length;
            var H = math.zeros(n, n, 'sparse');
            var t = this.tau2;

            var _iteratorNormalCompletion10 = true;
            var _didIteratorError10 = false;
            var _iteratorError10 = undefined;

            try {
                for (var _iterator10 = this.mesh.faces[Symbol.iterator](), _step10; !(_iteratorNormalCompletion10 = (_step10 = _iterator10.next()).done); _iteratorNormalCompletion10 = true) {
                    var face = _step10.value;

                    var i = face[0];
                    var j = face[1];
                    var k = face[2];

                    var l_k = math.subset(this.l, math.index(i, j));
                    var l_i = math.subset(this.l, math.index(j, k));
                    var l_j = math.subset(this.l, math.index(k, i));

                    var g_i = math.subset(this.gamma, math.index(i));
                    var g_j = math.subset(this.gamma, math.index(j));
                    var g_k = math.subset(this.gamma, math.index(k));

                    var th_i = this.angle(face, i);
                    var th_j = this.angle(face, j);
                    var th_k = this.angle(face, k);

                    var A = this.face_area(face);

                    var L = math.diag([l_i, l_j, l_k]);
                    var D = math.matrix([[0, t(l_i, g_j, g_k), t(l_i, g_k, g_j)], [t(l_j, g_i, g_k), 0, t(l_j, g_k, g_i)], [t(l_k, g_i, g_j), t(l_k, g_j, g_i), 0]]);

                    var Theta = math.cos(math.matrix([[Math.PI, th_k, th_j], [th_k, Math.PI, th_i], [th_j, th_i, Math.PI]]));

                    var Tijk = math.multiply(-.5 / A, math.multiply(math.multiply(math.multiply(L, Theta), math.inv(L)), D));

                    var _arr = [0, 1, 2];
                    for (var _i = 0; _i < _arr.length; _i++) {
                        var rowi = _arr[_i];
                        var a = [i, j, k][rowi];
                        var _arr2 = [0, 1, 2];
                        for (var _i2 = 0; _i2 < _arr2.length; _i2++) {
                            var coli = _arr2[_i2];
                            var b = [i, j, k][coli];
                            H.subset(math.index(a, b), math.subset(H, math.index(a, b)) + math.subset(Tijk, math.index(rowi, coli)));
                        }
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

            return H;
        }
    }], [{
        key: 'from_triangle_mesh',
        value: function from_triangle_mesh(mesh) {
            var n = mesh.verts.length;
            var gamma = mesh.verts.map(function () {
                return 1;
            }, this);
            var phi = math.zeros(n, n, 'sparse');

            var _iteratorNormalCompletion11 = true;
            var _didIteratorError11 = false;
            var _iteratorError11 = undefined;

            try {
                for (var _iterator11 = mesh.edges[Symbol.iterator](), _step11; !(_iteratorNormalCompletion11 = (_step11 = _iterator11.next()).done); _iteratorNormalCompletion11 = true) {
                    var edge = _step11.value;

                    phi.subset(math.index(edge[0], edge[1]), 0);
                    phi.subset(math.index(edge[1], edge[0]), 0);
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

            return new DiscreteRiemannMetric(mesh, gamma, phi);
        }
    }]);

    return DiscreteRiemannMetric;
}();

function adjacent_faces_and_edge(faces, face) {
    var adj = [];

    var _iteratorNormalCompletion12 = true;
    var _didIteratorError12 = false;
    var _iteratorError12 = undefined;

    try {
        for (var _iterator12 = faces[Symbol.iterator](), _step12; !(_iteratorNormalCompletion12 = (_step12 = _iterator12.next()).done); _iteratorNormalCompletion12 = true) {
            var aface = _step12.value;

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

    return adj;
}

function adj_or_faces_and_edge(faces, face) {
    var adj = [];

    var _iteratorNormalCompletion13 = true;
    var _didIteratorError13 = false;
    var _iteratorError13 = undefined;

    try {
        for (var _iterator13 = faces[Symbol.iterator](), _step13; !(_iteratorNormalCompletion13 = (_step13 = _iterator13.next()).done); _iteratorNormalCompletion13 = true) {
            var aface = _step13.value;

            if (aface.indexOf(face[0]) >= 0 && aface.indexOf(face[1]) >= 0) {
                adj.push([aface, [face[0], face[1]]]);
            } else if (aface.indexOf(face[1]) >= 0 && aface.indexOf(face[2]) >= 0) {
                adj.push([aface, [face[1], face[2]]]);
            } else if (aface.indexOf(face[2]) >= 0 && aface.indexOf(face[0]) >= 0) {
                adj.push([aface, [face[2], face[0]]]);
            }
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

    var _iteratorNormalCompletion14 = true;
    var _didIteratorError14 = false;
    var _iteratorError14 = undefined;

    try {
        for (var _iterator14 = adjacent_faces_and_edge(to_orient, f_0)[Symbol.iterator](), _step14; !(_iteratorNormalCompletion14 = (_step14 = _iterator14.next()).done); _iteratorNormalCompletion14 = true) {
            var adj_pair = _step14.value;

            adj_queue.add(adj_pair);
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

    var _failsafe = 0;
    while (to_orient.length > 0) {
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
        var _iteratorNormalCompletion15 = true;
        var _didIteratorError15 = false;
        var _iteratorError15 = undefined;

        try {
            for (var _iterator15 = adjacent_faces_and_edge(to_orient, [i, j, k])[Symbol.iterator](), _step15; !(_iteratorNormalCompletion15 = (_step15 = _iterator15.next()).done); _iteratorNormalCompletion15 = true) {
                var _adj_pair2 = _step15.value;

                adj_queue.add(_adj_pair2);
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
    }

    return oriented;
}

function u_theta(theta) {
    return [math.cos(theta), math.sin(theta)];
}

function embed_faces(g) {
    var mesh = g.mesh;
    var x = math.multiply(-30, math.ones(g.n, 2));
    var phi = {};
    var pi = Math.PI;

    var faces = orient_faces(mesh.faces);
    var to_embed = faces.slice();
    var embed_queue = new Set();

    var f_0 = to_embed.pop();
    var i = f_0[0],
        j = f_0[1],
        k = f_0[2];

    x.subset(math.index(i, [0, 1]), [0, 0]);
    x.subset(math.index(j, [0, 1]), [g.length([i, j]), 0]);

    var phi_ik = g.angle(f_0, i) % (2 * pi);

    x.subset(math.index(k, [0, 1]), math.multiply(g.length([i, k]), u_theta(phi_ik)));

    var phi_jk = (pi - g.angle(f_0, j)) % (2 * pi);

    phi[[i, j]] = 0;
    phi[[j, i]] = pi;
    phi[[j, k]] = phi_jk;
    phi[[k, j]] = pi + phi_jk;
    phi[[i, k]] = phi_ik;
    phi[[k, i]] = pi + phi_ik;

    var _iteratorNormalCompletion16 = true;
    var _didIteratorError16 = false;
    var _iteratorError16 = undefined;

    try {
        for (var _iterator16 = adj_or_faces_and_edge(to_embed, f_0)[Symbol.iterator](), _step16; !(_iteratorNormalCompletion16 = (_step16 = _iterator16.next()).done); _iteratorNormalCompletion16 = true) {
            var adj_pair = _step16.value;

            embed_queue.add(adj_pair);
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
            x.subset(math.index(_k, [0, 1]), math.add(x.subset(math.index(_i3, [0, 1])).toArray()[0], math.multiply(g.length([_i3, _k]), u_theta(phi_ik))));
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
        var _iteratorNormalCompletion17 = true;
        var _didIteratorError17 = false;
        var _iteratorError17 = undefined;

        try {
            for (var _iterator17 = adj_or_faces_and_edge(to_embed, F)[Symbol.iterator](), _step17; !(_iteratorNormalCompletion17 = (_step17 = _iterator17.next()).done); _iteratorNormalCompletion17 = true) {
                var _adj_pair4 = _step17.value;

                embed_queue.add(_adj_pair4);
            }
        } catch (err) {
            _didIteratorError17 = true;
            _iteratorError17 = err;
        } finally {
            try {
                if (!_iteratorNormalCompletion17 && _iterator17.return) {
                    _iterator17.return();
                }
            } finally {
                if (_didIteratorError17) {
                    throw _iteratorError17;
                }
            }
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

var meshDraw = new MeshDraw("#knot-draw");

function drawMap(sigma) {
    var cross_bend = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 8;

    var tstart = Date.now();
    meshDraw.clear();

    var m4v = new Map4v(sigma);

    var triangulation = m4v.triangulate();

    //console.log(triangulation);

    var testMesh = new TriangleMesh(triangulation[0], triangulation[1], triangulation[2], triangulation[3]);

    var cpmetric = DiscreteRiemannMetric.from_triangle_mesh(testMesh);

    var K = math.zeros(testMesh.verts.length, 1);
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
    K.subset(math.index(bdyCross, 0), bdyCross.map(function (bci) {
        return -Math.PI / fac;
    }));
    K.subset(math.index(bdyEdge, 0), bdyEdge.map(function (bei) {
        return (2 * Math.PI + bdyCross.length * Math.PI / fac) / bdyEdge.length;
    }));

    var flat_poly = cpmetric.newton(K, 1, 5e-2);

    var embedding = embed_faces(flat_poly);

    meshDraw.set_embedding(flat_poly, embedding, m4v);
    console.log("Computation finished in: " + (Date.now() - tstart) + " milliseconds");
}

// Trefoil
//let tref = new Map4v([[0, 6, 11, 5], [1, 8, 2, 7], [3, 9, 4, 10]]);
// 5_1
//let tref = new Map4v([[0, 10, 19, 9], [1, 12, 2, 11], [3, 13, 4, 14], [5, 16, 6, 15], [7, 17, 8, 18]]);
// Random 10x
//let tref = new Map4v([[1, 0, 2, 39], [3, 34, 4, 33], [5, 28, 6, 27], [7, 29, 8, 30], [9, 12, 10, 11], [13, 20, 14, 19], [15, 17, 16, 18], [21, 35, 22, 36], [23, 38, 24, 37], [25, 32, 26, 31]]);
//[[0, 6, 23, 5], [1, 16, 2, 15], [3, 17, 4, 18], [7, 14, 8, 13], [9, 19, 10, 20], [11, 22, 12, 21]]
// [[0, 25, 35, 26], [1, 27, 2, 28], [3, 18, 4, 17], [5, 15, 6, 16], [7, 30, 8, 29], [9, 23, 10, 24], [11, 33, 12, 34], [13, 20, 14, 19], [21, 32, 22, 31]]
//[[1, 52, 2, 51], [0, 54, 59, 53], [3, 10, 4, 9], [5, 47, 6, 48], [7, 50, 8, 49], [11, 33, 12, 34], [13, 32, 14, 31], [15, 58, 16, 57], [17, 55, 18, 56], [19, 46, 20, 45], [21, 39, 22, 40], [23, 42, 24, 41], [25, 43, 26, 44], [27, 38, 28, 37], [29, 35, 30, 36]]
//[[0, 6, 59, 5], [1, 20, 2, 19], [3, 21, 4, 22], [7, 38, 8, 37], [9, 35, 10, 36], [11, 41, 12, 42], [13, 44, 14, 43], [15, 45, 16, 46], [17, 32, 18, 31], [23, 58, 24, 57], [25, 52, 26, 51], [27, 49, 28, 50], [29, 55, 30, 56], [33, 40, 34, 39], [47, 54, 48, 53]]
//let sigma = [[0, 41, 99, 42], [1, 47, 2, 48], [3, 34, 4, 33], [5, 96, 6, 95], [7, 13, 8, 14], [9, 68, 10, 67], [11, 69, 12, 70], [15, 62, 16, 61], [17, 59, 18, 60], [19, 85, 20, 86], [21, 27, 22, 28], [23, 82, 24, 81], [25, 83, 26, 84], [29, 88, 30, 87], [31, 93, 32, 94], [35, 46, 36, 45], [37, 43, 38, 44], [39, 98, 40, 97], [49, 76, 50, 75], [51, 73, 52, 74], [53, 72, 54, 71], [55, 65, 56, 66], [57, 64, 58, 63], [77, 92, 78, 91], [79, 89, 80, 90]];
// [[1, 31, 2, 32], [0, 126, 199, 125], [3, 102, 4, 101], [5, 12, 6, 11], [7, 114, 8, 113], [9, 111, 10, 112], [13, 116, 14, 115], [15, 117, 16, 118], [17, 103, 18, 104], [19, 30, 20, 29], [21, 124, 22, 123], [23, 53, 24, 54], [25, 56, 26, 55], [27, 121, 28, 122], [33, 100, 34, 99], [35, 130, 36, 129], [37, 131, 38, 132], [39, 134, 40, 133], [41, 151, 42, 152], [43, 86, 44, 85], [45, 175, 46, 176], [47, 166, 48, 165], [49, 75, 50, 76], [51, 93, 52, 94], [57, 68, 58, 67], [59, 142, 60, 141], [61, 143, 62, 144], [63, 146, 64, 145], [65, 139, 66, 140], [69, 92, 70, 91], [71, 170, 72, 169], [73, 171, 74, 172], [77, 196, 78, 195], [79, 153, 80, 154], [81, 191, 82, 192], [83, 190, 84, 189], [87, 150, 88, 149], [89, 147, 90, 148], [95, 198, 96, 197], [97, 127, 98, 128], [105, 120, 106, 119], [107, 138, 108, 137], [109, 135, 110, 136], [155, 193, 156, 194], [157, 188, 158, 187], [159, 182, 160, 181], [161, 183, 162, 184], [163, 177, 164, 178], [167, 174, 168, 173], [179, 186, 180, 185]]
var sigma = [[1, 287, 2, 288], [0, 289, 299, 290], [3, 126, 4, 125], [5, 16, 6, 15], [7, 13, 8, 14], [9, 255, 10, 256], [11, 254, 12, 253], [17, 127, 18, 128], [19, 282, 20, 281], [21, 132, 22, 131], [23, 277, 24, 278], [25, 235, 26, 236], [27, 222, 28, 221], [29, 223, 30, 224], [31, 38, 32, 37], [33, 228, 34, 227], [35, 225, 36, 226], [39, 230, 40, 229], [41, 156, 42, 155], [43, 157, 44, 158], [45, 164, 46, 163], [47, 165, 48, 166], [49, 56, 50, 55], [51, 198, 52, 197], [53, 195, 54, 196], [57, 243, 58, 244], [59, 242, 60, 241], [61, 231, 62, 232], [63, 234, 64, 233], [65, 248, 66, 247], [67, 249, 68, 250], [69, 276, 70, 275], [71, 133, 72, 134], [73, 83, 74, 84], [75, 82, 76, 81], [77, 283, 78, 284], [79, 286, 80, 285], [85, 295, 86, 296], [87, 294, 88, 293], [89, 291, 90, 292], [91, 298, 92, 297], [93, 136, 94, 135], [95, 109, 96, 110], [97, 108, 98, 107], [99, 122, 100, 121], [101, 259, 102, 260], [103, 262, 104, 261], [105, 119, 106, 120], [111, 273, 112, 274], [113, 272, 114, 271], [115, 265, 116, 266], [117, 264, 118, 263], [123, 137, 124, 138], [129, 280, 130, 279], [139, 257, 140, 258], [141, 268, 142, 267], [143, 269, 144, 270], [145, 252, 146, 251], [147, 237, 148, 238], [149, 216, 150, 215], [151, 217, 152, 218], [153, 220, 154, 219], [159, 210, 160, 209], [161, 211, 162, 212], [167, 190, 168, 189], [169, 187, 170, 188], [171, 205, 172, 206], [173, 200, 174, 199], [175, 201, 176, 202], [177, 204, 178, 203], [179, 186, 180, 185], [181, 191, 182, 192], [183, 194, 184, 193], [207, 214, 208, 213], [239, 245, 240, 246]];
//[[0, 21, 27, 22], [1, 23, 2, 24], [3, 26, 4, 25], [5, 20, 6, 19], [7, 14, 8, 13], [9, 15, 10, 16], [11, 18, 12, 17]]);

drawMap(sigma);

document.getElementById("map_submit").onclick = function (ev) {
    document.querySelector(ev.target.dataset.errbox).textContent = "";
    try {
        drawMap(JSON.parse(document.querySelector(ev.target.dataset.target).value));
    } catch (err) {
        document.querySelector(ev.target.dataset.errbox).textContent = "Error";
        console.log("Error:", err);
    }
};