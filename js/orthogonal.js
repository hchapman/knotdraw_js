'use strict';

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

//let lab = new Lalolab();

function least_squares(X /* : Matrix */, Y /* : Matrix */) /* : least_squares */{
    var betaHat = solve(mul(X, transpose(X)), mul(X, Y));

    return betaHat;
}

var MeshEdge = function () {
    function MeshEdge(parent, start, stop) {
        _classCallCheck(this, MeshEdge);

        this.parent = parent;
        this.start = start;
        this.stop = stop;

        this.svg = this.parent.edgeG.line(start.x(), start.y(), stop.x(), stop.y());
        this.svg.addClass("scaffold");
    }

    _createClass(MeshEdge, [{
        key: 'set_nodes',
        value: function set_nodes(start, end) {
            var anim_ms = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0;

            if (anim_ms == 0) {
                this.svg.attr({ 'x1': start.x(), 'y1': start.y(), 'x2': end.x(), 'y2': end.y() });
            } else {
                this.svg.animate({ 'x1': start.x(), 'y1': start.y(), 'x2': end.x(), 'y2': end.y() }, anim_ms);
            }
        }
    }]);

    return MeshEdge;
}();

var MeshNode = function () {
    function MeshNode(parent, x, y) {
        _classCallCheck(this, MeshNode);

        this.parent = parent;
        this.svg = this.parent.nodeG.circle(x, y, .3);
        this._x = x;
        this._y = y;
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
        key: 'move',
        value: function move(x, y) {
            var new_r = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : undefined;
            var anim_ms = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 0;

            this._x = x;
            this._y = y;
            if (anim_ms == 0) {
                // Do not animate
                this.svg.attr({ 'cx': x, 'cy': y, 'r': new_r });
            } else {
                // Animate the motion
                this.svg.animate({ 'cx': x, 'cy': y, 'r': new_r }, anim_ms);
            }
        }
    }, {
        key: 'set_r',
        value: function set_r(r) {
            var anim_ms = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;

            if (anim_ms == 0) {
                this.svg.attr({ 'r': r });
            } else {
                this.svg.animate({ 'r': r }, anim_ms);
            }
        }
    }, {
        key: 'cur_x',
        value: function cur_x() {
            return this.svg.attr('cx');
        }
    }, {
        key: 'x',
        value: function x() {
            return this._x;
        }
    }, {
        key: 'cur_y',
        value: function cur_y() {
            return this.svg.attr('cy');
        }
    }, {
        key: 'y',
        value: function y() {
            return this._y;
        }
    }, {
        key: 'onClick',
        value: function onClick(e) {
            //console.log(this);
            //console.log(this.obj);
            if (this.obj.length == 2) {
                this.parent.delete_face(this.i);
            }
        }
    }]);

    return MeshNode;
}();

var CompEdgeNode = function () {
    function CompEdgeNode(parent, x, y, r) {
        _classCallCheck(this, CompEdgeNode);

        this.parent = parent;
        this.svg = this.parent.knotG.circle(x, y, r);
        this._x = x;
        this._y = y;
        this._r = r;

        this.svg.node.addEventListener('click', this.onClick.bind(this));

        this.dragging = false;
    }

    _createClass(CompEdgeNode, [{
        key: 'set_obj',
        value: function set_obj(i, obj) {
            this.i = i;
            this.obj = obj;
        }
    }, {
        key: 'move',
        value: function move(x, y, r) {
            var anim_ms = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 0;

            this._x = x;
            this._y = y;
            this._r = r;
            if (anim_ms == 0) {
                // Do not animate
                this.svg.attr({ 'cx': x, 'cy': y, 'r': r });
            } else {
                // Animate the motion
                this.svg.animate({ 'cx': x, 'cy': y, 'r': r }, anim_ms);
            }
        }
    }, {
        key: 'set_r',
        value: function set_r(r) {
            var anim_ms = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;

            if (anim_ms == 0) {
                this.svg.attr({ 'r': r });
            } else {
                this.svg.animate({ 'r': r }, anim_ms);
            }
        }
    }, {
        key: 'cur_x',
        value: function cur_x() {
            return this.svg.attr('cx');
        }
    }, {
        key: 'x',
        value: function x() {
            return this._x;
        }
    }, {
        key: 'cur_y',
        value: function cur_y() {
            return this.svg.attr('cy');
        }
    }, {
        key: 'y',
        value: function y() {
            return this._y;
        }
    }, {
        key: 'onClick',
        value: function onClick(e) {
            //console.log(this);
            //console.log(this.obj);
            if (this.obj.length == 2) {
                this.parent.delete_face(this.i);
            }
        }
    }]);

    return CompEdgeNode;
}();

var MeshDraw = function () {
    function MeshDraw(div) {
        _classCallCheck(this, MeshDraw);

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

    _createClass(MeshDraw, [{
        key: 'clear',
        value: function clear() {
            this.nodes = [];
            this.edges = {};
            this.comps = [];

            this.edgeG.clear();
            this.nodeG.clear();
            this.knotG.clear();
        }
    }, {
        key: 'set_link_diagram',
        value: function set_link_diagram(ld) {
            var pts = array2mat(ld.verts);
            //console.log(pts);
            var min_x = min(get(pts, range(), 0));
            var min_y = min(get(pts, range(), 1));

            var max_x = max(get(pts, range(), 0));
            var max_y = max(get(pts, range(), 1));

            var wid = max_x - min_x;
            var hgt = max_y - min_y;

            var dx = wid * 0.05;
            var dy = hgt * 0.05;

            //console.log(min_x, max_x);

            //this.nodeG.clear();
            //this.edgeG.clear();

            this.draw.attr({ viewBox: [min_x - dx, min_y - dy, wid + 2 * dx, hgt + 2 * dy].join(",") });

            var i = 0;
            //console.log(ld.verts, "!!");
            var _iteratorNormalCompletion = true;
            var _didIteratorError = false;
            var _iteratorError = undefined;

            try {
                for (var _iterator = ld.verts[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
                    var vert = _step.value;

                    //this.nodeG.circle(vert[0], vert[1], .25);
                    if (this.nodes[i] === undefined) {
                        var t = this.nodeG.text(vert[0], vert[1], i.toString());
                        t.attr({ "style": "font-size: .5px;" });
                        this.nodes[i] = t;
                    } else {
                        this.nodes[i].attr({ 'x': vert[0], 'y': vert[1] });
                    }
                    i++;
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
                for (var _iterator2 = ld.edges[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                    var edge = _step2.value;
                    var _ref = [ld.verts[edge[0]], ld.verts[edge[1]]],
                        a = _ref[0],
                        b = _ref[1];

                    if (this.edges[edge] === undefined) {
                        this.edges[edge] = this.edgeG.line(a[0], a[1], b[0], b[1]).addClass("edge");
                    } else {
                        this.edges[edge].attr({ 'x1': a[0], 'y1': a[1], 'x2': b[0], 'y2': b[1] });
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

            var _iteratorNormalCompletion3 = true;
            var _didIteratorError3 = false;
            var _iteratorError3 = undefined;

            try {
                for (var _iterator3 = conv.comps[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
                    var comp = _step3.value;

                    for (var pi = 0; pi < comp.length; pi++) {
                        var edge = [comp[pi], comp[(pi + 1) % comp.length]];
                        //console.log(edge);
                        if (edge[0] in this.nodes && edge[1] in this.nodes) {
                            var mesh_edge = this.add_edge(edge[0], edge[1]);
                            mesh_edge.svg.addClass('edge');
                        }
                    }
                }

                //console.log(this.map4v.faces);
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

            for (var idx in conv.faces) {
                //console.log(fi);
                if (conv.faces[idx].length > 0) {
                    var mesh_face = this.add_face(conv.faces[idx][0], parseInt(idx));
                }
            }

            var ci = 0;
            var _iteratorNormalCompletion4 = true;
            var _didIteratorError4 = false;
            var _iteratorError4 = undefined;

            try {
                for (var _iterator4 = conv.comps[Symbol.iterator](), _step4; !(_iteratorNormalCompletion4 = (_step4 = _iterator4.next()).done); _iteratorNormalCompletion4 = true) {
                    var component = _step4.value;

                    this.comps[ci] = this.add_component(component, map4v, points, conv);
                    this.comps[ci].addClass('q' + ci + "-9");
                    ci += 1;
                }
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
        }
    }, {
        key: 'update_embedding',
        value: function update_embedding(g /*metric*/, embedding, map4v /*original map*/, conv) {
            var _this = this;

            var anim_ms = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : 0;

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

            Snap.animate(this.draw.attr("viewBox").vb.split(" "), [min_x - dx, min_y - dy, wid + 2 * dx, hgt + 2 * dy], function (v) {
                _this.draw.attr("viewBox", v.join(" "));
            }, this.anim_ms);

            //console.log(points);
            for (var i = 0; i < points.m; i++) {
                var node = this.nodes[i];
                node.move(points.val[i * points.n], points.val[i * points.n + 1], this.g.gamma[i], this.anim_ms);
            }

            var _iteratorNormalCompletion5 = true;
            var _didIteratorError5 = false;
            var _iteratorError5 = undefined;

            try {
                for (var _iterator5 = conv.comps[Symbol.iterator](), _step5; !(_iteratorNormalCompletion5 = (_step5 = _iterator5.next()).done); _iteratorNormalCompletion5 = true) {
                    var comp = _step5.value;

                    for (var pi = 0; pi < comp.length; pi++) {
                        var edge = [comp[pi], comp[(pi + 1) % comp.length]];
                        //console.log(edge);
                        if (edge[0] in this.nodes && edge[1] in this.nodes) {
                            this.update_edge(edge[0], edge[1]);
                        }
                    }
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

            var ci = 0;
            var _iteratorNormalCompletion6 = true;
            var _didIteratorError6 = false;
            var _iteratorError6 = undefined;

            try {
                for (var _iterator6 = conv.comps[Symbol.iterator](), _step6; !(_iteratorNormalCompletion6 = (_step6 = _iterator6.next()).done); _iteratorNormalCompletion6 = true) {
                    var component = _step6.value;

                    this.update_component(ci, component, map4v, points, conv);
                    ci += 1;
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
        key: 'update_edge',
        value: function update_edge(i, j) {
            this.edges[[i, j]].set_nodes(this.nodes[i], this.nodes[j], this.anim_ms);
        }
    }, {
        key: 'add_face',
        value: function add_face(i, fi) {
            var face_node = this.nodes[i];
            face_node.svg.addClass("face");
            face_node.set_obj(i, this.map4v.faces[fi]);
        }
    }, {
        key: 'component_path_gen',
        value: function component_path_gen(component, map4v, points) {
            var path = [];
            var pts = points;

            var root_arc = component[0];

            var last_ep = undefined;

            for (var i = 0; i < component.length; i++) {
                var j = (i + 1) % component.length;

                var p_i = component[i];
                var p_j = component[j];
                var p_a = [pts.val[p_i * pts.n], pts.val[p_i * pts.n + 1]];
                var p_b = [pts.val[p_j * pts.n], pts.val[p_j * pts.n + 1]];

                var dp = sub(p_b, p_a);
                var whole_l = norm(dp);
                var mid_dp = mul(this.g.gamma[p_i] / whole_l, dp);
                var p_mid = add(p_a, mid_dp);

                path.push(p_mid); // midpoint is anchor
                path.push(p_b); // p_b is quadratic bezier control
            }

            // Close the path
            path.push(path[0]);

            //console.log(path);

            var pathStr = "M";
            pathStr += path[0].join(",");

            var idx = 0;
            var _iteratorNormalCompletion7 = true;
            var _didIteratorError7 = false;
            var _iteratorError7 = undefined;

            try {
                for (var _iterator7 = path.slice(1)[Symbol.iterator](), _step7; !(_iteratorNormalCompletion7 = (_step7 = _iterator7.next()).done); _iteratorNormalCompletion7 = true) {
                    var pt = _step7.value;

                    if (idx % 2 == 0) {
                        pathStr += "Q";
                    } else {
                        pathStr += " ";
                    }
                    idx += 1;
                    pathStr += pt.join(",");
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

            return [path, pathStr];
        }
    }, {
        key: 'quadratic_segments',
        value: function quadratic_segments(path) {
            var slices = [];
            for (var i = 0; i < path.length - 4; i += 2) {
                slices.push(path.slice(i, i + 3));
            }
            slices.unshift(path.slice(path.length - 3, path.length));

            return slices;
        }
    }, {
        key: 'quadratic_segments_to_cubic',
        value: function quadratic_segments_to_cubic(segs) {
            var csegs = [];
            var _iteratorNormalCompletion8 = true;
            var _didIteratorError8 = false;
            var _iteratorError8 = undefined;

            try {
                for (var _iterator8 = segs[Symbol.iterator](), _step8; !(_iteratorNormalCompletion8 = (_step8 = _iterator8.next()).done); _iteratorNormalCompletion8 = true) {
                    var seg = _step8.value;

                    var _seg = _slicedToArray(seg, 3),
                        q0 = _seg[0],
                        q1 = _seg[1],
                        q2 = _seg[2];

                    csegs.push([q0, add(q0, mul(2 / 3, sub(q1, q0))), add(q2, mul(2 / 3, sub(q1, q2))), q2]);
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

            return csegs;
        }
    }, {
        key: 'add_component',
        value: function add_component(component, map4v, points, conv) {
            var _this2 = this;

            // A path of the form anchor, control, anchor...
            var _component_path_gen = this.component_path_gen(component, map4v, points),
                _component_path_gen2 = _slicedToArray(_component_path_gen, 2),
                path = _component_path_gen2[0],
                pathStr = _component_path_gen2[1];

            var segs = this.quadratic_segments_to_cubic(this.quadratic_segments(path));

            var _loop = function _loop(si) {
                var tri_i = component[si];
                if (!conv.edges.some(function (e) {
                    return e.includes(tri_i);
                })) {
                    return 'continue';
                }

                var seg = segs[si];

                var p = Snap.path.findDotsAtSegment(seg[0][0], seg[0][1], seg[1][0], seg[1][1], seg[2][0], seg[2][1], seg[3][0], seg[3][1], .5);
                _this2.comp_edgenodes[tri_i] = new CompEdgeNode(_this2, p.x, p.y, _this2.g.gamma[tri_i] / 2);
            };

            for (var si = 0; si < segs.length; si++) {
                var _ret = _loop(si);

                if (_ret === 'continue') continue;
            }

            //console.log(pathStr);
            var knot = this.knotG.path(pathStr + "Z");
            knot.addClass('knot');

            return knot;
        }
    }, {
        key: 'update_component',
        value: function update_component(ci, component, map4v, points, conv) {
            var _this3 = this;

            // A path of the form anchor, control, anchor...
            var _component_path_gen3 = this.component_path_gen(component, map4v, points),
                _component_path_gen4 = _slicedToArray(_component_path_gen3, 2),
                path = _component_path_gen4[0],
                pathStr = _component_path_gen4[1];

            var segs = this.quadratic_segments_to_cubic(this.quadratic_segments(path));

            var _loop2 = function _loop2(si) {
                var tri_i = component[si];
                if (!conv.edges.some(function (e) {
                    return e.includes(tri_i);
                })) {
                    return 'continue';
                }

                var seg = segs[si];

                var p = Snap.path.findDotsAtSegment(seg[0][0], seg[0][1], seg[1][0], seg[1][1], seg[2][0], seg[2][1], seg[3][0], seg[3][1], .5);

                _this3.comp_edgenodes[tri_i].move(p.x, p.y, _this3.g.gamma[tri_i] / 2, _this3.anim_ms);
            };

            for (var si = 0; si < segs.length; si++) {
                var _ret2 = _loop2(si);

                if (_ret2 === 'continue') continue;
            }

            //console.log(pathStr);
            this.comps[ci].stop();
            this.comps[ci].animate({ "d": pathStr + "Z" }, this.anim_ms);
        }
    }]);

    return MeshDraw;
}();

var meshDraw = new MeshDraw("#knot-draw");

var cpWorker = new Worker("js/cp_worker.js");

function drawMapAsync(sigma) {
    var cross_bend = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 8;

    cpWorker.postMessage({
        function: "setLinkDiagram",
        arguments: [sigma, cross_bend]
    });
}

var cpWorkerFunctions = {
    setEmbedding: function setEmbedding(flat_poly, embedding, m4v, conv) {
        meshDraw.clear();
        meshDraw.set_embedding(flat_poly, embedding, m4v, conv);
    },

    updateEmbedding: function updateEmbedding(flat_poly, embedding, m4v, conv) {
        meshDraw.update_embedding(flat_poly, embedding, m4v, conv);
    },

    setLinkDiagram: function setLinkDiagram(link_diagram) {
        meshDraw.set_link_diagram(link_diagram);
    }
};

cpWorker.onmessage = function (ev) {
    //console.log("hello", ev);
    cpWorkerFunctions[ev.data.function].apply(this, ev.data.arguments);
};

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

//let sigma = [[1, 36, 2, 35], [0, 38, 75, 37], [3, 18, 4, 17], [5, 28, 6, 27], [7, 21, 8, 22], [9, 31, 10, 32], [11, 34, 12, 33], [13, 24, 14, 23], [15, 25, 16, 26], [19, 30, 20, 29], [39, 73, 40, 74], [41, 72, 42, 71], [43, 66, 44, 65], [45, 51, 46, 52], [47, 61, 48, 62], [49, 60, 50, 59], [53, 63, 54, 64], [55, 69, 56, 70], [57, 68, 58, 67]];

var sigma = [[1, 3, 2, 4], [0, 60, 59, 119], [5, 15, 6, 16], [7, 84, 8, 83], [9, 101, 10, 102], [11, 92, 12, 91], [13, 89, 14, 90], [17, 62, 18, 61], [19, 79, 20, 80], [21, 70, 22, 69], [23, 42, 24, 41], [25, 28, 26, 27], [29, 32, 30, 31], [33, 52, 34, 51], [35, 37, 36, 38], [39, 49, 40, 50], [43, 54, 44, 53], [45, 48, 46, 47], [55, 67, 56, 68], [57, 82, 58, 81], [63, 77, 64, 78], [65, 72, 66, 71], [73, 75, 74, 76], [85, 88, 86, 87], [93, 99, 94, 100], [95, 98, 96, 97], [103, 118, 104, 117], [105, 112, 106, 111], [107, 109, 108, 110], [113, 115, 114, 116]];

// Complicated 2-link
//let sigma = [[1, 4, 2, 3], [0, 66, 295, 65], [5, 64, 6, 63], [7, 49, 8, 50], [9, 12, 10, 11], [13, 43, 14, 44], [15, 46, 16, 45], [17, 19, 18, 20], [21, 24, 22, 23], [25, 32, 26, 31], [27, 30, 28, 29], [33, 47, 34, 48], [35, 38, 36, 37], [39, 42, 40, 41], [51, 54, 52, 53], [55, 294, 56, 293], [57, 188, 58, 187], [59, 289, 60, 290], [61, 175, 62, 176], [67, 174, 68, 173], [69, 140, 70, 139], [71, 278, 72, 277], [73, 131, 74, 132], [75, 78, 76, 77], [79, 117, 80, 118], [81, 87, 82, 88], [83, 86, 84, 85], [89, 104, 90, 103], [91, 94, 92, 93], [95, 98, 96, 97], [99, 101, 100, 102], [105, 108, 106, 107], [109, 111, 110, 112], [113, 256, 114, 255], [115, 257, 116, 258], [119, 122, 120, 121], [123, 134, 124, 133], [125, 275, 126, 276], [127, 137, 128, 138], [129, 244, 130, 243], [135, 269, 136, 270], [141, 280, 142, 279], [143, 145, 144, 146], [147, 241, 148, 242], [149, 151, 150, 152], [153, 155, 154, 156], [157, 192, 158, 191], [159, 169, 160, 170], [161, 163, 162, 164], [165, 168, 166, 167], [171, 189, 172, 190], [177, 184, 178, 183], [179, 181, 180, 182], [185, 292, 186, 291], [193, 227, 194, 228], [195, 198, 196, 197], [199, 285, 200, 286], [201, 208, 202, 207], [203, 222, 204, 221], [205, 211, 206, 212], [209, 224, 210, 223], [213, 219, 214, 220], [215, 217, 216, 218], [225, 284, 226, 283], [229, 287, 230, 288], [231, 234, 232, 233], [235, 298, 236, 297], [237, 299, 238, 296], [239, 282, 240, 281], [245, 268, 246, 267], [247, 249, 248, 250], [251, 265, 252, 266], [253, 260, 254, 259], [261, 264, 262, 263], [271, 274, 272, 273]];

// More complicated
//let sigma = [[0, 41, 99, 42], [1, 47, 2, 48], [3, 34, 4, 33], [5, 96, 6, 95], [7, 13, 8, 14], [9, 68, 10, 67], [11, 69, 12, 70], [15, 62, 16, 61], [17, 59, 18, 60], [19, 85, 20, 86], [21, 27, 22, 28], [23, 82, 24, 81], [25, 83, 26, 84], [29, 88, 30, 87], [31, 93, 32, 94], [35, 46, 36, 45], [37, 43, 38, 44], [39, 98, 40, 97], [49, 76, 50, 75], [51, 73, 52, 74], [53, 72, 54, 71], [55, 65, 56, 66], [57, 64, 58, 63], [77, 92, 78, 91], [79, 89, 80, 90]];

//let sigma = [[1, 44, 2, 59], [0, 58, 43, 57], [3, 41, 4, 42], [5, 16, 6, 15], [7, 25, 8, 26], [9, 32, 10, 31], [11, 54, 12, 53], [13, 47, 14, 48], [17, 36, 18, 35], [19, 37, 20, 38], [21, 40, 22, 39], [23, 33, 24, 34], [27, 50, 28, 49], [29, 51, 30, 52], [45, 56, 46, 55]];

//let sigma = [[1, 3, 2, 4], [0, 48, 47, 59], [5, 11, 6, 12], [7, 9, 8, 10], [13, 43, 14, 44], [15, 25, 16, 26], [17, 24, 18, 23], [19, 22, 20, 21], [27, 42, 28, 41], [29, 31, 30, 32], [33, 35, 34, 36], [37, 40, 38, 39], [45, 58, 46, 57], [49, 51, 50, 52], [53, 55, 54, 56]];

// Even moreso
//let sigma = [[1, 287, 2, 288], [0, 289, 299, 290], [3, 126, 4, 125], [5, 16, 6, 15], [7, 13, 8, 14], [9, 255, 10, 256], [11, 254, 12, 253], [17, 127, 18, 128], [19, 282, 20, 281], [21, 132, 22, 131], [23, 277, 24, 278], [25, 235, 26, 236], [27, 222, 28, 221], [29, 223, 30, 224], [31, 38, 32, 37], [33, 228, 34, 227], [35, 225, 36, 226], [39, 230, 40, 229], [41, 156, 42, 155], [43, 157, 44, 158], [45, 164, 46, 163], [47, 165, 48, 166], [49, 56, 50, 55], [51, 198, 52, 197], [53, 195, 54, 196], [57, 243, 58, 244], [59, 242, 60, 241], [61, 231, 62, 232], [63, 234, 64, 233], [65, 248, 66, 247], [67, 249, 68, 250], [69, 276, 70, 275], [71, 133, 72, 134], [73, 83, 74, 84], [75, 82, 76, 81], [77, 283, 78, 284], [79, 286, 80, 285], [85, 295, 86, 296], [87, 294, 88, 293], [89, 291, 90, 292], [91, 298, 92, 297], [93, 136, 94, 135], [95, 109, 96, 110], [97, 108, 98, 107], [99, 122, 100, 121], [101, 259, 102, 260], [103, 262, 104, 261], [105, 119, 106, 120], [111, 273, 112, 274], [113, 272, 114, 271], [115, 265, 116, 266], [117, 264, 118, 263], [123, 137, 124, 138], [129, 280, 130, 279], [139, 257, 140, 258], [141, 268, 142, 267], [143, 269, 144, 270], [145, 252, 146, 251], [147, 237, 148, 238], [149, 216, 150, 215], [151, 217, 152, 218], [153, 220, 154, 219], [159, 210, 160, 209], [161, 211, 162, 212], [167, 190, 168, 189], [169, 187, 170, 188], [171, 205, 172, 206], [173, 200, 174, 199], [175, 201, 176, 202], [177, 204, 178, 203], [179, 186, 180, 185], [181, 191, 182, 192], [183, 194, 184, 193], [207, 214, 208, 213], [239, 245, 240, 246]];

drawMapAsync(sigma);

document.getElementById("map_submit").onclick = function (ev) {
    document.querySelector(ev.target.dataset.errbox).textContent = "";
    try {
        meshDraw.clear();
        drawMapAsync(JSON.parse(document.querySelector(ev.target.dataset.target).value));
    } catch (err) {
        document.querySelector(ev.target.dataset.errbox).textContent = "Error";
        console.log("Error:", err);
    }
};
//# sourceMappingURL=orthogonal.js.map