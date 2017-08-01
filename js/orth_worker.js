'use strict';

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _shadow = require('lib/shadow.js');

var _shadow2 = _interopRequireDefault(_shadow);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

importScripts('lalolib/lalolib.js');
importScripts('lodash.min.js');

function leastSquares(X /* : Matrix */, Y /* : Matrix */) /* : leastSquares */{
    //console.log("X", X, inv(X), det(X), qr(X, true), Y);
    var betaHat = solve(mul(X, transpose(X)), mul(X, Y));

    return betaHat;
}

var OrthogonalFace = function () {
    function OrthogonalFace(link, arcs) {
        var _this = this;

        var exterior = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : false;

        _classCallCheck(this, OrthogonalFace);

        this.link = link;
        this.arcs = arcs;
        this.exterior = exterior;
        this.edges = this.arcs.map(function (a) {
            return _this.link.edges[a.edge].map(function (b) {
                return b.index;
            });
        });

        this.turns = this.arcs.map(function (a) {
            return 1;
        });
    }

    _createClass(OrthogonalFace, [{
        key: 'edgeOfIntersection',
        value: function edgeOfIntersection(other) {
            return this.edges.filter(function (e) {
                return other.edges.some(function (oe) {
                    return oe[0] == e[0] && oe[1] == e[1];
                });
            });
        }
    }, {
        key: 'sourceCapacity',
        value: function sourceCapacity() {
            return this.exterior ? 0 : Math.max(4 - this.arcs.length, 0);
        }
    }, {
        key: 'sinkCapacity',
        value: function sinkCapacity() {
            return this.exterior ? this.arcs.length + 4 : Math.max(this.arcs.length - 4, 0);
        }
    }, {
        key: 'bend',
        value: function bend(arc, turns) {
            var i = this.arcs.indexOf(arc);
            turns.reverse();

            var _iteratorNormalCompletion = true;
            var _didIteratorError = false;
            var _iteratorError = undefined;

            try {
                for (var _iterator = turns[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
                    var t = _step.value;

                    var nArc = this.link.verts[arc.vert][(arc.vertpos - 1) % 2];
                    var oArc = this.link.edges[nArc.edge][(nArc.edgepos + 1) % 2];

                    this.arcs.splice(i, 0, oArc);
                    this.turns.splice(i, 0, t);
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
    }]);

    return OrthogonalFace;
}();

var DiGraph = function () {
    /* Lifted from networkx's network_simplex routines */

    function DiGraph() {
        _classCallCheck(this, DiGraph);

        this.nodes = new Map();
        this.edges = new Map();
    }

    _createClass(DiGraph, [{
        key: 'copy',
        value: function copy() {
            var G = new DiGraph();
            G.nodes = new Map(this.nodes);
            G.edges = new Map(this.edges);
            return G;
        }
    }, {
        key: 'addNode',
        value: function addNode(name) {
            var attrs = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : undefined;

            if (attrs === undefined) {
                attrs = {};
            }
            if (this.nodes.has(name)) {
                var _iteratorNormalCompletion2 = true;
                var _didIteratorError2 = false;
                var _iteratorError2 = undefined;

                try {
                    for (var _iterator2 = Object.entries(attrs)[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                        var _ref = _step2.value;

                        var _ref2 = _slicedToArray(_ref, 2);

                        var k = _ref2[0];
                        var _v = _ref2[1];

                        this.nodes.get(name)[k] = _v;
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
            } else {
                this.nodes.set(name, { name: name });
                var _iteratorNormalCompletion3 = true;
                var _didIteratorError3 = false;
                var _iteratorError3 = undefined;

                try {
                    for (var _iterator3 = Object.entries(attrs)[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
                        var _ref3 = _step3.value;

                        var _ref4 = _slicedToArray(_ref3, 2);

                        var _k = _ref4[0];
                        var _v2 = _ref4[1];

                        this.nodes.get(name)[_k] = _v2;
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
            }
            if (!this.edges.has(name)) {
                this.edges.set(name, new Map());
            }
        }
    }, {
        key: 'removeNode',
        value: function removeNode(name) {
            this.nodes.delete(name);
            this.edges.delete(name);
            var _iteratorNormalCompletion4 = true;
            var _didIteratorError4 = false;
            var _iteratorError4 = undefined;

            try {
                for (var _iterator4 = this.edges[Symbol.iterator](), _step4; !(_iteratorNormalCompletion4 = (_step4 = _iterator4.next()).done); _iteratorNormalCompletion4 = true) {
                    var _ref5 = _step4.value;

                    var _ref6 = _slicedToArray(_ref5, 2);

                    var source = _ref6[0];
                    var sinks = _ref6[1];

                    sinks.delete(name);
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
        key: 'addEdge',
        value: function addEdge(source, sink) {
            var attrs = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : undefined;

            if (!this.nodes.has(source)) {
                this.addNode(source);
            }
            if (!this.nodes.has(sink)) {
                this.addNode(sink);
            }
            this.edges.get(source).set(sink, { source: source, sink: sink });

            if (attrs === undefined) {
                return;
            }

            var _iteratorNormalCompletion5 = true;
            var _didIteratorError5 = false;
            var _iteratorError5 = undefined;

            try {
                for (var _iterator5 = Object.entries(attrs)[Symbol.iterator](), _step5; !(_iteratorNormalCompletion5 = (_step5 = _iterator5.next()).done); _iteratorNormalCompletion5 = true) {
                    var _ref7 = _step5.value;

                    var _ref8 = _slicedToArray(_ref7, 2);

                    var k = _ref8[0];
                    var _v3 = _ref8[1];

                    this.edges.get(source).get(sink)[k] = _v3;
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
    }, {
        key: 'removeEdge',
        value: function removeEdge(source, sink) {
            this.edges.get(source).delete(sink);
        }
    }, {
        key: 'getEdge',
        value: function getEdge(source, sink) {
            return this.edges.get(source).get(sink);
        }
    }, {
        key: 'hasEdge',
        value: function hasEdge(source, sink) {
            if (!this.edges.has(source)) {
                return false;
            }
            if (!this.edges.get(source).has(sink)) {
                return false;
            }
            return true;
        }
    }, {
        key: 'getReverseEdges',
        value: function getReverseEdges() {
            var rev_edges = new Map();
            var _iteratorNormalCompletion6 = true;
            var _didIteratorError6 = false;
            var _iteratorError6 = undefined;

            try {
                for (var _iterator6 = this.nodes[Symbol.iterator](), _step6; !(_iteratorNormalCompletion6 = (_step6 = _iterator6.next()).done); _iteratorNormalCompletion6 = true) {
                    var _ref9 = _step6.value;

                    var _ref10 = _slicedToArray(_ref9, 2);

                    var _v4 = _ref10[0];
                    var d = _ref10[1];

                    rev_edges.set(_v4, new Map());
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

            var _iteratorNormalCompletion7 = true;
            var _didIteratorError7 = false;
            var _iteratorError7 = undefined;

            try {
                for (var _iterator7 = this.edges[Symbol.iterator](), _step7; !(_iteratorNormalCompletion7 = (_step7 = _iterator7.next()).done); _iteratorNormalCompletion7 = true) {
                    var _ref11 = _step7.value;

                    var _ref12 = _slicedToArray(_ref11, 2);

                    var _u = _ref12[0];
                    var sinks = _ref12[1];
                    var _iteratorNormalCompletion8 = true;
                    var _didIteratorError8 = false;
                    var _iteratorError8 = undefined;

                    try {
                        for (var _iterator8 = sinks[Symbol.iterator](), _step8; !(_iteratorNormalCompletion8 = (_step8 = _iterator8.next()).done); _iteratorNormalCompletion8 = true) {
                            var _ref13 = _step8.value;

                            var _ref14 = _slicedToArray(_ref13, 2);

                            var _v5 = _ref14[0];
                            var _d2 = _ref14[1];

                            if (!rev_edges.has(_v5)) {
                                rev_edges.set(_v5, new Map());
                            }
                            rev_edges.get(_v5).set(_u, {});
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

            return rev_edges;
        }
    }, {
        key: 'treePath',
        value: function treePath(start, stop) {
            // Requires that this is a tree to avoid real programming...

            // this.edges is a map source -> sink
            var rev_edges = this.getReverseEdges();

            var search_nodes = new Set();
            var _iteratorNormalCompletion9 = true;
            var _didIteratorError9 = false;
            var _iteratorError9 = undefined;

            try {
                for (var _iterator9 = this.nodes.keys()[Symbol.iterator](), _step9; !(_iteratorNormalCompletion9 = (_step9 = _iterator9.next()).done); _iteratorNormalCompletion9 = true) {
                    var _u2 = _step9.value;

                    search_nodes.add(_u2);
                }

                // DFS for now, since it's easier
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

            var path = void 0;
            var search_stack = [[[start], start]];

            do {
                console.log(search_stack, search_nodes);

                var _search_stack$pop = search_stack.pop(),
                    _search_stack$pop2 = _slicedToArray(_search_stack$pop, 2),
                    _path = _search_stack$pop2[0],
                    node = _search_stack$pop2[1];

                search_nodes.delete(node);
                var _iteratorNormalCompletion10 = true;
                var _didIteratorError10 = false;
                var _iteratorError10 = undefined;

                try {
                    for (var _iterator10 = this.edges.get(node)[Symbol.iterator](), _step10; !(_iteratorNormalCompletion10 = (_step10 = _iterator10.next()).done); _iteratorNormalCompletion10 = true) {
                        var _ref15 = _step10.value;

                        var _ref16 = _slicedToArray(_ref15, 2);

                        var _v6 = _ref16[0];
                        var d = _ref16[1];

                        if (_v6 == stop) {
                            return _path.concat([_v6]);
                        }
                        if (search_nodes.has(_v6)) {
                            search_stack.push([_path.concat([_v6]), _v6]);
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

                console.log(node, rev_edges, this.edges);
                var _iteratorNormalCompletion11 = true;
                var _didIteratorError11 = false;
                var _iteratorError11 = undefined;

                try {
                    for (var _iterator11 = rev_edges.get(node)[Symbol.iterator](), _step11; !(_iteratorNormalCompletion11 = (_step11 = _iterator11.next()).done); _iteratorNormalCompletion11 = true) {
                        var _ref17 = _step11.value;

                        var _ref18 = _slicedToArray(_ref17, 2);

                        var _v7 = _ref18[0];
                        var _d3 = _ref18[1];

                        if (_v7 == stop) {
                            return _path.concat([_v7]);
                        }
                        if (search_nodes.has(_v7)) {
                            search_stack.push([_path.concat([_v7]), _v7]);
                        }
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
            } while (search_nodes.size > 0);

            return path;
        }
    }, {
        key: 'connectedComponentSubgraphs',
        value: function connectedComponentSubgraphs() {
            // Requires that this is a tree to avoid real programming...
            // CAVEAT: We only care about nodes so "meh" to edges

            // this.edges is a map source -> sink
            var rev_edges = this.getReverseEdges();

            var search_nodes = new Set();
            var _iteratorNormalCompletion12 = true;
            var _didIteratorError12 = false;
            var _iteratorError12 = undefined;

            try {
                for (var _iterator12 = this.nodes.keys()[Symbol.iterator](), _step12; !(_iteratorNormalCompletion12 = (_step12 = _iterator12.next()).done); _iteratorNormalCompletion12 = true) {
                    var _u3 = _step12.value;

                    search_nodes.add(_u3);
                }

                // DFS for now, since it's easier
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

            var search_stack = [];
            var components = [];

            while (search_nodes) {
                search_stack = [search_nodes.values().next()];

                var G = new DiGraph();
                while (search_stack) {
                    var node = search_stack.pop();
                    search_nodes.delete(node);
                    G.addNode(node);
                    var _iteratorNormalCompletion13 = true;
                    var _didIteratorError13 = false;
                    var _iteratorError13 = undefined;

                    try {
                        for (var _iterator13 = this.edges.get(node)[Symbol.iterator](), _step13; !(_iteratorNormalCompletion13 = (_step13 = _iterator13.next()).done); _iteratorNormalCompletion13 = true) {
                            var _ref19 = _step13.value;

                            var _ref20 = _slicedToArray(_ref19, 2);

                            var _v8 = _ref20[0];
                            var d = _ref20[1];

                            if (search_nodes.has(_v8)) {
                                search_stack.push(_v8);
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

                    var _iteratorNormalCompletion14 = true;
                    var _didIteratorError14 = false;
                    var _iteratorError14 = undefined;

                    try {
                        for (var _iterator14 = rev_edges.get(node)[Symbol.iterator](), _step14; !(_iteratorNormalCompletion14 = (_step14 = _iterator14.next()).done); _iteratorNormalCompletion14 = true) {
                            var _ref21 = _step14.value;

                            var _ref22 = _slicedToArray(_ref21, 2);

                            var _v9 = _ref22[0];
                            var _d4 = _ref22[1];

                            if (search_nodes.has(_v9)) {
                                search_stack.push(_v9);
                            }
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
                components.push(G);
            }

            return components;
        }
    }, {
        key: 'minCostFlow',
        value: function minCostFlow() {
            var _initialTreeSolution2 = this._initialTreeSolution(),
                _initialTreeSolution3 = _slicedToArray(_initialTreeSolution2, 6),
                H = _initialTreeSolution3[0],
                T = _initialTreeSolution3[1],
                y = _initialTreeSolution3[2],
                artificialEdges = _initialTreeSolution3[3],
                flowCost = _initialTreeSolution3[4],
                r = _initialTreeSolution3[5];

            console.log(H, T);

            var c = [];
            var _iteratorNormalCompletion15 = true;
            var _didIteratorError15 = false;
            var _iteratorError15 = undefined;

            try {
                for (var _iterator15 = H.edges[Symbol.iterator](), _step15; !(_iteratorNormalCompletion15 = (_step15 = _iterator15.next()).done); _iteratorNormalCompletion15 = true) {
                    var _ref23 = _step15.value;

                    var _ref24 = _slicedToArray(_ref23, 2);

                    var _u4 = _ref24[0];
                    var sinks = _ref24[1];
                    var _iteratorNormalCompletion22 = true;
                    var _didIteratorError22 = false;
                    var _iteratorError22 = undefined;

                    try {
                        for (var _iterator22 = sinks[Symbol.iterator](), _step22; !(_iteratorNormalCompletion22 = (_step22 = _iterator22.next()).done); _iteratorNormalCompletion22 = true) {
                            var _ref31 = _step22.value;

                            var _ref32 = _slicedToArray(_ref31, 2);

                            var _v15 = _ref32[0];
                            var _d5 = _ref32[1];

                            c[[_u4, _v15]] = _.get(_d5, 'weight', 0);
                        }
                    } catch (err) {
                        _didIteratorError22 = true;
                        _iteratorError22 = err;
                    } finally {
                        try {
                            if (!_iteratorNormalCompletion22 && _iterator22.return) {
                                _iterator22.return();
                            }
                        } finally {
                            if (_didIteratorError22) {
                                throw _iteratorError22;
                            }
                        }
                    }
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

            var reverse = void 0;
            while (true) {
                var newEdge = H._findEnteringEdge(c);
                if (newEdge.length != 2) {
                    break;
                }

                var cycleCost = Math.abs(c[newEdge]);

                console.log("!!!", r, newEdge);
                var path1 = T.treePath(r, newEdge[0]);
                var path2 = T.treePath(r, newEdge[1]);

                var join = r;

                for (var _i = 1; _i < path1.length; _i++) {
                    var node = path1[_i];
                    if (_i + 1 < path2.length && node == path2[_i + 1]) {
                        join = node;
                    } else {
                        break;
                    }
                }

                path1 = path1.slice(path1.indexOf(join));
                path2 = path2.slice(path2.indexOf(join));
                var cycle = [];

                if (_.get(H.getEdge.apply(H, _toConsumableArray(newEdge)), "flow", 0) == 0) {
                    reverse = false;
                    path2.reverse();
                    cycle = path1.concat(path2);
                } else {
                    reverse = true;
                    path1.reverse();
                    cycle = path2.concat(path1);
                }

                var _H$_findLeavingEdge = H._findLeavingEdge(T, cycle, newEdge, reverse),
                    _H$_findLeavingEdge2 = _slicedToArray(_H$_findLeavingEdge, 2),
                    leavingEdge = _H$_findLeavingEdge2[0],
                    eps = _H$_findLeavingEdge2[1];

                if (eps != 0) {
                    flowCost -= cycleCost * eps;
                    if (cycle.length == 3) {
                        if (reverse) {
                            eps = -eps;
                        }

                        var _newEdge = _slicedToArray(newEdge, 2),
                            _u5 = _newEdge[0],
                            _v10 = _newEdge[1];

                        _.set(H.getEdge(_u5, _v10), _.get(H.getEdge(_u5, _v10), "flow", 0) + eps);
                        _.set(H.getEdge(_v10, _u5), _.get(H.getEdge(_v10, _u5), "flow", 0) + eps);
                    } else {
                        for (var j = 0; j < cycle.length - 1; j++) {
                            v = cycle[i + 1];
                            if ([u, v] == newEdge || T.hasEdge(u, v)) {
                                _.set(H.getEdge(u, v), _.get(H.getEdge(u, v), "flow", 0) + eps);
                            } else {
                                _.set(H.getEdge(v, u), _.get(H.getEdge(v, u), "flow", 0) - eps);
                            }
                        }
                    }
                }

                T.addEdge.apply(T, _toConsumableArray(newEdge));
                T.removeEdge.apply(T, _toConsumableArray(leavingEdge));

                if (newEdge != leavingEdge) {
                    var forest = T.copy();
                    forest.removeEdge.apply(forest, _toConsumableArray(newEdge));

                    var _forest$connectedComp = forest.connectedComponentSubgraphs(),
                        _forest$connectedComp2 = _slicedToArray(_forest$connectedComp, 2),
                        R = _forest$connectedComp2[0],
                        notR = _forest$connectedComp2[1];

                    if (notR.nodes.has(r)) {
                        var tmp = R;
                        R = notR;
                        notR = tmp;
                    }

                    if (R.nodes.has(newEdge[0])) {
                        var _iteratorNormalCompletion16 = true;
                        var _didIteratorError16 = false;
                        var _iteratorError16 = undefined;

                        try {
                            for (var _iterator16 = notR.nodes.keys()[Symbol.iterator](), _step16; !(_iteratorNormalCompletion16 = (_step16 = _iterator16.next()).done); _iteratorNormalCompletion16 = true) {
                                var _v11 = _step16.value;

                                y[_v11] += c[newEdge];
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
                    } else {
                        var _iteratorNormalCompletion17 = true;
                        var _didIteratorError17 = false;
                        var _iteratorError17 = undefined;

                        try {
                            for (var _iterator17 = notR.nodes.keys()[Symbol.iterator](), _step17; !(_iteratorNormalCompletion17 = (_step17 = _iterator17.next()).done); _iteratorNormalCompletion17 = true) {
                                var _v12 = _step17.value;

                                y[_v12] -= c[newEdge];
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

                    var _iteratorNormalCompletion18 = true;
                    var _didIteratorError18 = false;
                    var _iteratorError18 = undefined;

                    try {
                        for (var _iterator18 = H.edges[Symbol.iterator](), _step18; !(_iteratorNormalCompletion18 = (_step18 = _iterator18.next()).done); _iteratorNormalCompletion18 = true) {
                            var _ref25 = _step18.value;

                            var _ref26 = _slicedToArray(_ref25, 2);

                            var _u6 = _ref26[0];
                            var _sinks = _ref26[1];
                            var _iteratorNormalCompletion19 = true;
                            var _didIteratorError19 = false;
                            var _iteratorError19 = undefined;

                            try {
                                for (var _iterator19 = _sinks[Symbol.iterator](), _step19; !(_iteratorNormalCompletion19 = (_step19 = _iterator19.next()).done); _iteratorNormalCompletion19 = true) {
                                    var _ref27 = _step19.value;

                                    var _ref28 = _slicedToArray(_ref27, 2);

                                    var _v13 = _ref28[0];
                                    var d = _ref28[1];

                                    if (notR.nodes.has(_u6) || notR.nodes.has(_v13)) {
                                        c[[_u6, _v13]] = _.get(H.getEdge(_u6, _v13), "weight", 0) + y[_u6] - y[_v13];
                                    }
                                }
                            } catch (err) {
                                _didIteratorError19 = true;
                                _iteratorError19 = err;
                            } finally {
                                try {
                                    if (!_iteratorNormalCompletion19 && _iterator19.return) {
                                        _iterator19.return();
                                    }
                                } finally {
                                    if (_didIteratorError19) {
                                        throw _iteratorError19;
                                    }
                                }
                            }
                        }
                    } catch (err) {
                        _didIteratorError18 = true;
                        _iteratorError18 = err;
                    } finally {
                        try {
                            if (!_iteratorNormalCompletion18 && _iterator18.return) {
                                _iterator18.return();
                            }
                        } finally {
                            if (_didIteratorError18) {
                                throw _iteratorError18;
                            }
                        }
                    }
                }
            }

            var _iteratorNormalCompletion20 = true;
            var _didIteratorError20 = false;
            var _iteratorError20 = undefined;

            try {
                for (var _iterator20 = artificialEdges[Symbol.iterator](), _step20; !(_iteratorNormalCompletion20 = (_step20 = _iterator20.next()).done); _iteratorNormalCompletion20 = true) {
                    var _ref29 = _step20.value;

                    var _ref30 = _slicedToArray(_ref29, 2);

                    var _u7 = _ref30[0];
                    var _v14 = _ref30[1];

                    H.removeEdge(_u7, _v14);
                }
            } catch (err) {
                _didIteratorError20 = true;
                _iteratorError20 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion20 && _iterator20.return) {
                        _iterator20.return();
                    }
                } finally {
                    if (_didIteratorError20) {
                        throw _iteratorError20;
                    }
                }
            }

            var _iteratorNormalCompletion21 = true;
            var _didIteratorError21 = false;
            var _iteratorError21 = undefined;

            try {
                for (var _iterator21 = H.nodes.keys()[Symbol.iterator](), _step21; !(_iteratorNormalCompletion21 = (_step21 = _iterator21.next()).done); _iteratorNormalCompletion21 = true) {
                    var _u8 = _step21.value;

                    if (!this.nodes.has(_u8)) {
                        H.removeNode(_u8);
                    }
                }
            } catch (err) {
                _didIteratorError21 = true;
                _iteratorError21 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion21 && _iterator21.return) {
                        _iterator21.return();
                    }
                } finally {
                    if (_didIteratorError21) {
                        throw _iteratorError21;
                    }
                }
            }

            var flowDict = this._createFlowDict(H);

            return [flowCost, flowDict];
        }
    }, {
        key: '_createFlowDict',
        value: function _createFlowDict(H) {
            var flowDict = new Map();

            var _iteratorNormalCompletion23 = true;
            var _didIteratorError23 = false;
            var _iteratorError23 = undefined;

            try {
                for (var _iterator23 = this.edges[Symbol.iterator](), _step23; !(_iteratorNormalCompletion23 = (_step23 = _iterator23.next()).done); _iteratorNormalCompletion23 = true) {
                    var _ref33 = _step23.value;

                    var _ref34 = _slicedToArray(_ref33, 2);

                    var _u9 = _ref34[0];
                    var sinks = _ref34[1];

                    flowDict.set(_u9, new Map());
                    var _iteratorNormalCompletion24 = true;
                    var _didIteratorError24 = false;
                    var _iteratorError24 = undefined;

                    try {
                        for (var _iterator24 = sinks[Symbol.iterator](), _step24; !(_iteratorNormalCompletion24 = (_step24 = _iterator24.next()).done); _iteratorNormalCompletion24 = true) {
                            var _ref35 = _step24.value;

                            var _ref36 = _slicedToArray(_ref35, 2);

                            var _v16 = _ref36[0];
                            var d = _ref36[1];

                            if (H.hasEdge(_u9, _v16)) {
                                flowDict.get(_u9).set(_v16, _.get(d, "flow", 0));
                            } else {
                                flowDict.get(_u9).set(_v16, 0);
                            }
                        }
                    } catch (err) {
                        _didIteratorError24 = true;
                        _iteratorError24 = err;
                    } finally {
                        try {
                            if (!_iteratorNormalCompletion24 && _iterator24.return) {
                                _iterator24.return();
                            }
                        } finally {
                            if (_didIteratorError24) {
                                throw _iteratorError24;
                            }
                        }
                    }
                }
            } catch (err) {
                _didIteratorError23 = true;
                _iteratorError23 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion23 && _iterator23.return) {
                        _iterator23.return();
                    }
                } finally {
                    if (_didIteratorError23) {
                        throw _iteratorError23;
                    }
                }
            }

            return flowDict;
        }
    }, {
        key: '_initialTreeSolution',
        value: function _initialTreeSolution() {
            var H = new DiGraph();

            var maxWeight = 0;
            var _iteratorNormalCompletion25 = true;
            var _didIteratorError25 = false;
            var _iteratorError25 = undefined;

            try {
                for (var _iterator25 = this.edges[Symbol.iterator](), _step25; !(_iteratorNormalCompletion25 = (_step25 = _iterator25.next()).done); _iteratorNormalCompletion25 = true) {
                    var _ref37 = _step25.value;

                    var _ref38 = _slicedToArray(_ref37, 2);

                    var _u10 = _ref38[0];
                    var sinks = _ref38[1];
                    var _iteratorNormalCompletion28 = true;
                    var _didIteratorError28 = false;
                    var _iteratorError28 = undefined;

                    try {
                        for (var _iterator28 = sinks[Symbol.iterator](), _step28; !(_iteratorNormalCompletion28 = (_step28 = _iterator28.next()).done); _iteratorNormalCompletion28 = true) {
                            var _ref43 = _step28.value;

                            var _ref44 = _slicedToArray(_ref43, 2);

                            var _v18 = _ref44[0];
                            var _d7 = _ref44[1];

                            if (_.get(_d7, 'capacity', 1) > 0) {
                                H.addEdge(_u10, _v18, _d7);
                                maxWeight = Math.max(maxWeight, _.get(_d7, 'weight', 0));
                            }
                        }
                    } catch (err) {
                        _didIteratorError28 = true;
                        _iteratorError28 = err;
                    } finally {
                        try {
                            if (!_iteratorNormalCompletion28 && _iterator28.return) {
                                _iterator28.return();
                            }
                        } finally {
                            if (_didIteratorError28) {
                                throw _iteratorError28;
                            }
                        }
                    }
                }
            } catch (err) {
                _didIteratorError25 = true;
                _iteratorError25 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion25 && _iterator25.return) {
                        _iterator25.return();
                    }
                } finally {
                    if (_didIteratorError25) {
                        throw _iteratorError25;
                    }
                }
            }

            var _iteratorNormalCompletion26 = true;
            var _didIteratorError26 = false;
            var _iteratorError26 = undefined;

            try {
                for (var _iterator26 = this.nodes[Symbol.iterator](), _step26; !(_iteratorNormalCompletion26 = (_step26 = _iterator26.next()).done); _iteratorNormalCompletion26 = true) {
                    var _ref39 = _step26.value;

                    var _ref40 = _slicedToArray(_ref39, 2);

                    var _n = _ref40[0];
                    var d = _ref40[1];

                    if (_.get(d, 'demand', 0) != 0) {
                        H.addNode(_n, d);
                    }
                }
            } catch (err) {
                _didIteratorError26 = true;
                _iteratorError26 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion26 && _iterator26.return) {
                        _iterator26.return();
                    }
                } finally {
                    if (_didIteratorError26) {
                        throw _iteratorError26;
                    }
                }
            }

            var nodeIter = H.nodes.entries();

            var _nodeIter$next$value = _slicedToArray(nodeIter.next().value, 2),
                r = _nodeIter$next$value[0],
                _d = _nodeIter$next$value[1];

            var T = new DiGraph();
            var y = { r: 0 };
            var artificialEdges = [];
            var flowCost = 0;

            var n = H.nodes.size;
            var hugeWeight = 1 + n * maxWeight;

            var _iteratorNormalCompletion27 = true;
            var _didIteratorError27 = false;
            var _iteratorError27 = undefined;

            try {
                for (var _iterator27 = nodeIter[Symbol.iterator](), _step27; !(_iteratorNormalCompletion27 = (_step27 = _iterator27.next()).done); _iteratorNormalCompletion27 = true) {
                    var _ref41 = _step27.value;

                    var _ref42 = _slicedToArray(_ref41, 2);

                    var _v17 = _ref42[0];
                    var _d6 = _ref42[1];

                    var vDemand = _.get(_d6, 'demand', 0);
                    if (vDemand >= 0) {
                        if (!H.hasEdge(r, _v17)) {
                            H.addEdge(r, _v17, { weight: hugeWeight, flow: vDemand });
                            artificialEdges.push([r, _v17]);
                            y[_v17] = hugeWeight;
                            T.addEdge(r, _v17);
                            flowCost += vDemand * hugeWeight;
                        } else {
                            if (!"capacity" in H.getEdge(r, _v17) || vDemand <= H.getEdge(r, _v17).capacity) {
                                H.getEdge(r, _v17).flow = vDemand;
                                y[_v17] = _.get(H.getEdge(r, _v17), "weight", 0);
                                T.addEdge(r, _v17);
                                flowCost += vDemand * y[_v17];
                            } else {
                                var newLabel = -H.nodes.size;
                                console.log(newLabel);
                                H.addEdge(r, newLabel, { weight: hugeWeight, flow: vDemand });
                                H.addEdge(newLabel, _v17, { weight: hugeWeight, flow: vDemand });
                                artificialEdges.push([r, newLabel]);
                                artificialEdges.push([newLabel, _v17]);
                                y[_v17] = 2 * hugeWeight;
                                y[newLabel] = hugeWeight;
                                T.addEdge(r, newLabel);
                                T.addEdge(newLabel, _v17);
                                flowCost += vDemand * y[_v17];
                            }
                        }
                    } else {
                        // vDemand < 0
                        if (!H.hasEdge(_v17, r)) {
                            H.addEdge(_v17, r, { weight: hugeWeight, flow: -vDemand });
                            artificialEdges.push([_v17, r]);
                            y[_v17] = -hugeWeight;
                            T.addEdge(_v17, r);
                            flowCost += vDemand * hugeWeight;
                        } else {
                            if (!"capacity" in H.getEdge(_v17, r) || -vDemand <= H.getEdge(_v17, r).capacity) {
                                H.getEdge(_v17, r).flow = vDemand;
                                y[_v17] = -_.get(H.getEdge(_v17, r), "weight", 0);
                                T.add_edge(_v17, r);
                                flowCost += vDemand * y[_v17];
                            } else {
                                var _newLabel = -H.nodes.size;
                                H.addEdge(_v17, _newLabel, { weight: hugeWeight, flow: -vDemand });
                                H.addEdge(_newLabel, r, { weight: hugeWeight, flow: -vDemand });
                                artificialEdges.push([_v17, _newLabel]);
                                artificialEdges.push([_newLabel, r]);
                                y[_v17] = -2 * hugeWeight;
                                y[_newLabel] = -hugeWeight;
                                T.addEdge(_v17, _newLabel);
                                T.addEdge(_newLabel, r);
                                flowCost += vDemand * y[_v17];
                            }
                        }
                    }

                    return [H, T, y, artificialEdges, flowCost, r];
                }
            } catch (err) {
                _didIteratorError27 = true;
                _iteratorError27 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion27 && _iterator27.return) {
                        _iterator27.return();
                    }
                } finally {
                    if (_didIteratorError27) {
                        throw _iteratorError27;
                    }
                }
            }
        }
    }, {
        key: '_findEnteringEdge',
        value: function _findEnteringEdge(c) {
            var newEdge = [];

            console.log("-->", this.edges);
            var _iteratorNormalCompletion29 = true;
            var _didIteratorError29 = false;
            var _iteratorError29 = undefined;

            try {
                for (var _iterator29 = this.edges[Symbol.iterator](), _step29; !(_iteratorNormalCompletion29 = (_step29 = _iterator29.next()).done); _iteratorNormalCompletion29 = true) {
                    var _ref45 = _step29.value;

                    var _ref46 = _slicedToArray(_ref45, 2);

                    var _u11 = _ref46[0];
                    var sinks = _ref46[1];
                    var _iteratorNormalCompletion30 = true;
                    var _didIteratorError30 = false;
                    var _iteratorError30 = undefined;

                    try {
                        for (var _iterator30 = sinks[Symbol.iterator](), _step30; !(_iteratorNormalCompletion30 = (_step30 = _iterator30.next()).done); _iteratorNormalCompletion30 = true) {
                            var _ref47 = _step30.value;

                            var _ref48 = _slicedToArray(_ref47, 2);

                            var _v19 = _ref48[0];
                            var d = _ref48[1];

                            if (_.get(d, 'flow', 0) == 0) {
                                if (c[[_u11, _v19]] < 0) {
                                    newEdge = [_u11, _v19];
                                    break;
                                }
                            } else {
                                if ("capacity" in d && _.get(d, 'flow', 0) == d.capacity && c[[_u11, _v19]] > 0) {
                                    newEdge = [_u11, _v19];
                                    break;
                                }
                            }
                        }
                    } catch (err) {
                        _didIteratorError30 = true;
                        _iteratorError30 = err;
                    } finally {
                        try {
                            if (!_iteratorNormalCompletion30 && _iterator30.return) {
                                _iterator30.return();
                            }
                        } finally {
                            if (_didIteratorError30) {
                                throw _iteratorError30;
                            }
                        }
                    }
                }
            } catch (err) {
                _didIteratorError29 = true;
                _iteratorError29 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion29 && _iterator29.return) {
                        _iterator29.return();
                    }
                } finally {
                    if (_didIteratorError29) {
                        throw _iteratorError29;
                    }
                }
            }

            return newEdge;
        }
    }, {
        key: '_findLeavingEdge',
        value: function _findLeavingEdge(T, cycle, newEdge, reverse) {
            var eps = false;
            var leavingEdge = [];

            if (cycle.length == 3) {
                var _newEdge2 = _slicedToArray(newEdge, 2),
                    _u12 = _newEdge2[0],
                    _v20 = _newEdge2[1];

                if (reverse) {
                    if (_.get(this.getEdge(_u12, _v20), "flow", 0) > _.get(this.getEdge(_v20, _u12), "flow", 0)) {
                        return [[_v20, _u12], _.get(this.getEdge(_v20, _u12), "flow", 0)];
                    } else {
                        return [[_u12, _v20], _.get(this.getEdge(_u12, _v20), "flow", 0)];
                    }
                } else {
                    var uv_res = _.get(this.getEdge(_u12, _v20), "capacity", 0) - _.get(this.getEdge(_u12, _v20), "flow", 0);
                    var vu_res = _.get(this.getEdge(_v20, _u12), "capacity", 0) - _.get(this.getEdge(_v20, _u12), "flow", 0);

                    if (uv_res > vu_res) {
                        return [[_v20, _u12], vu_res];
                    } else {
                        return [[_u12, _v20], uv_res];
                    }
                }
            }

            for (var _i2 = 0; _i2 < cycle.length - 1; _i2++) {
                var _u13 = cycle[_i2];

                var edgeCapacity = false;
                var edge = [];
                var _v21 = cycle[_i2 + 1];

                if ([_u13, _v21] == newEdge || T.hasEdge(_u13, _v21)) {
                    if ("capacity" in this.getEdge(_u13, _v21)) {
                        edgeCapacity = this.getEdge(_u13, _v21).capacity - _.get(this.getEdge(_u13, _v21), "flow", 0);
                        edge = [_u13, _v21];
                    }
                } else {
                    edgeCapacity = _.get(this.getEdge(_v21, _u13), "flow", 0);
                    edge = [_v21, _u13];
                }

                if (edge) {
                    if (leavingEdge) {
                        if (edgeCapacity < eps) {
                            eps = edgeCapacity;
                            leavingEdge = edge;
                        }
                    } else {
                        eps = edgeCapacity;
                        leavingEdge = edge;
                    }
                }
            }

            return [leavingEdge, eps];
        }
    }]);

    return DiGraph;
}();

var OrthogonalDiagramEmbedding = function () {
    function OrthogonalDiagramEmbedding(shadow) {
        _classCallCheck(this, OrthogonalDiagramEmbedding);

        this.shadow = shadow;

        this.faces = shadow.faces.map(function (f) {
            return new OrthogonalFace(shadow, f);
        });
        var F = this.faces.reduce(function (bigF, f) {
            return bigF.length >= f.length ? bigF : f;
        });
        F.exterior = true;

        this.faceNetwork = this.flowNetwork();
        this.bend();
        this.orientEdges();

        //this.edges =
        this.repairComponents();
    }

    _createClass(OrthogonalDiagramEmbedding, [{
        key: 'flowNetwork',
        value: function flowNetwork() {
            var G = new DiGraph();

            var sourceDemand = this.faces.reduce(function (tot, f) {
                return tot + f.sourceCapacity();
            }, 0);
            G.addNode('s', { demand: sourceDemand });
            for (var fi = 0; fi < this.faces.length; fi++) {
                var sourceCapacity = this.faces[fi].sourceCapacity();
                if (sourceCapacity > 0) {
                    G.addEdge('s', fi, { weight: 0, capacity: sourceCapacity });
                }
            }

            var sinkDemand = this.faces.reduce(function (tot, f) {
                return tot + f.sinkCapacity();
            }, 0);
            G.addNode('t', { demand: sinkDemand });
            for (var _fi = 0; _fi < this.faces.length; _fi++) {
                var sinkCapacity = this.faces[_fi].sinkCapacity();
                if (sinkCapacity > 0) {
                    G.addEdge(_fi, 't', { weight: 0, capacity: sinkCapacity });
                }
            }

            for (var ai = 0; ai < this.faces.length; ai++) {
                for (var bi = 0; bi < this.faces.length; bi++) {
                    if (ai != bi && this.faces[ai].edgeOfIntersection(this.faces[bi])) {
                        G.addEdge(ai, bi, { weight: 1 });
                    }
                }
            }

            return G;
        }
    }, {
        key: 'bend',
        value: function bend() {
            var flow = this.faceNetwork.minCostFlow()[1];
            console.log(flow);

            var _iteratorNormalCompletion31 = true;
            var _didIteratorError31 = false;
            var _iteratorError31 = undefined;

            try {
                for (var _iterator31 = flow[Symbol.iterator](), _step31; !(_iteratorNormalCompletion31 = (_step31 = _iterator31.next()).done); _iteratorNormalCompletion31 = true) {
                    var _ref49 = _step31.value;

                    var _ref50 = _slicedToArray(_ref49, 2);

                    var a = _ref50[0];
                    var flows = _ref50[1];
                    var _iteratorNormalCompletion32 = true;
                    var _didIteratorError32 = false;
                    var _iteratorError32 = undefined;

                    try {
                        for (var _iterator32 = flows[Symbol.iterator](), _step32; !(_iteratorNormalCompletion32 = (_step32 = _iterator32.next()).done); _iteratorNormalCompletion32 = true) {
                            var _ref51 = _step32.value;

                            var _ref52 = _slicedToArray(_ref51, 2);

                            var b = _ref52[0];
                            var w_a = _ref52[1];

                            if (!w_a || a == 's' || b == 's' || a == 't' || b == 't') {
                                continue;
                            }

                            var w_b = flow.get(b).get(a);

                            var _ref53 = [this.faces[a], this.faces[b]],
                                A = _ref53[0],
                                B = _ref53[1];

                            console.log(a, b, A, B);

                            var _A$edgeOfIntersection = A.edgeOfIntersection(B),
                                _A$edgeOfIntersection2 = _slicedToArray(_A$edgeOfIntersection, 2),
                                e_a = _A$edgeOfIntersection2[0],
                                e_b = _A$edgeOfIntersection2[1];

                            var turnsA = new Array(w_a).map(function (x) {
                                return 1;
                            }).concat(new Array(w_b).map(function (x) {
                                return -1;
                            }));
                            var turnsB = new Array(w_b).map(function (x) {
                                return 1;
                            }).concat(new Array(w_a).map(function (x) {
                                return -1;
                            }));

                            A.bend(e_a, turnsA);
                            B.bend(e_b, turnsB);
                        }
                    } catch (err) {
                        _didIteratorError32 = true;
                        _iteratorError32 = err;
                    } finally {
                        try {
                            if (!_iteratorNormalCompletion32 && _iterator32.return) {
                                _iterator32.return();
                            }
                        } finally {
                            if (_didIteratorError32) {
                                throw _iteratorError32;
                            }
                        }
                    }
                }
            } catch (err) {
                _didIteratorError31 = true;
                _iteratorError31 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion31 && _iterator31.return) {
                        _iterator31.return();
                    }
                } finally {
                    if (_didIteratorError31) {
                        throw _iteratorError31;
                    }
                }
            }
        }
    }, {
        key: 'repairComponents',
        value: function repairComponents() {
            //
        }
    }, {
        key: 'orientEdges',
        value: function orientEdges() {
            //
        }
    }]);

    return OrthogonalDiagramEmbedding;
}();

var ForceLinkDiagram = function () {
    /* Link diagram embedding improved by ImPrEd */
    function ForceLinkDiagram(verts, edges, faces) {
        _classCallCheck(this, ForceLinkDiagram);

        this.verts = verts;
        this.edges = edges;
        this.faces = faces;

        //console.log("+++++++");
        //console.log(verts);
        //console.log(edges);
        //console.log(faces);

        this.adjMap = {};
        var _iteratorNormalCompletion33 = true;
        var _didIteratorError33 = false;
        var _iteratorError33 = undefined;

        try {
            for (var _iterator33 = edges[Symbol.iterator](), _step33; !(_iteratorNormalCompletion33 = (_step33 = _iterator33.next()).done); _iteratorNormalCompletion33 = true) {
                var edge = _step33.value;

                var _edge = _slicedToArray(edge, 2),
                    a = _edge[0],
                    b = _edge[1];

                if (a in this.adjMap) {
                    this.adjMap[a].push(b);
                } else {
                    this.adjMap[a] = [b];
                }

                if (b in this.adjMap) {
                    this.adjMap[b].push(a);
                } else {
                    this.adjMap[b] = [a];
                }
            }
        } catch (err) {
            _didIteratorError33 = true;
            _iteratorError33 = err;
        } finally {
            try {
                if (!_iteratorNormalCompletion33 && _iterator33.return) {
                    _iterator33.return();
                }
            } finally {
                if (_didIteratorError33) {
                    throw _iteratorError33;
                }
            }
        }

        this.delta = 2;
        this.gamma = 5;

        this.dbar = 3 * this.delta;

        this.aExp = 1;
        this.erExp = 2;

        this.calculateSurroundingEdges();
    }

    _createClass(ForceLinkDiagram, [{
        key: 'distance',
        value: function distance(u, v) {
            return norm(sub(u, v));
        }
    }, {
        key: 'forceAvert',
        value: function forceAvert(u, v) {
            //console.log(this.distance(u,v));
            return mul(Math.pow(this.distance(u, v) / this.delta, this.aExp), sub(v, u));
        }
    }, {
        key: 'forceRvert',
        value: function forceRvert(u, v) {
            var d = this.distance(u, v);
            return mul(Math.pow(this.delta / d, this.erExp), sub(u, v));
        }
    }, {
        key: 'computeVe',
        value: function computeVe(v, a, b) {
            var m = (a[1] - b[1]) / (a[0] - b[0]);
            var n = -1 / m;
            var c = a[1] - m * a[0];
            var d = v[1] - n * v[0];

            var x = (d - c) / (m - n);
            return [x, m * x + c];
        }
    }, {
        key: 'veOnEdge',
        value: function veOnEdge(ve, a, b) {
            return (ve[0] <= a[0] && ve[0] >= b[0] || ve[0] <= b[0] && ve[0] >= a[0]) && (ve[1] <= a[1] && ve[1] >= b[1] || ve[1] <= b[1] && ve[1] >= a[1]);
        }
    }, {
        key: 'forceRedge',
        value: function forceRedge(u, a, b, ve) {
            var d = this.distance(u, ve);
            if (d >= this.gamma) {
                // node and "virtual edge" too far
                return [0, 0];
            }

            return mul(-Math.pow(this.gamma - d, this.erExp) / d, sub(ve, u));
        }
    }, {
        key: 'surroundingEdges',
        value: function surroundingEdges(ui) {
            var _this2 = this;

            // calculate the surrounding edges SUi
            var edges = [];

            var _loop = function _loop(face) {
                if (face.includes(ui)) {
                    if (ui == 63) {
                        console.log(face);
                    }

                    var _loop2 = function _loop2(_i3) {
                        console.assert(_this2.edges.filter(function (e) {
                            return e[0] == face[_i3] && e[1] == face[_i3 + 1] || e[1] == face[_i3] && e[0] == face[_i3 + 1];
                        }).length > 0);
                        edges.push([face[_i3], face[_i3 + 1]]);
                    };

                    for (var _i3 = 0; _i3 < face.length - 1; _i3++) {
                        _loop2(_i3);
                    }
                    edges.push([face[face.length - 1], face[0]]);
                }
            };

            var _iteratorNormalCompletion34 = true;
            var _didIteratorError34 = false;
            var _iteratorError34 = undefined;

            try {
                for (var _iterator34 = this.faces[Symbol.iterator](), _step34; !(_iteratorNormalCompletion34 = (_step34 = _iterator34.next()).done); _iteratorNormalCompletion34 = true) {
                    var face = _step34.value;

                    _loop(face);
                }
            } catch (err) {
                _didIteratorError34 = true;
                _iteratorError34 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion34 && _iterator34.return) {
                        _iterator34.return();
                    }
                } finally {
                    if (_didIteratorError34) {
                        throw _iteratorError34;
                    }
                }
            }

            return edges;
        }
    }, {
        key: 'calculateSurroundingEdges',
        value: function calculateSurroundingEdges() {
            this.surrEdges = [];
            for (var _i4 = 0; _i4 < this.verts.length; _i4++) {
                this.surrEdges[_i4] = this.surroundingEdges(_i4);
            }
        }
    }, {
        key: 'move',
        value: function move(ui, FUx, FUy, MU) {
            var i = void 0;
            if (FUx >= 0) {
                if (FUy >= 0) {
                    if (FUx >= FUy) {
                        i = 0;
                    } else {
                        i = 1;
                    }
                } else {
                    if (FUx >= -FUy) {
                        i = 7;
                    } else {
                        i = 6;
                    }
                }
            } else {
                if (FUy >= 0) {
                    if (-FUx >= FUy) {
                        i = 3;
                    } else {
                        i = 2;
                    }
                } else {
                    if (-FUx >= -FUy) {
                        i = 4;
                    } else {
                        i = 5;
                    }
                }
            }

            var FU = [FUx, FUy];

            var fU = norm(FU);
            var du = void 0;
            if (fU <= MU[i]) {
                du = FU;
            } else {
                du = mul(MU[i] / fU, FU);
            }

            //if (ui == 9) {
            //    console.log(ui, i, du, FU, MU[i]);
            //}

            //console.log(this.verts[ui], du, FU, MU);
            this.verts[ui][0] += du[0];
            this.verts[ui][1] += du[1];
            //console.log(FU, MU[i], du, this.verts[ui]);
        }
    }, {
        key: 'update',
        value: function update() {
            //console.log(this.adjMap);
            var FX = zeros(this.verts.length);
            var FY = zeros(this.verts.length);
            var M = [];
            for (var _i5 = 0; _i5 < this.verts.length; _i5++) {
                M.push([this.dbar, this.dbar, this.dbar, this.dbar, this.dbar, this.dbar, this.dbar, this.dbar]);
            }

            var barycenter = mul(1 / this.verts.length, sum(this.verts, 2));

            for (var ui = 0; ui < this.verts.length; ui++) {
                // Calculate gravity force
                var db = sub(barycenter, this.verts[ui]);
                var nDb = norm(db);
                FX[ui] += db[0] / nDb;
                FY[ui] += db[1] / nDb;

                // Calculate total node-node repulsive force
                for (var vi = 0; vi < this.verts.length; vi++) {
                    if (ui != vi) {
                        if (this.distance(ui, vi) >= 3 * this.delta) {
                            continue;
                        }

                        if (this.adjMap[ui].length == 2) {
                            if (this.adjMap[ui].includes(vi)) {
                                continue;
                            }
                        }

                        var F = this.forceRvert(this.verts[ui], this.verts[vi]);
                        //console.log("Fnnr", F);
                        //console.log(ui, vi, this.verts[ui], this.verts[vi], "Fnnr", F);
                        if (!isNaN(F[0])) {
                            FX[ui] += F[0];
                            FY[ui] += F[1];
                        }
                    }
                }

                // calculate edge attractive force
                var _iteratorNormalCompletion35 = true;
                var _didIteratorError35 = false;
                var _iteratorError35 = undefined;

                try {
                    for (var _iterator35 = this.adjMap[ui][Symbol.iterator](), _step35; !(_iteratorNormalCompletion35 = (_step35 = _iterator35.next()).done); _iteratorNormalCompletion35 = true) {
                        var _vi = _step35.value;

                        var _F = this.forceAvert(this.verts[ui], this.verts[_vi]);

                        FX[ui] += _F[0];
                        FY[ui] += _F[1];
                    }

                    // calculate node-edge repulsive force
                } catch (err) {
                    _didIteratorError35 = true;
                    _iteratorError35 = err;
                } finally {
                    try {
                        if (!_iteratorNormalCompletion35 && _iterator35.return) {
                            _iterator35.return();
                        }
                    } finally {
                        if (_didIteratorError35) {
                            throw _iteratorError35;
                        }
                    }
                }

                var _iteratorNormalCompletion36 = true;
                var _didIteratorError36 = false;
                var _iteratorError36 = undefined;

                try {
                    for (var _iterator36 = this.surrEdges[ui][Symbol.iterator](), _step36; !(_iteratorNormalCompletion36 = (_step36 = _iterator36.next()).done); _iteratorNormalCompletion36 = true) {
                        var edge = _step36.value;

                        var _edge3 = _slicedToArray(edge, 2),
                            ai = _edge3[0],
                            bi = _edge3[1];

                        if (ui == ai || ui == bi) {
                            continue;
                        }
                        var ve = this.computeVe(this.verts[ui], this.verts[ai], this.verts[bi]);

                        if (this.veOnEdge(ve, this.verts[ai], this.verts[bi])) {
                            var _F2 = this.forceRedge(this.verts[ui], this.verts[ai], this.verts[bi], ve);
                            if (!isNaN(_F2[0])) {
                                FX[ui] += _F2[0];
                                FY[ui] += _F2[1];
                            }
                        }
                    }
                } catch (err) {
                    _didIteratorError36 = true;
                    _iteratorError36 = err;
                } finally {
                    try {
                        if (!_iteratorNormalCompletion36 && _iterator36.return) {
                            _iterator36.return();
                        }
                    } finally {
                        if (_didIteratorError36) {
                            throw _iteratorError36;
                        }
                    }
                }

                var MU = M[ui];
                //console.log("Surr:", this.surrEdges);

                var _iteratorNormalCompletion37 = true;
                var _didIteratorError37 = false;
                var _iteratorError37 = undefined;

                try {
                    for (var _iterator37 = this.surrEdges[ui][Symbol.iterator](), _step37; !(_iteratorNormalCompletion37 = (_step37 = _iterator37.next()).done); _iteratorNormalCompletion37 = true) {
                        var _edge2 = _step37.value;

                        var _edge4 = _slicedToArray(_edge2, 2),
                            ai = _edge4[0],
                            bi = _edge4[1];

                        if (ui == ai || ui == bi) {
                            continue;
                        }
                        var ve = this.computeVe(this.verts[ui], this.verts[ai], this.verts[bi]);

                        var cv = void 0;

                        if (ui == 0 && ai == 5 && bi == 2) {}
                        //console.log(this.verts[ai], this.verts[bi], ve);

                        //console.log("v-e", ui, ai, bi);
                        if (this.veOnEdge(ve, this.verts[ai], this.verts[bi])) {
                            cv = sub(ve, this.verts[ui]);

                            var _i6 = void 0;
                            if (cv[0] >= 0) {
                                if (cv[1] >= 0) {
                                    if (cv[0] >= cv[1]) {
                                        _i6 = 0;
                                    } else {
                                        _i6 = 1;
                                    }
                                } else {
                                    if (cv[0] >= -cv[1]) {
                                        _i6 = 7;
                                    } else {
                                        _i6 = 6;
                                    }
                                }
                            } else {
                                if (cv[1] >= 0) {
                                    if (-cv[0] >= cv[1]) {
                                        _i6 = 3;
                                    } else {
                                        _i6 = 2;
                                    }
                                } else {
                                    if (-cv[0] >= -cv[1]) {
                                        _i6 = 4;
                                    } else {
                                        _i6 = 5;
                                    }
                                }
                            }

                            var maxR = norm(cv) / 2.1;
                            //console.log("???", cv);

                            //console.log(MU, maxR, Math.cos(Math.atan2(cv[1], cv[0])));
                            var ell = (_i6 + 4) % 8;
                            for (var j = 0; j < MU.length; j++) {
                                if ((_i6 - j + 8) % 8 == 0) {
                                    MU[j] = min(MU[j], maxR);
                                } else if ((_i6 - j + 8) % 8 == 1 || (_i6 - j + 8) % 8 == 2) {
                                    MU[j] = min(MU[j], maxR / Math.cos(Math.atan2(cv[1], cv[0]) - (j + 1) * Math.PI / 4));
                                } else if ((_i6 - j + 8) % 8 == 6 || (_i6 - j + 8) % 8 == 7) {
                                    MU[j] = min(MU[j], maxR / Math.cos(Math.atan2(cv[1], cv[0]) - j * Math.PI / 4));
                                }
                            }

                            for (var _j = 0; _j < MU.length; _j++) {
                                if ((ell - _j + 8) % 8 == 0) {
                                    M[ai][_j] = min(M[ai][_j], maxR);
                                } else if ((ell - _j + 8) % 8 == 1 || (ell - _j + 8) % 8 == 2) {
                                    M[ai][_j] = min(M[ai][_j], maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - (_j + 1) * Math.PI / 4));
                                } else if ((ell - _j + 8) % 8 == 6 || (ell - _j + 8) % 8 == 7) {
                                    M[ai][_j] = min(M[ai][_j], maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - _j * Math.PI / 4));
                                }
                            }

                            for (var _j2 = 0; _j2 < MU.length; _j2++) {
                                if ((ell - _j2 + 8) % 8 == 0) {
                                    M[bi][_j2] = min(M[bi][_j2], maxR);
                                } else if ((ell - _j2 + 8) % 8 == 1 || (ell - _j2 + 8) % 8 == 2) {
                                    M[bi][_j2] = min(M[bi][_j2], maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - (_j2 + 1) * Math.PI / 4));
                                } else if ((ell - _j2 + 8) % 8 == 6 || (ell - _j2 + 8) % 8 == 7) {
                                    M[bi][_j2] = min(M[bi][_j2], maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - _j2 * Math.PI / 4));
                                }
                            }
                        } else {
                            var va = sub(this.verts[ai], this.verts[ui]);
                            var vb = sub(this.verts[bi], this.verts[ui]);
                            if (norm(va) < norm(vb)) {
                                cv = va;
                            } else {
                                cv = vb;
                            }

                            var _i7 = void 0;
                            if (cv[0] >= 0) {
                                if (cv[1] >= 0) {
                                    if (cv[0] >= cv[1]) {
                                        _i7 = 0;
                                    } else {
                                        _i7 = 1;
                                    }
                                } else {
                                    if (cv[0] >= -cv[1]) {
                                        _i7 = 7;
                                    } else {
                                        _i7 = 6;
                                    }
                                }
                            } else {
                                if (cv[1] >= 0) {
                                    if (-cv[0] >= cv[1]) {
                                        _i7 = 3;
                                    } else {
                                        _i7 = 2;
                                    }
                                } else {
                                    if (-cv[0] >= -cv[1]) {
                                        _i7 = 4;
                                    } else {
                                        _i7 = 5;
                                    }
                                }
                            }

                            var _maxR = norm(cv) / 2.1;
                            //console.log("???", cv);

                            //console.log(MU, maxR, Math.cos(Math.atan2(cv[1], cv[0])));
                            var _ell = (_i7 + 4) % 8;
                            for (var _j3 = 0; _j3 < MU.length; _j3++) {
                                if ((_i7 - _j3 + 8) % 8 == 0) {
                                    MU[_j3] = min(MU[_j3], _maxR);
                                } else if ((_i7 - _j3 + 8) % 8 == 1 || (_i7 - _j3 + 8) % 8 == 2) {
                                    MU[_j3] = min(MU[_j3], _maxR / Math.cos(Math.atan2(cv[1], cv[0]) - (_j3 + 1) * Math.PI / 4));
                                } else if ((_i7 - _j3 + 8) % 8 == 6 || (_i7 - _j3 + 8) % 8 == 7) {
                                    MU[_j3] = min(MU[_j3], _maxR / Math.cos(Math.atan2(cv[1], cv[0]) - _j3 * Math.PI / 4));
                                }
                            }

                            var m = cv[1] / cv[0]; // Slope of cv
                            var n = -1 / m; // Slope of l

                            for (var _j4 = 0; _j4 < MU.length; _j4++) {
                                if ((_ell - _j4 + 8) % 8 == 0) {
                                    M[ai][_j4] = min(M[ai][_j4], _maxR);
                                } else if ((_ell - _j4 + 8) % 8 == 1 || (_ell - _j4 + 8) % 8 == 2) {
                                    M[ai][_j4] = min(M[ai][_j4], _maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - (_j4 + 1) * Math.PI / 4));
                                } else if ((_ell - _j4 + 8) % 8 == 6 || (_ell - _j4 + 8) % 8 == 7) {
                                    M[ai][_j4] = min(M[ai][_j4], _maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - _j4 * Math.PI / 4));
                                }
                            }

                            for (var _j5 = 0; _j5 < MU.length; _j5++) {
                                if ((_ell - _j5 + 8) % 8 == 0) {
                                    M[bi][_j5] = min(M[bi][_j5], _maxR);
                                } else if ((_ell - _j5 + 8) % 8 == 1 || (_ell - _j5 + 8) % 8 == 2) {
                                    M[bi][_j5] = min(M[bi][_j5], _maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - (_j5 + 1) * Math.PI / 4));
                                } else if ((_ell - _j5 + 8) % 8 == 6 || (_ell - _j5 + 8) % 8 == 7) {
                                    M[bi][_j5] = min(M[bi][_j5], _maxR / Math.cos(Math.atan2(-cv[1], -cv[0]) - _j5 * Math.PI / 4));
                                }
                            }
                        }
                    }
                    //if (ui == 0) console.log("MU0", MU, FX[ui], FY[ui]);
                } catch (err) {
                    _didIteratorError37 = true;
                    _iteratorError37 = err;
                } finally {
                    try {
                        if (!_iteratorNormalCompletion37 && _iterator37.return) {
                            _iterator37.return();
                        }
                    } finally {
                        if (_didIteratorError37) {
                            throw _iteratorError37;
                        }
                    }
                }
            }

            //console.log("Fx", FX);
            for (var _ui in this.verts) {

                this.move(_ui, FX[_ui], FY[_ui], M[_ui]);
            }
        }
    }]);

    return ForceLinkDiagram;
}();

function sleep(millis) {
    var date = new Date();
    var curDate = null;
    do {
        curDate = new Date();
    } while (curDate - date < millis);
}

var workerFunctions = {
    setLinkDiagram: function setLinkDiagram(sigma, crossBend) {
        self.shadow = new _shadow2.default(sigma);
        self.orthshadow = new OrthogonalDiagramEmbedding(self.shadow);

        workerFunctions.embedDiagram();
    },

    embedDiagram: function embedDiagram() {
        var tstart = Date.now();

        var thresh = 5e-10;
    }
};

onmessage = function onmessage(e) {
    workerFunctions[e.data.function].apply(workerFunctions, _toConsumableArray(e.data.arguments));
};
//# sourceMappingURL=orth_worker.js.map