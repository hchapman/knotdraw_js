class Arc {
    constructor(index) {
        this.index = index;
    }

    toString() {
        return this.index.toString();
    }

    valueOf() {
        return this.index;
    }
}

export default class LinkShadow {
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
            this.newArc(2*ei);
            this.newArc(2*ei+1);
            //console.log(this.arcs)

            this.setEdge(ei, [2*ei, 2*ei+1]);
        }

        // Set verts by looking through verts
        for (let vi in verts) {
            this.setVert(vi, verts[vi]);
        }

        this.generateFaces();
        this.generateComponents();
        this.connectArcs();
    }

    copy() {
        let link = new LinkShadow(this.verts.map(
            v => v.map(a => a.index)
        ));

        return link;
    }

    generateFaces() {
        let leftArcs = new Set(this.arcs);
        this.faces = [];

        while(leftArcs.size > 0) {
            let startArc = Array.from(leftArcs).pop();

            let face = [];
            let arc = startArc;
            let _failsafe = 0;
            do {
                leftArcs.delete(arc);
                face.push(arc);
                arc.face = this.faces.length;
                //console.log(arc);

                let oArc = this.edges[arc.edge][(arc.edgepos+1)%2];
                arc = this.verts[oArc.vert][(oArc.vertpos+1)%4];
                //console.log(arc);
                _failsafe += 1;
                if (_failsafe > 500) {
                    console.log("Failure");
                    return this.faces;
                }
            } while (arc != startArc)

            //face.reverse();
            this.faces.push(face);
        }
        return this.faces;
    }

    generateComponents(oneOrient=true) {
        let leftArcs = new Set(this.arcs);

        this.components = [];
        while (leftArcs.size > 0) {
            let startArc = Array.from(leftArcs).pop();

            let component = this.component(startArc);
            for (let arc of component) {
                leftArcs.delete(arc);
                if (oneOrient) {
                    // Delete the other arc edge-opposite this one
                    leftArcs.delete(this.edges[arc.edge][(arc.edgepos+1)%2]);
                }
            }

            this.components.push(component);
        }
        return this.components;
    }

    component(arc) {
        let startArc = arc;

        let component = [];
        do {
            component.push(arc);

            let oArc = this.edges[arc.edge][(arc.edgepos+1)%2];
            arc = this.verts[oArc.vert][(oArc.vertpos+2)%4];
        } while (arc != startArc)

        return component;
    }

    outVertI(arc) {
        return arc.vert;
    }

    inVertI(arc) {
        return this.edgeOpposite(arc).vert;
    }

    edgeOpposite(arc) {
        return this.edges[arc.edge][(arc.edgepos+1)%2];
    }
    vertNext(arc) {
        return this.verts[arc.vert][(arc.vertpos+1)%4];
    }
    vertPrev(arc) {
        return this.verts[arc.vert][(arc.vertpos+3)%4];
    }

    newArc(idx) {
        if (idx === undefined) { idx = this.arcs.length; }

        this.arcs[idx] = new Arc(idx);
        // {index: idx, edge:undefined, edgepos:undefined, vert:undefined, vertpos:undefined};
        return this.arcs[idx];
    }

    connectArcs() {
        /* Hook up the arcs so they know, locally, their linkings */

        for (let arc of this.arcs) {
            arc.edgeOpposite = this.edgeOpposite(arc);
            arc.vertNext = this.vertNext(arc);
            arc.vertPrev = this.vertPrev(arc);
        }
    }

    setEdge(idx, ais) {
        this.edges[idx] = ais.map((ai) => {return this.arcs[ai];}, this);

        for (let i in ais) {
            this.arcs[ais[i]].edge = idx;
            this.arcs[ais[i]].edgepos = parseInt(i);
        }
    }

    setVert(idx, ais) {
        this.verts[idx] = ais.map((ai) => {return this.arcs[ai];}, this);

        for (let i in ais) {
            if (ais[i] === undefined) { continue; }

            this.arcs[ais[i]].vert = parseInt(idx);
            this.arcs[ais[i]].vertpos = parseInt(i);
        }
    }
}
