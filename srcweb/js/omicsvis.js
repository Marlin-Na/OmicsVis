'use strict';

// Assuming global tnt and d3v5 variable.
class TrackView {
    constructor(div) {
        this.dom = div;
        this.board = tnt.board();
        this.board.allow_drag(false);
        this.data_src = null;
    }
    init_vis() {
        // width is determined from the width of div
        let width = d3v5.select(this.dom).node().getBoundingClientRect().width/1.1;
        this.board.from(0).to(1000).max(1000).width(width);

        // Initialize tracks

        let axis_track = tnt.board.track()
            .id("axis")
            .height(20)
            .color("white")
            .display(tnt.board.track.feature.axis().orientation("top"));

        let contig_track = tnt.board.track()
            .id("contig")
            .height(10)
            .color("white")
            .data(tnt.board.track.data.sync().retriever(() => []))
            .display(tnt.board.track.feature.block().color("#082A46"));

        let gene_track = tnt.board.track()
            .id("gene")
            .height(30)
            .color("white")
            .data(tnt.board.track.data.sync().retriever(() => []))
            .display(
                tnt.board.track.feature.block().color("#AD9274")
                    .on("click", function(d) {
                        tnt.tooltip.table()
                            .width(300)
                            .call(this, {
                                header: "EGGNOG Annotation",
                                rows: [
                                    {"label": "RBS Motif", "value": d.rbs_motif},
                                    {"label": "Type", "value": d.eggnog.X13},
                                ]
                            });
                    })
            );

        // Initialize the board
        this.board(this.dom);
        this.board
            .add_track(axis_track)
            .add_track(contig_track)
            .add_track(gene_track);

        this.board.start();
        return this;
    }
    update_vis() {
        let board = this.board;
        let contig_track = board.find_track("contig");
        let gene_track = board.find_track("gene");
        let _this = this;

        if (this.data_src === null) {
            console.error("data_src is not set");
            return this;
        }

        // Attaching callbacks to the promise
        this.data_src = this.data_src
            .then(data => {
                let seqlen = data.seqlen;
                let seqname = data.seqname;
                let gtrack_data = data.gene_track;

                board.from(0).to(seqlen+1)
                    .max(seqlen+1)
                    .zoom_out(seqlen+1); // Maximum extent

                let gene_data_retriever = tnt.board.track.data.sync()
                    .retriever(function(loc) {
                        // We are using start and end
                        return gtrack_data;
                    });
                gene_track.data(gene_data_retriever);

                let contig_data_retriever = tnt.board.track.data.sync()
                    .retriever(function(loc) {
                        return [{
                            start: 0,
                            end: seqlen + 1
                        }];
                    });
                contig_track.data(contig_data_retriever);

                board.start();
                return data;
            })
            .catch(err => console.log(err));
        return this;
    }
    set_data_src(name) {
        // such as MAG06_80
        let data_src = fetch("/sample_data/" + name + ".json")
            .then(res => {
                if (!res.ok)
                    throw new Error("HTTP error: " + res.status);
                return res.json();
            })
            .catch(err => console.log(err));
        this.data_src = data_src;

        // TODO: call update_vis?
        return this;
    }
}

async function load_contig_list() {
    let res;
    res = await fetch("/sample_data/index.txt");
    if (res.ok)
        res = await res.text();
    else
        throw "Network error";
    res = res.split("\n").filter(d => ! (d === ""));
    let dom_contig_list = document.getElementById("contig-list");
    let dom_ul = document.createElement("ul");
    dom_contig_list.appendChild(dom_ul);
    function on_click() {
        let is_checked = this.checked;
        let contig_id = this.dataset.contigId;
        let vis_dom = document.getElementById("vis");
        if (is_checked) {
            let container = document.createElement("div");
            container.id = "vis-" + contig_id;
            vis_dom.appendChild(container);
            let view = new TrackView(container);
            view.set_data_src(contig_id);
            view.init_vis();
            view.update_vis();
        }
        else
            document.getElementById("vis-" + contig_id).remove();
    }
    for (let each of res) {
        let el = document.createElement("li");
        let checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        checkbox.onclick = on_click;
        checkbox.dataset.contigId = each;
        el.appendChild(checkbox);
        el.appendChild(document.createTextNode(each));
        dom_ul.appendChild(el);
    }
}

async function load_data_table() {
    // Alternative to load_contig_list
    let res;
    res = await fetch("/sample_data/index.txt");
    if (res.ok)
        res = await res.text();
    else
        throw "Network error";
    res = res.split("\n").filter(d => ! (d === ""));
    // Convert to json
    res = res.map(d => {return {contig_id: d}});

    let columnDefs = [
        {headerName: "Contig ID", field: "contig_id"}
    ];
    let rowData = res;

    function onRowSelected(event) {
        console.log(event);
        window.event = event;
        let contig_id = event.data.contig_id;
        let is_checked = event.node.isSelected();

        let vis_dom = document.getElementById("vis");
        if (is_checked) {
            let container = document.createElement("div");
            container.id = "vis-" + contig_id;
            vis_dom.appendChild(container);
            let view = new TrackView(container);
            view.set_data_src(contig_id);
            view.init_vis();
            view.update_vis();
        }
        else
            document.getElementById("vis-" + contig_id).remove();
    }

    let gridOptions = {
        columnDefs: columnDefs,
        rowData: rowData,
        rowSelection: "multiple",
        rowMultiSelectWithClick: true,
        onRowSelected: onRowSelected
    };
    let dom_table = document.getElementById("contig-table");
    let table = new agGrid.Grid(dom_table, gridOptions);
}


async function main() {
    //await load_contig_list();
    await load_data_table();
}
main();

