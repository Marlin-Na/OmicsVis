'use strict';

// Assuming global tnt and d3v5 variable.
class TrackView {
    constructor(div) {
        this.dom = div;
        this.board = tnt.board();
        this.board.allow_drag(false);
        this.data_src = null;
    }
    init_vis(concise = false) {
        // width is determined from the width of div
        // let width = d3v5.select(this.dom).node().getBoundingClientRect().width/1.1;
        let width = concise? 300 : 800; // Fixed width
        this.board.from(0).to(1000).max(1000).width(width);

        // Initialize tracks

        let axis_track = tnt.board.track()
            .id("axis")
            .height(5)
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
            .height(20)
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

        if (!concise)
            this.board.add_track(axis_track);

        this.board
            .add_track(contig_track)
            .add_track(gene_track);

        this.board.start();
        return this;
    }
    // resize_vis() {
    //     let width = d3v5.select(this.dom).node().getBoundingClientRect().width/1.1;
    //     this.board.width(width);
    // }
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

async function load_data_table() {
    // Alternative to load_contig_list
    let res;
    res = await fetch("/sample_data/index.json");
    if (res.ok)
        res = await res.json();
    else
        throw "Network error";



    let colorScale_contiglength = d3v5.scaleSequential(d3v5.interpolateReds)
        .domain(res.map(d => d.length));
    let colorScale_NGenes = d3v5.scaleSequential(d3v5.interpolateReds)
        .domain(res.map(d => d.number_gene));

    let columnDefs = [
        {
            headerName: "Contig ID",
            field: "contig",
            sortable: true,
        },
        {
            headerName: "Len",
            headerTooltip: "Contig Length",
            field: "length",
            filter: "agNumberColumnFilter",
            type: "numericColumn",
            sortable: true,
            width: 100,
            cellStyle: function(params) {
                let color = colorScale_contiglength(params.value);
                return {backgroundColor: color, color: "grey"};
            }
        },
        {
            headerName: "NGenes",
            headerTooltip: "Number of Genes",
            field: "number_gene",
            filter: "agNumberColumnFilter",
            type: "numericColumn",
            sortable: true,
            width: 100,
            cellStyle: function(params) {
                let color = colorScale_NGenes(params.value);
                return {backgroundColor: color, color: "grey"};
            }
        },
    ];
    let rowData = res;

    function onRowSelected(event) {
        let contig_id = event.data.contig;
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
        pagination: true,
        paginationAutoPageSize: true,
        rowSelection: "multiple",
        rowMultiSelectWithClick: true,
        onRowSelected: onRowSelected
    };
    let dom_table = document.getElementById("contig-table");
    let table = new agGrid.Grid(dom_table, gridOptions);
}


async function main() {
    await load_data_table();
}
main();

