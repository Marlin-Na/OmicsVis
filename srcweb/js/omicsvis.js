'use strict';

// Assuming global tnt and d3v5 variable.

class TrackView {

    constructor(div) {
        this.dom = div;
        this.opt = {
            width: 800,
        };
        this.board = null;
        this.data_src = null;
        this.data = null;
    }

    init_vis() {
        // width is determined from the width of div
        // let width = d3v5.select(this.dom).node().getBoundingClientRect().width/1.1;
        let opt = this.opt;
        let width = opt.width; // Fixed width

        if (this.board !== null)
            console.error("This board has already been initialized.");
        this.board = tnt.board();
        this.board.allow_drag(false);
        this.board.from(0).to(1000).max(1000).width(width);

        let _this = this;

        // Initialize tracks

        let axis_track = tnt.board.track()
            .id("axis")
            .height(5)
            .display(tnt.board.track.feature.axis().orientation("top"));

        let contig_track = tnt.board.track()
            .id("contig")
            .height(10)
            .display(tnt.board.track.feature.block().color("#082A46"));

        let gene_track = tnt.board.track()
            .id("gene")
            .height(20)
            .color("white")
            .display(
                tnt.board.track.feature.genome.gene().color("#AD9274") // default color
                    .on("click", function(d) {
                        console.log(d);
                        tnt.tooltip.table()
                            .width(300)
                            .call(this, {
                                header: `Gene: ${d.gene_ID}`,
                                rows: [
                                    {label: "Gene Strand", value: d.gene_strand},
                                    {label: "RBS Motif", value: d.gene_rbs_motif},
                                    {label: "Type", value: d.eggnog_X13},
                                ]
                            });
                    })
            );
        gene_track.display().layout()
            .fixed_slot_type("collapsed") //.fixed_slot_type("expanded")
            .keep_slots(false)
            .on_layout_run(function(types, current) {
                // With some extra space
                let needed_height = types.collapsed.needed_slots * types.collapsed.slot_height + 20;
                if (needed_height !== gene_track.height()) {
                    gene_track.height(needed_height);
                    _this.board.tracks(_this.board.tracks());
                }
            });

        let diamond_track = tnt.board.track()
            .id("diamond")
            .height(60)
            .color("white")
            .display(
                tnt.board.track.feature.genome.gene().color("green") // default color
                    .on("click", function(d) {
                        console.log(d);

                        // ad hoc
                        let mapped_to = d.eggnog_pos_X2.split(".");
                        let taxonomy = mapped_to[0];
                        let gene = mapped_to[1];
                        taxonomy = `<a target="_blank" href="https://www.uniprot.org/taxonomy/${taxonomy}">${taxonomy}</a>`;
                        gene = `<a target="_blank" href="https://www.uniprot.org/uniprot/?query=+${gene}">${gene}</a>`;

                        tnt.tooltip.table()
                            .width(300)
                            .call(this, {
                                header: `${d.gene_ID} Mapping`,
                                rows: [
                                    {label: "Taxonomy", value: taxonomy},
                                    {label: "Gene", value: gene},
                                ]
                            });
                    })
            );
        // Refer http://bl.ocks.org/emepyc/7c73519ee7a1300eb68a
        diamond_track.display().layout()
            .fixed_slot_type("expanded")
            .keep_slots(false)
            .on_layout_run(function(types, current) {
                let needed_height = types.expanded.needed_slots * types.expanded.slot_height;
                if (needed_height !== diamond_track.height()) {
                    diamond_track.height(needed_height);
                    //genome.tracks(genome.tracks());
                    _this.board.tracks(_this.board.tracks());
                }
            });

        // Initialize the board
        this.board(this.dom);

        this.board.add_track(axis_track);

        this.board
            .add_track(contig_track)
            .add_track(gene_track)
            .add_track(diamond_track);

        // Do not start
        //this.board.start();
        return this;
    }
    // resize_vis() {
    //     let width = d3v5.select(this.dom).node().getBoundingClientRect().width/1.1;
    //     this.board.width(width);
    // }

    async fetch_data_src() {
        if (this.data_src === null) {
            console.error("data_src is not set");
            return;
        }
        if (this.data !== null) {
            // i.e. we should assume it has been already fetched
            return;
        }
        let data = await this.data_src;
        this.data = data;
        return;
    }

    async update_vis() {
        if (this.board === null)
            console.error("The board has not been initialized");

        await this.fetch_data_src();
        let data = this.data;
        let board = this.board;

        let contig_track = board.find_track("contig");
        let gene_track = board.find_track("gene");
        let diamond_track = board.find_track("diamond");

        let _this = this;

        // Modify the board range
        board
            .from(-1)
            .to(data.meta.seqlen+1);
        
        // Add contig track data
        contig_track.data(tnt.board.track.data.sync().retriever(
            function(loc) {
                return [{
                    start: 1, end: data.meta.seqlen
                }];
            }
        ));

        // Add gene track data
        gene_track.data(tnt.board.track.data.sync().retriever(
            function(loc) {
                let gtrack_data = data.gene_track;
                gtrack_data.forEach(e => {
                    e.start = e.gene_start;
                    e.end = e.gene_end;
                    e.id = e.gene_ID;
                    e.display_label = "";
                });
                return gtrack_data;
            }
        ));

        // Add diamond track data
        diamond_track.data(tnt.board.track.data.sync().retriever(
            function(loc) {
                let dtrack_data = data.diamond_track;
                dtrack_data.forEach(e => {
                    e.start = e.eggnog_pos_start;
                    e.end = e.eggnog_pos_end;
                    e.id = e.gene_ID;
                    e.display_label = e.gene_ID;
                    //e.display_label = "";
                    //e.color = "green";
                })
                return dtrack_data;
            }
        ));

        board.start();

        // tmp
        this.set_color_by("strand");
        board.start();
    }

    set_color_by(what) {
        let _this = this;
        if (this.data === null) {
            console.error("The board is not yet loaded");
            return;
        }
        if (what === "strand") {
            let scale = d3v5.scaleOrdinal();
            scale.domain(["-", "+"]);
            scale.range(d3v5.schemeAccent.slice(6, 8));
            this.data.gene_track.forEach(e => {
                e.color = scale(e.gene_strand);
            });
        }
    }

    set_data_src(name) {
        // such as MAG06_80
        let data_src = fetch("sample_data/" + name + ".json")
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

class TrackViewPanel {
    constructor(div) {
        this.dom = div;
        this.ActiveViews = new Map();
        window.ActiveViews = this.ActiveViews;
    }
    add_view(contig_id) {
        let vis_dom = d3v5.select(this.dom).node();
        let container = document.createElement("div");
        container.id = "vis-" + contig_id;
        vis_dom.appendChild(container);
        let view = new TrackView(container);
        view.set_data_src(contig_id);
        view.init_vis();
        this.ActiveViews.set(contig_id, view);
        view.update_vis();
    }
    remove_view(contig_id) {
        document.getElementById("vis-" + contig_id).remove();
        ActiveViews.delete(contig_id);
    }
}

class IndexTable {
    constructor(div) {
        this.dom = div;
        this.gridOptions = null;
        this.option_filterSelected = false;

        // TODO: move outside of this class.
        this.vispanel = new TrackViewPanel("#vis");
    }
    async load() {
        let _this = this;
        let res;
        res = await fetch("sample_data/index.json");
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
                checkboxSelection: true,
                //headerCheckboxSelection: true,
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

        this.gridOptions = {
            columnDefs: columnDefs,
            rowData: res,
            pagination: true,
            paginationAutoPageSize: true,
            rowSelection: "multiple",
            rowMultiSelectWithClick: true,
            onRowSelected: onRowSelected,
            isExternalFilterPresent: isExternalFilterPresent,
            doesExternalFilterPass: doesExternalFilterPass,
            onCellMouseOver: onCellMouseOver,
            onCellMouseOut: onCellMouseOut
        };




        let dom_table = d3v5.select(this.dom).node();
        let table = new agGrid.Grid(dom_table, this.gridOptions);


        function onRowSelected(event) {
            let contig_id = event.data.contig;
            let is_checked = event.node.isSelected();

            if (is_checked) {
                _this.vispanel.add_view(contig_id);
            }
            else {
                _this.vispanel.remove_view(contig_id);
            }
        }
        function isExternalFilterPresent() {
            if (_this.option_filterSelected === false)
                return false;
            else
                return true;
        }
        function doesExternalFilterPass(node) {
            if (_this.option_filterSelected)
                return node.isSelected();
        }
        // TODO: move to the TrackViewPanel class
        function onCellMouseOver(event) {
            if (event.node.selected) {
                let contig_id = event.data.contig;
                let vis_dom = document.getElementById("vis-" + contig_id);
                let the_board = ActiveViews.get(contig_id).board;
                d3v5.select(vis_dom).classed("tntboard-highlight", true);
            }
        }
        function onCellMouseOut(event) {
            let contig_id = event.data.contig;
            if (ActiveViews.has(contig_id)) {
                let vis_dom = document.getElementById("vis-" + contig_id);
                let the_board = ActiveViews.get(contig_id).board;
                d3v5.select(vis_dom).classed("tntboard-highlight", false);
            }
        }

    }
}

let table;

async function main() {
    table = new IndexTable("#contig-table");
    await table.load();
}
main();

// Header filter handler
function binding_filterSelected(node) {
    table.option_filterSelected = node.checked;
    table.gridOptions.api.onFilterChanged();
}


