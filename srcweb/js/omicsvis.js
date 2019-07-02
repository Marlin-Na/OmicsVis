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
            .display(tnt.board.track.feature.axis().orientation("top"));

        let contig_track = tnt.board.track()
            .id("contig")
            .height(10)
            .data(tnt.board.track.data.sync().retriever(() => []))
            .display(tnt.board.track.feature.block().color("#082A46"));

        let gene_track = tnt.board.track()
            .id("gene")
            .height(20)
            .data(tnt.board.track.data.sync().retriever(() => []))
            .display(
                tnt.board.track.feature.block().color("#AD9274")
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

        let diamond_track = tnt.board.track()
            .id("diamond")
            .height(60)
            .data(tnt.board.track.data.sync().retriever(() => []))
            .display(
                tnt.board.track.feature.genome.gene().color("green")
                    .on("click", function(d) {
                        console.log(d);

                        // ad hoc
                        let mapped_to = d.eggnog_pos_X2.split(".");
                        let taxonomy = mapped_to[0];
                        let gene = mapped_to[1];
                        taxonomy = `<a target="_blank" href="https://www.uniprot.org/taxonomy/${taxonomy}">${taxonomy}</a>`
                        gene = `<a target="_blank" href="https://www.uniprot.org/uniprot/?query=+${gene}">${gene}</a>`
                        
                        tnt.tooltip.table()
                            .width(300)
                            .call(this, {
                                header: `${d.gene_ID} Mapping`,
                                rows: [
                                    {label: "Taxonomy", value: taxonomy},
                                    {label: "Gene", value: gene},
                                ]
                            })
                    })
            );

        // Initialize the board
        this.board(this.dom);

        if (!concise)
            this.board.add_track(axis_track);

        this.board
            .add_track(contig_track)
            .add_track(gene_track)
            .add_track(diamond_track);

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
        let diamond_track = board.find_track("diamond");
        let _this = this;

        if (this.data_src === null) {
            console.error("data_src is not set");
            return this;
        }

        // Attaching callbacks to the promise
        this.data_src = this.data_src
            .then(data => {
                let seqlen = data.meta.seqlen;
                let seqname = data.meta.seqname;
                let gtrack_data = data.gene_track;
                let dtrack_data = data.diamond_track;

                board.from(0).to(seqlen+1)
                    .max(seqlen+1)
                    .zoom_out(seqlen+1); // Maximum extent

                let gene_data_retriever = tnt.board.track.data.sync()
                    .retriever(function(loc) {
                        // We are using start and end
                        gtrack_data.forEach(e => {
                            e.start = e.gene_start;
                            e.end = e.gene_end;
                        });
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

                let diamond_data_retriever = tnt.board.track.data.sync()
                    .retriever(function(loc) {
                        dtrack_data.forEach(e => {
                            e.start = e.eggnog_pos_start;
                            e.end = e.eggnog_pos_end;
                            e.id = e.gene_ID;
                            e.display_label = e.gene_ID;
                        });
                        return dtrack_data;
                    });
                diamond_track.data(diamond_data_retriever);

                board.start();
                return data;
            })
            .catch(err => console.log(err));
        return this;
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

let gridOptions;
let option_filterSelected = false;

// Header filter handler
function binding_filterSelected(node) {
    option_filterSelected = node.checked;
    gridOptions.api.onFilterChanged();
}

async function load_data_table() {
    // Alternative to load_contig_list
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
    let rowData = res;

    let ActiveViews = new Map();

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
            ActiveViews.set(contig_id, view);
        }
        else {
            document.getElementById("vis-" + contig_id).remove();
            ActiveViews.delete(contig_id);
        }
    }
    function isExternalFilterPresent() {
        if (option_filterSelected === false)
            return false;
        else
            return true;
    }
    function doesExternalFilterPass(node) {
        if (option_filterSelected)
            return node.isSelected();
    }
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

    gridOptions = {
        columnDefs: columnDefs,
        rowData: rowData,
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
    let dom_table = document.getElementById("contig-table");
    let table = new agGrid.Grid(dom_table, gridOptions);
}


async function main() {
    await load_data_table();
}
main();

