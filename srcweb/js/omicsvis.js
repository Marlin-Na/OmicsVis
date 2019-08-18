'use strict';

// ad hoc configuration
const INDEX_PATH = "/data/top_abundance_index.json";
const CONTIG_DATA_PATH  = "/data/contigs";
const EXPERIMENT_INFO_PATH = "/data/experiment_info.json";
const EXPERIMENT_DATA_PATH = "/data/experiment";

// Currently this file assumes global tnt (by importing d3 version 3, tnt.genome, tnt.tooltip)
// and d3v5 variable are available.

// I intend to use TrackView class to represent a TnT Board instance for one contig.
// They are initialized and removed on the fly according to the "selected" states in
// the table.
class TrackView {

    constructor(div) {
        this.dom = div;
        this.opt = {
            width: 800,
        };
        // The tnt board instance
        this.board = null;
        // "data_src" is set by set_data_src method, it will replace null with a
        // promise that resolves to the actual data.
        this.data_src = null;
        // Place to store the actual data
        this.data = null;
        // // tnt actually does not accept a color callback for individual features (really?),
        // // thus we embed this callback to the data retriever.
        // this.gene_colorgen = function(d) {
        //     return undefined;
        // };
    }

    // init_vis should only called once.
    init_vis() {
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
                        let tooltip_data = d.tooltip_data;
                        let tooltip_rows = [];
                        for (let key of Object.keys(tooltip_data)) {
                            tooltip_rows.push({
                                label: key,
                                value: tooltip_data[key]
                            });
                        }
                        console.log(tooltip_data);
                        console.log(tooltip_rows);
                        tnt.tooltip.table()
                            .width(300)
                            .call(this, {
                                header: `Gene: ${tooltip_data.gene_id}`,
                                // rows: [
                                //     {label: "Gene Strand", value: d.gene_strand},
                                //     {label: "RBS Motif", value: d.gene_rbs_motif},
                                //     {label: "Type", value: d.eggnog_X13},
                                // ]
                                rows: tooltip_rows
                            });
                    })
            );
        // Set fixed "collapsed" layout and allow variable height
        // Refer http://bl.ocks.org/emepyc/7c73519ee7a1300eb68a
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

        // let diamond_track = tnt.board.track()
        //     .id("diamond")
        //     .height(60)
        //     .color("white")
        //     .display(
        //         tnt.board.track.feature.genome.gene().color("green") // default color
        //             .on("click", function(d) {
        //                 // Specify the tooltip with link to external website
        //                 console.log(d);
        //                 // ad hoc
        //                 let mapped_to = d.eggnog_pos_X2.split(".");
        //                 let taxonomy = mapped_to[0];
        //                 let gene = mapped_to[1];
        //                 taxonomy = `<a target="_blank" href="https://www.uniprot.org/taxonomy/${taxonomy}">${taxonomy}</a>`;
        //                 gene = `<a target="_blank" href="https://www.uniprot.org/uniprot/?query=+${gene}">${gene}</a>`;

        //                 tnt.tooltip.table()
        //                     .width(300)
        //                     .call(this, {
        //                         header: `${d.gene_ID} Mapping`,
        //                         rows: [
        //                             {label: "Taxonomy", value: taxonomy},
        //                             {label: "Gene", value: gene},
        //                         ]
        //                     });
        //             })
        //     );

        // diamond_track.display().layout()
        //     .fixed_slot_type("expanded")
        //     .keep_slots(false)
        //     .on_layout_run(function(types, current) {
        //         let needed_height = types.expanded.needed_slots * types.expanded.slot_height;
        //         if (needed_height !== diamond_track.height()) {
        //             diamond_track.height(needed_height);
        //             //genome.tracks(genome.tracks());
        //             _this.board.tracks(_this.board.tracks());
        //         }
        //     });

        // Bind thd dom
        this.board(this.dom);

        // Add tracks
        this.board.add_track(axis_track);

        this.board
            .add_track(contig_track)
            .add_track(gene_track);

        // Do not start
        //this.board.start();
        return this;
    }
    // resize_vis() {
    //     let width = d3v5.select(this.dom).node().getBoundingClientRect().width/1.1;
    //     this.board.width(width);
    // }

    // Set "data_src" as a promise that will be resolved to the actual data
    set_data_src(name) {
        // such as MAG06_80
        let data_src = fetch(CONTIG_DATA_PATH + "/" + name + ".json")
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

    // This will resolves the promise in "data_src" to the actual data and assign to "data"
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

    // It will modify the default track parameters with actual data.
    // Then it will start rendering the board.
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
            .to(data.contig_length+1);

        // Add contig track data
        contig_track.data(tnt.board.track.data.sync().retriever(
            function(loc) {
                return [{
                    start: 1, end: data.contig_length
                }];
            }
        ));

        // Add gene track data
        gene_track.data(tnt.board.track.data.sync().retriever(
            function(loc) {
                // let gtrack_data = data.genes;
                // gtrack_data.forEach(e => {
                //     e.start = e.gene_start;
                //     e.end = e.gene_end;
                //     e.id = e.gene_id;
                //     e.display_label = "";
                //     // Set color
                //     e.color = _this.gene_colorgen(e);
                // });
                // return gtrack_data;


                let gtrack_src_data = data.genes;
                let gtrack_data = gtrack_src_data.map(d => {
                    let ans = {};
                    ans.start = d.gene_start;
                    ans.end   = d.gene_end;
                    ans.id    = d.gene_id;
                    ans.display_label = "";
                    // Set color
                    //ans.color = _this.gene_colorgen(d);
                    ans.tooltip_data = d;
                    return ans;
                });
                return gtrack_data;
            }
        ));

        // // Add diamond track data
        // diamond_track.data(tnt.board.track.data.sync().retriever(
        //     function(loc) {
        //         let dtrack_data = data.diamond_track;
        //         dtrack_data.forEach(e => {
        //             e.start = e.eggnog_pos_start;
        //             e.end = e.eggnog_pos_end;
        //             e.id = e.gene_ID;
        //             e.display_label = e.gene_ID;
        //             //e.display_label = "";
        //             //e.color = "green";
        //         })
        //         return dtrack_data;
        //     }
        // ));

        board.start();
    }

    // "Restart" the board. Needed when "gene_colorgen" etc are modified
    reload() {
        if (this.board === null) {
            return; // do nothing
        }
        this.board.start();
    }

    // // Modify the color scale of features, i.e. by strand/metric
    // set_gene_colorby(what) {
    //     if (what === null || what === "null") {
    //         this.gene_colorgen = function(d) {
    //             return undefined;
    //         };
    //         this.reload();
    //         return;
    //     }
    //     if (what === "strand") {
    //         let scale = d3v5.scaleOrdinal();
    //         scale.domain(["-", "+"]);
    //         scale.range(d3v5.schemeAccent.slice(6, 8));
    //         this.gene_colorgen = function(d) {
    //             return scale(d.gene_strand);
    //         };
    //         this.reload();
    //         return;
    //     }
    //     // e.g. metric_1
    //     let scale = d3v5.scaleSequential(d3v5.interpolateBlues);
    //     scale.domain([0, 100]);
    //     this.gene_colorgen = function(d) {
    //         if (d[what] === undefined)
    //             return "grey";
    //         return scale(d[what]);
    //     };
    //     this.reload();
    // }

}

// I intend to use this class to handle the list of active tnt board instances
class TrackViewPanel {
    constructor(div) {
        this.dom = div;
        // Store the available tnt board instances with contig_id as the key
        this.ActiveViews = new Map();
        //this.opt_gene_colorby = null;
    }
    set_gene_colorby(what) {
        // Calls the set_gene_colorby method of TrackView
        this.opt_gene_colorby = what;
        this.ActiveViews.forEach(v => {
            v.set_gene_colorby(what);
        });
    }
    // "add_view" will append an div container inside the vis_dom and initialize
    // the tracks. "remove_view" will remove the dom corresponding to the contig_id.
    add_view(contig_id) {
        let vis_dom = d3v5.select(this.dom).node();

        let container_template = (contig_id) => {
            let container = document.createElement("div");
            container.id = "viscontainer-" + contig_id;
            container.className = "viscontainer";
            container.innerHTML = `
                 <label>${contig_id}</label>
            `;
            let tntdiv = document.createElement("div");
            tntdiv.id = "vis-" + contig_id;
            container.appendChild(tntdiv);
            return [container, tntdiv];
        };
        let [container, tntdiv] = container_template(contig_id);
        vis_dom.appendChild(container);
        let view = new TrackView(tntdiv);
        view.set_data_src(contig_id);
        // // Set gene color callback
        // if (this.opt_gene_colorby !== null)
        //     view.set_gene_colorby(this.opt_gene_colorby);
        view.init_vis();
        this.ActiveViews.set(contig_id, view);
        view.update_vis();
    }
    remove_view(contig_id) {
        document.getElementById("viscontainer-" + contig_id).remove();
        this.ActiveViews.delete(contig_id);
    }
}

class ExperimentPanel {
    constructor(div) {
        this.div = div;
        this.svg = null;
        this.experiment_info = fetch(EXPERIMENT_INFO_PATH)
            .then(res => {
                if (res.ok)
                    return res.json();
                throw new Error("Network response was not okay");
            })
            .catch(err => console.error(err));
        this.updateloop = new Promise(resolve => resolve("ok"));
        this._init();
        this.data = [];
    }
    _init() {
        d3v5.select(this.div)
            .append("svg")
            .attr("width", this._get_div_width())
            .attr("height", 0);
        this.svg = d3v5.select(this.div).selectAll("svg");
    }
    _get_div_width() {
        return d3v5.select(this.div).node().clientWidth;
    }
    add_contig(contig_id) {
        let _this = this;
        this.updateloop = this.updateloop.then(res => {
            return _this._add_contig(contig_id);
        }).catch(err => console.log(err));
        this.update();
        return this;
    }
    async _add_contig(contig_id) {
        let _this = this;
        let res = await fetch(EXPERIMENT_DATA_PATH + "/" + contig_id + ".json");
        // Many contig does not have experiment data, so that the request will fail
        if (!res.ok) {
            console.log("experiment data not found for", contig_id);
            return;
        }
        res = await res.json();
        console.log(res);
        let exp_datas = res.experimentsbygene;
        exp_datas.forEach(dat => {
            _this.data.push(dat);
        });
        console.log(this.data);
        return;
    }
    remove_contig(contig_id) {
        let _this = this;
        this.updateloop = this.updateloop.then(res => {
            return _this._remove_contig(contig_id);
        }).catch(err => console.log(err));
        this.update();
        return this;
    }
    async _remove_contig(contig_id) {
        console.log("remove contig", contig_id);
        let _this = this;
        this.data = this.data.filter(d => d.contig_id != contig_id);
        console.log(this.data);
        return;
    }
    update() {
        let _this = this;
        this.updateloop = this.updateloop.then(res => {
            return _this._update();
        }).catch(err => console.log(err));
        return this;
    }
    async _update() {
        const svg = this.svg;
        const svg_width  = this._get_div_width();
        const experiment_info = await this.experiment_info;
        const experiment_fields = experiment_info.field_names;
        //console.log(experiment_info);
        const data = this.data;

        const margin = {left: 140, right: 20, top: 20, bottom: 10};
        const width = svg_width - margin.left - margin.right;
        const cell_height = 20;
        const height = cell_height * data.length;
        const svg_height = height + margin.top + margin.bottom;

        // // No data available
        // if (this.data.length === 0) {
        //     svg.remove();
        //     return;
        // }

        // Resize the svg container
        svg
            .attr("width", svg_width)
            .attr("height", svg_height);

        // Let's draw something
        let scaleYByGeneId = d3v5.scaleBand()
            .range([margin.top, margin.top + height])
            .domain(data.map(d => d.gene_id))
            .padding(0.01);
        let scaleXByField = d3v5.scaleBand()
            .range([margin.left, margin.left + width])
            .domain(experiment_fields)
            .padding(0.01);

        let rect_data = [];
        data.forEach(row => {
            let gene_id = row.gene_id;
            row.experiment_data.forEach((d, i) => {
                let field_name = experiment_fields[i];
                let value = d;
                rect_data.push({
                    "value":   value,
                    "field":   field_name,
                    "gene_id": gene_id
                });
            });
        });

        // TODO: thinking more about this color scale, there are some extreme values here,
        // so I am currently using quantile color scale.
        // Also there are many null values
        let _scaleColorByValue = d3v5.scaleQuantile()
            .domain(rect_data.map(d => d.value).filter(d => d !== null))
            .range(['blue', 'beige', 'red']);
        let scaleColorByValue = function(v) {
            if (v === null) return "grey";
            return _scaleColorByValue(v);
        };
        let rectSel = svg.selectAll("rect.heatmap_rect")
            .data(rect_data);
        rectSel
            .enter()
            .append("rect")
            .classed("heatmap_rect", true)
            .merge(rectSel)
            .attr("x", d => scaleXByField(d.field))
            .attr("y", d => scaleYByGeneId(d.gene_id))
            .attr("width", scaleXByField.bandwidth())
            .attr("height", scaleYByGeneId.bandwidth())
            .attr("fill", d => scaleColorByValue(d.value));

        rectSel.exit().remove();

        let topaxisSel = svg.selectAll("g.topaxis")
            .data(["onlyone"]);
        topaxisSel.enter()
            .append("g")
            .classed("topaxis", true)
            .merge(topaxisSel)
            .attr("transform", `translate(0,${margin.top})`)
            .call(d3v5.axisTop(scaleXByField));

        let leftaxisSel = svg.selectAll("g.leftaxis")
            .data(["onlyoneyes"]);
        leftaxisSel.enter()
            .append("g")
            .classed("leftaxis", true)
            .merge(leftaxisSel)
            .attr("transform", `translate(${margin.left}, 0)`)
            .call(d3v5.axisLeft(scaleYByGeneId));

        // //let heatmap_frame = svg.selectAll("");

        // let hhhrect = svg.selectAll("rect.hhh")
        //     .data(['wwwwwww']);
        // hhhrect
        //     .enter()
        //     .append("rect")
        //     .attr("class", "hhh")
        //     .merge(hhhrect)
        //     .attr("width", width)
        //     .attr("height", height)
        //     .attr("x", margin.left)
        //     .attr("y", margin.top)
        //     .style("fill", "green");

        // return;
    }
}

// This class currently handles the table panel and top-level logic for the application
class IndexTable {
    constructor(div) {
        this.dom = div;
        // Options for the table
        this.gridOptions = null;
        this.option_filterSelected = false;

        // TODO: move outside of this class.
        this.vispanel = new TrackViewPanel("#vis");
        this.experimentpanel = new ExperimentPanel("#experiment-vis");
    }
    // This function will fetch the "index.json" and load the table
    async load() {
        let _this = this;
        let res;
        res = await fetch(INDEX_PATH);
        if (res.ok)
            res = await res.json();
        else
            throw "Network error";

        let colorScale_contiglength = d3v5.scaleSequential(d3v5.interpolateReds)
            .domain(res.map(d => d.contig_length));
        let colorScale_abundance = d3v5.scaleSequential(d3v5.interpolateReds)
            .domain(res.map(d => d.contig_abundance));

        let columnDefs = [
            {
                headerName: "Contig ID",
                field: "contig_id",
                sortable: true,
                checkboxSelection: true,
                //headerCheckboxSelection: true,
            },
            {
                headerName: "Length",
                headerTooltip: "Contig Length",
                field: "contig_length",
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
                headerName: "Abundance",
                headerTooltip: "Abundance",
                field: "contig_abundance",
                filter: "agNumberColumnFilter",
                type: "numericColumn",
                sortable: true,
                width: 100,
                cellStyle: function(params) {
                    let color = colorScale_abundance(params.value);
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

        // Call add_view or remove_view
        function onRowSelected(event) {
            let contig_id = event.data.contig_id;
            let is_checked = event.node.isSelected();

            if (is_checked) {
                _this.vispanel.add_view(contig_id);
                _this.experimentpanel.add_contig(contig_id);
            }
            else {
                _this.vispanel.remove_view(contig_id);
                _this.experimentpanel.remove_contig(contig_id);
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
                let contig_id = event.data.contig_id;
                let vis_dom = document.getElementById("viscontainer-" + contig_id);
                let the_board = _this.vispanel.ActiveViews.get(contig_id).board;
                d3v5.select(vis_dom).classed("tntboard-highlight", true);
            }
        }
        function onCellMouseOut(event) {
            let contig_id = event.data.contig_id;
            if (_this.vispanel.ActiveViews.has(contig_id)) {
                let vis_dom = document.getElementById("viscontainer-" + contig_id);
                let the_board = _this.vispanel.ActiveViews.get(contig_id).board;
                d3v5.select(vis_dom).classed("tntboard-highlight", false);
            }
        }

    }

    // // Calls set_gene_colorby method from TrackViewPanel
    // set_gene_colorby(what) {
    //     let _this = this;
    //     let vispanel = this.vispanel;
    //     vispanel.set_gene_colorby(what);
    // }
}

let table;

async function main() {
    table = new IndexTable("#contig-table");
    await table.load();
}
main();


//// Callback used in buttons, checkboxes, etc.
function binding_filterSelected(node) {
    table.option_filterSelected = node.checked;
    table.gridOptions.api.onFilterChanged();
}

// function binding_gene_colorby(node) {
//     let value = node.value;
//     console.log("binding colorby", value);
//     table.set_gene_colorby(value);
// }

