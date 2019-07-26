
## Multi-Omics Vis

### Data processing

The `sample_data` directory contains the raw data files, and
`scripts/prepare-sampledata.R` is an ad-hoc script to transform
the raw data files and create json files for visualization at
`srcweb/sample_data/`.

In these json files, `index.json` is the file that contains data to
render the table. It currently looks like following, corresponding to
the three columns in the table:

```
[
  {
    "contig": "MAG06_1",
    "length": 4095,
    "number_gene": 6
  },
  ...
]
```

Other json files contain the data to render the track panels for each contig,
they look like the following:

```
{
  "meta": {
    "seqname": "MAG06_10",
    "seqlen": 15686
  },
  "gene_track": [
    {
      "gene_start": 33,
      ...
    },
    ...
  ],
  "diamond_track": [
    {
      "gene_ID": "MAG06_10_1",
      ...
    },
    ...
  ],
}
```

In future, these files should be replaced by a restful API backed by
a web server.

### TnT Libraries

TnT is a set of libraries for track and tree based visualization. These links
might be of interest:

- Documentation of TnT Board http://tntvis.github.io/tnt.board/
- Documentation of TnT Genome http://tntvis.github.io/tnt.genome/
- Documentation of TnT Tooltip http://tntvis.github.io/tnt.tooltip/

This application is using the "gene" feature implemented in TnT Genome to show
the tracks. The exact data format needed for gene track is not documented in TnT,
but there is an [example](https://github.com/tntvis/tnt.genome/blob/master/examples/custom-coordinates/index.html)
showing how to create a "gene" track or "transcripts" track with custom data.

### Ag-Grid

Currently I am using Ag-Grid to display the table. In future, it should support data
from a remote source and allow searching, filtering, etc.

However, one consideration is that it seems that remote data source is an enterprise
feature in ag-grid. In that case, we may switch to use DataTable library to display
the table and implement corresponding searching api on the server.

### srcweb

`srcweb` directory contains source of the web application. Currently it is merely static
files. However, it might be helpful to add an extra build step so that we can support
scss, make the code more modular, etc.

### omicsvis.js

`srcweb/js/omicsvis.js` contains the main logic for the application.
See the inline comments.

