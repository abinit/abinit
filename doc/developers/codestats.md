---
plotly: true
---

<div class="plotly-graph-div" id="abinit_stats_plot" style="width:90%;height:750px;"></div>

<script>
$(function() {
    Plotly.d3.json("../statistics.json", function(stats) {
        var mode = "lines+markers";
        var x = stats.dates;
        var trace1 = {x: x, y: stats.num_f90lines, mode: mode, name: "Number of F90 lines"};
        var trace2 = {x: x, y: stats.num_f90files, mode: mode, name: "Number of F90 files", yaxis: 'y2'};
        var trace3 = {x: x, y: stats.num_tests, mode: mode, name: "Number of tests", yaxis: 'y3'};
        var trace4 = {x: x, y: stats.targz_sizes, mode: mode, name: "Tarball size [Mb]",  yaxis: 'y4'};
        var data = [trace1, trace2, trace3, trace4];

        var layout = {
          //title: "Date released: date mentioned in the release notes",
          //xaxis: {tickvals: stats.date, ticktex: stats.version},
          legend: {traceorder: 'reversed'},
          yaxis: {domain: [0, 0.25]},
          yaxis2: {domain: [0.25, 0.5]},
          yaxis3: {domain: [0.5, 0.75]},
          yaxis4: {domain: [0.75, 1.0]}
        };

        Plotly.newPlot(document.getElementById('abinit_stats_plot'), data, layout, {showLink: false});
    });
});
</script>
