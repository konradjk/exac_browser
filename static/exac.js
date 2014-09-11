if (typeof google !== 'undefined') {
    google.load('visualization', '1.0', {'packages': ['corechart']});

    function ChartMarker(options) {
        this.setValues(options);

        this.$inner = $('<div>').css({
            position: 'relative',
            left: '-50%', top: '-50%',
            width: options.width,
            height: options.height,
            fontSize: '1px',
            lineHeight: '1px',
            backgroundColor: 'transparent',
            cursor: 'default'
        });

        this.$div = $('<div>')
            .append(this.$inner)
            .css({
                position: 'absolute',
                display: 'none'
            });
    };
    ChartMarker.prototype = new google.maps.OverlayView;
    ChartMarker.prototype.onAdd = function () {
        $(this.getPanes().overlayMouseTarget).append(this.$div);
    };
    ChartMarker.prototype.onRemove = function () {
        this.$div.remove();
    };

    ChartMarker.prototype.draw = function () {
        var marker = this;
        var projection = this.getProjection();
        var position = projection.fromLatLngToDivPixel(this.get('position'));

        this.$div.css({
            left: position.x,
            top: position.y - 30,
            display: 'block'
        })

        this.chart = new google.visualization.PieChart(this.$inner[0]);
        this.chart.draw(this.get('chartData'), this.get('chartOptions'));
    };
}


function draw_histogram_d3(chart_data) {
    var margin = {top: 10, right: 30, bottom: 30, left: 50},
        width = 500 - margin.left - margin.right,
        height = 250 - margin.top - margin.bottom;

    var x = d3.scale.linear()
        .domain([0, d3.max(chart_data)])
        .range([0, width]);

    // Generate a histogram using twenty uniformly-spaced bins.
    var data = d3.layout.histogram()
        .bins(x.ticks(20))
        (chart_data);
    //console.log('Initial data: ', data);
    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d.y; })])
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#quality_display_container').append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var bar = svg.selectAll(".bar")
        .data(data)
      .enter().append("g")
        .attr("class", "bar")
        .attr("transform", function(d) { return "translate(" + x(d.x) + "," + y(d.y) + ")"; });

    bar.append("rect")
        .attr("x", 1)
        .attr("width", x(data[0].dx) - 1)
        .attr("height", function(d) { return height - y(d.y); });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);
}


function change_histogram(raw_chart_data, plot_object, variant_only) {
    var margin = {top: 10, right: 30, bottom: 30, left: 50},
        width = 500 - margin.left - margin.right,
        height = 250 - margin.top - margin.bottom;

    var chart_data = [];
    if (variant_only) {
        $.each(raw_chart_data[plot_object], function(i, x) {
            if (!_.isEqual(raw_chart_data['genotypes'][i], ['0', '0']) && !_.isEqual(raw_chart_data['genotypes'][i], ['.', '.'])) {
                chart_data.push(x);
            }
        });
    } else {
        chart_data = raw_chart_data[plot_object];
    }
    //console.log(chart_data.length);
    var x = d3.scale.linear()
        .domain([0, d3.max(chart_data)])
        .range([0, width]);

    // Generate a histogram using twenty uniformly-spaced bins.
    var data = d3.layout.histogram()
        .bins(x.ticks(20))
        (chart_data);

//    console.log('Old data: ', d3.select('#quality_display_container').select('svg').select('#inner_graph').selectAll('rect').data());
//    console.log('New data:', data);
    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d.y; })])
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#quality_display_container').select('svg').select('#inner_graph');
    old_data_len = svg.selectAll('rect').data().length;

    svg.select(".x.axis")
        .transition()
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    svg.select(".y.axis")
        .transition()
        .call(yAxis);

    var bar = svg.selectAll(".bar")
        .data(data)
        .transition()
        .duration(500)
        .attr("transform", function(d) { return "translate(" + x(d.x) + "," + y(d.y) + ")"; });

    var rects = svg.selectAll('rect')
        .data(data)
        .transition()
        .duration(500)
        .attr('x', 1)
        .attr("width", x(data[0].dx) - 1)
        .attr("height", function(d) {
            return height - y(d.y);
        });
    if (old_data_len < data.length) {
        for (var i = old_data_len; i < data.length; i++) {
//            console.log('Adding: ', i);
            svg.selectAll('.bar').select('g')
                .data([data[i]]).enter()
                .append('g')
                .attr('class', 'bar')
                .append('rect')
                .attr('fill', 'steelblue')
                .transition()
                .duration(500)
                .attr('x', x(data[i].x) + 1)
                .attr('y', y(data[i].y))
                .attr('width', x(data[0].dx) - 1)
                .attr("height", height - y(data[i].y));
        }
    } else if (old_data_len > data.length) {
        for (var i = old_data_len-1; i >= data.length; i--) {
//            console.log('Removing: ', i);
            svg.selectAll('.bar')
                .data(data)
                .exit()
                .transition()
                .duration(500)
                .remove();
        }
    }
    //console.log('Final data: ', svg.selectAll('rect').data());
}


if (typeof google !== 'undefined') {
    function map_initialize() {
        var center = new google.maps.LatLng(30, 140);

        var map = new google.maps.Map($('#frequency_map_container')[0], {
            zoom: 1,
            center: center,
            mapTypeId: google.maps.MapTypeId.ROADMAP,
            panControl: false,
            zoomControl: false,
            scrollwheel: false,
            disableDoubleClickZoom: true,
            mapTypeControl: false,
            draggable: false,
            streetViewControl: false
        });

        var mapStyle = [
            {
                featureType: "administrative",
                elementType: "geometry.fill",
                stylers: [
                    { visibility: "off" }
                ]
            },
            {
                featureType: "administrative",
                elementType: "geometry.stroke",
                stylers: [
                    { visibility: "off" }
                ]
            },
            {
                featureType: "administrative",
                elementType: "labels",
                stylers: [
                    { visibility: "off" }
                ]
            }
        ];
        var styledMap = new google.maps.StyledMapType(mapStyle);
        map.mapTypes.set('myCustomMap', styledMap);
        map.setMapTypeId('myCustomMap');

        redraw_map(map);
        return map;
    }

    function redraw_map(map) {
        var options = {
            fontSize: 8,
            backgroundColor: 'transparent',
            legend: 'none'
        };


//    $.each(all_data, function(data) {
        var latLng = new google.maps.LatLng(40.708762, -74.006731);
        var data = google.visualization.arrayToDataTable([
            [ 'Allele', 'Count' ],
            [ 'Ref', Math.random() * 100 ],
            [ 'Alt', Math.random() * 100 ]
        ]);

        var marker = new ChartMarker({
            map: map,
            position: latLng,
            width: '60px',
            height: '60px',
            chartData: data,
            chartOptions: options
        });
//    });

    }
}

function gene_chart(data, exon_data, variant_data) {
    var margin = {top: 10, right: 30, bottom: 30, left: 50},
        margin_lower = {top: 5, right: margin.right, bottom: 5, left: margin.left},
        width = 1100 - margin.left - margin.right;

    var lower_graph_height = 50 - margin_lower.top - margin_lower.bottom,
        graph_height = 300 - margin.top - margin.bottom - lower_graph_height - margin_lower.top - margin_lower.bottom;

    var transcript = exon_data[0].transcript_id;
    var padding = 20;
    var total_exon_length = 0;
    var total_exon_length_padded = 0;
    var running_exon_length = [];
    $.each(exon_data, function(i, x) {
        total_exon_length += (x.stop - x.start + 1);
        total_exon_length_padded += (x.stop - x.start + 1 + padding);
        running_exon_length.push(total_exon_length);
    });
//    console.log("Total length: ", total_exon_length);
//    console.log("Total length padded: ", total_exon_length_padded);
//    console.log("Running lengths: ", running_exon_length);

    var start_pos = exon_data[0].start;
    var exon_x_scale = d3.scale.linear()
        .domain([0, total_exon_length_padded - padding])
        .range([0, width]);

    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d[0]; })])
        .range([graph_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(exon_x_scale)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#gene_plot_container').append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", graph_height + margin.top + margin.bottom)
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var new_data = _.filter(data, function(x) { return x[1] >= 0; });

//    console.log('Exon data: ', exon_data);
    svg.selectAll("bar")
        .data(new_data)
        .enter()
        .append("rect")
        .attr('class', 'main_plot_bars')
        .style("fill", "steelblue")
        .attr("x", function(d, i) {
//            console.log("i: ", i);
//            console.log("d: ", d);
            var relative_start_pos = i;
            if (d[1] > 0) {
                relative_start_pos += (d[1])*padding;
            }
            return exon_x_scale(relative_start_pos);
        })
        .attr("width", 1)
        .attr("y", function(d) { return y(d[0]); })
        .attr("height", function(d) { return graph_height - y(d[0]); });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + graph_height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);

    var svg_outer = d3.select('#gene_plot_container').append("svg")
        .attr("width", width + margin_lower.left + margin_lower.right)
        .attr("height", lower_graph_height)
        .append("g")
        .attr('id', 'track')
        .attr("transform", "translate(" + margin_lower.left + "," + 0 + ")");

    var exon_color = "lightsteelblue";
    svg_outer.append("line")
        .attr("y1", lower_graph_height/2)
        .attr("y2", lower_graph_height/2)
        .attr("x1", 0)
        .attr("x2", exon_x_scale(data.length))
        .attr("stroke-width", 10)
        .attr("stroke", exon_color);

    svg_outer.selectAll("bar")
        .data(exon_data)
        .enter()
        .append("rect")
        .attr('class', 'track_bar')
        .style("fill", exon_color)
        .attr("x", function(d, i) {
            var relative_start_pos;
            if (i == 0) {
                relative_start_pos = 0;
            } else {
                relative_start_pos = running_exon_length[i-1] + i*padding;
            }
            return exon_x_scale(relative_start_pos);
        })
        .attr("y", 0)
        .attr("rx", 6)
        .attr("ry", 6)
        .attr("width", function(d, i) { return exon_x_scale(d.stop-d.start+1); })
        .attr("height", lower_graph_height);

//    console.log("Variant data", variant_data);

    var variant_size_scale = d3.scale.log()
        .domain([d3.min(variant_data, function(d) { return d.allele_freq; }), d3.max(variant_data, function(d) { return d.allele_freq; })])
        //Circle/Ellipse
        .range([lower_graph_height/3, 2]);
        //Rectangle
//        .range([lower_graph_height, 2]);

    svg_outer.selectAll("bar")
        .data(variant_data)
        .enter()
        .append("a")
        .attr('class', 'track_variant_link')
        .attr("xlink:href", function(d, i) { return "/variant/" + d.chrom + "-" + d.pos + "-" + d.ref + "-" + d.alt; })
        .attr("data-toggle", "tooltip")
        .attr("title", function(d) {
            return d.vep_annotations[0]['Consequence'];
        })
        //Circle
//        .append("circle")
        //Ellipse
        .append("ellipse")
        .attr("class", "track_variant")
        .style("fill", "darkred")
        .style("opacity", 0.5)
        .attr("cx", function(d, i) {
            var tx_coord = d.transcript_coordinates[transcript];
            if (tx_coord == 0) {
                return -1000;
            } else {
                var variant_exon_number = d.vep_annotations[0]['EXON'].split('/')[0] - 1;
                return exon_x_scale(tx_coord + variant_exon_number*padding);
            }
        })
        .attr("cy", lower_graph_height/2)
        //Circle
//        .attr("r", function(d, i) { return variant_size_scale(d.allele_freq); })
        //Ellipse
        .attr("rx", 2)
        .attr("ry", function(d, i) { return variant_size_scale(d.allele_freq); });
        //Rectangle
//        .append("rect")
//        .attr("class", "track_variant")
//        .style("fill", "darkred")
//        .style("opacity", 0.5)
//        .attr("x", function(d, i) {
//            var tx_coord = d.transcript_coordinates[transcript];
//            if (tx_coord == 0) {
//                return -1000;
//            } else {
//                var variant_exon_number = d.vep_annotations[0]['EXON'].split('/')[0] - 1;
//                return exon_x_scale(tx_coord + variant_exon_number*padding);
//            }
//        })
//        .attr("y", function(d, i) { return lower_graph_height/2 - variant_size_scale(d.allele_freq)/2; } )
//        .attr("width", 2)
//        .attr("height", function(d, i) { return variant_size_scale(d.allele_freq); })
//        .attr("rx", 6)
//        .attr("ry", 6);
}

function change_variant_size(change_to) {
    var svg_outer = d3.select('#gene_plot_container').select('g');

    var variant_size_scale;
    if (change_to == 'proportional') {
        variant_size_scale = d3.scale.log()
            .domain([d3.min(variant_data, function (d) {
                return d.allele_freq;
            }), d3.max(variant_data, function (d) {
                return d.allele_freq;
            })])
            .range([lower_graph_height / 3, 2]);
    } else {
        variant_size_scale = d3.scale.log()
            .domain([d3.min(variant_data, function (d) {
                return d.allele_freq;
            }), d3.max(variant_data, function (d) {
                return d.allele_freq;
            })])
            .range([2, lower_graph_height / 3]);
    }
    svg_outer.selectAll("a")
        .selectAll("ellipse")
        .attr("ry", function(d, i) { return variant_size_scale(d.allele_freq); });
}

function gene_chart_separate(data, exon_data, variant_data) {
    var margin = {top: 10, right: 30, bottom: 30, left: 50},
        margin_lower = {top: 5, right: margin.right, bottom: 5, left: margin.left},
        width = 1100 - margin.left - margin.right;

    var lower_graph_height = 50 - margin_lower.top - margin_lower.bottom,
        graph_height = 300 - margin.top - margin.bottom - lower_graph_height - margin_lower.top - margin_lower.bottom;

    var x = d3.scale.linear()
        .domain([0, data.length])
        .range([0, width]);

    var y = d3.scale.linear()
        .domain([0, d3.max(data)])
        .range([graph_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#gene_plot_container').append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", graph_height + margin.top + margin.bottom)
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    svg.selectAll("bar")
        .data(data)
        .enter()
        .append("rect")
        .attr('class', 'main_plot_bars')
        .style("fill", "steelblue")
        .attr("x", function(d, i) { return x(i); })
        .attr("width", 1)
        .attr("y", function(d) { return y(d); })
        .attr("height", function(d) { return graph_height - y(d); });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + graph_height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);

    var svg_outer = d3.select('#gene_plot_container').append("svg")
        .attr("width", width + margin_lower.left + margin_lower.right)
        .attr("height", lower_graph_height)
        .append("g")
        .attr('id', 'track')
        .attr("transform", "translate(" + margin_lower.left + "," + 0 + ")");

    var exon_color = "lightsteelblue";
    svg_outer.append("line")
        .attr("y1", lower_graph_height/2)
        .attr("y2", lower_graph_height/2)
        .attr("x1", 0)
        .attr("x2", x(data.length))
        .attr("stroke-width", 10)
        .attr("stroke", exon_color);

    var transcript = exon_data[0].transcript_id;
    var padding = 20;
    var total_exon_length = 0;
    var total_exon_length_padded = 0;
    var running_exon_length = [];
    $.each(exon_data, function(i, x) {
        total_exon_length += (x.stop - x.start + 1);
        total_exon_length_padded += (x.stop - x.start + 1 + padding);
        running_exon_length.push(total_exon_length);
    });
    console.log("Total length: ", total_exon_length);
    console.log("Total length padded: ", total_exon_length_padded);
    console.log("Running lengths: ", running_exon_length);

    var start_pos = exon_data[0].start;
    var exon_x_scale = d3.scale.linear()
        .domain([0, total_exon_length_padded - padding])
        .range([0, width]);

//    var exon_data = [0, x(data.length/2)];
    console.log("Exon data", exon_data);
    svg_outer.selectAll("bar")
        .data(exon_data)
        .enter()
        .append("rect")
        .attr('class', 'track_bar')
        .style("fill", exon_color)
        .attr("x", function(d, i) {
            var relative_start_pos;
            if (i == 0) {
                relative_start_pos = 0;
            } else {
                relative_start_pos = running_exon_length[i-1] + i*padding;
            }
            console.log(relative_start_pos);
            return exon_x_scale(relative_start_pos);
        })
        .attr("y", 0)
        .attr("rx", 6)
        .attr("ry", 6)
        .attr("width", function(d, i) { return exon_x_scale(d.stop-d.start+1); })
        .attr("height", lower_graph_height);

    console.log("Variant data", variant_data);

    var variant_size_scale = d3.scale.log()
        .domain([d3.min(variant_data, function(d) { return d.allele_freq; }), d3.max(variant_data, function(d) { return d.allele_freq; })])
        .range([lower_graph_height/3, 2]);

    svg_outer.selectAll("bar")
        .data(variant_data)
        .enter()
        .append("a")
        .attr('class', 'track_variant_link')
        .attr("xlink:href", function(d, i) { return "/variant/" + d.chrom + "-" + d.pos + "-" + d.ref + "-" + d.alt; })
        .attr("data-toggle", "tooltip")
        .attr("title", function(d) {
            return d.vep_annotations[0]['Consequence'];
        })
        .append("circle")
        .attr("class", "track_variant")
        .style("fill", "darkred")
        .style("opacity", 0.5)
        .attr("r", function(d, i) { return variant_size_scale(d.allele_freq); })
        .attr("cx", function(d, i) {
            var tx_coord = d.transcript_coordinates[transcript];
            if (tx_coord == 0) {
                return -1000;
            } else {
                var variant_exon_number = d.vep_annotations[0]['EXON'].split('/')[0] - 1;
                return exon_x_scale(tx_coord + variant_exon_number*padding);
            }
        })
        .attr("cy", lower_graph_height/2);

}