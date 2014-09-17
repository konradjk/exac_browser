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

quality_chart_margin = {top: 10, right: 30, bottom: 30, left: 50},
    quality_chart_width = 500 - quality_chart_margin.left - quality_chart_margin.right,
    quality_chart_height = 250 - quality_chart_margin.top - quality_chart_margin.bottom;


function draw_histogram_d3(chart_data) {
    var x = d3.scale.linear()
        .domain([0, d3.max(chart_data)])
        .range([0, quality_chart_width]);

    // Generate a histogram using twenty uniformly-spaced bins.
    var data = d3.layout.histogram()
        .bins(x.ticks(20))
        (chart_data);
    //console.log('Initial data: ', data);
    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d.y; })])
        .range([quality_chart_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#quality_display_container').append("svg")
        .attr("width", quality_chart_width + quality_chart_margin.left + quality_chart_margin.right)
        .attr("height", quality_chart_height + quality_chart_margin.top + quality_chart_margin.bottom)
      .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + quality_chart_margin.left + "," + quality_chart_margin.top + ")");

    var bar = svg.selectAll(".bar")
        .data(data)
      .enter().append("g")
        .attr("class", "bar")
        .attr("transform", function(d) { return "translate(" + x(d.x) + "," + y(d.y) + ")"; });

    bar.append("rect")
        .attr("x", 1)
        .attr("width", x(data[0].dx) - 1)
        .attr("height", function(d) { return quality_chart_height - y(d.y); });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + quality_chart_height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);
}


function change_histogram(raw_chart_data, plot_object, variant_only) {
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
        .range([0, quality_chart_width]);

    // Generate a histogram using twenty uniformly-spaced bins.
    var data = d3.layout.histogram()
        .bins(x.ticks(20))
        (chart_data);

//    console.log('Old data: ', d3.select('#quality_display_container').select('svg').select('#inner_graph').selectAll('rect').data());
//    console.log('New data:', data);
    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d.y; })])
        .range([quality_chart_height, 0]);

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
        .attr("transform", "translate(0," + quality_chart_height + ")")
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
            return quality_chart_height - y(d.y);
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
                .attr("height", quality_chart_height - y(data[i].y));
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
gene_chart_margin = {top: 10, right: 30, bottom: 30, left: 80},
    gene_chart_margin_lower = {top: 5, right: gene_chart_margin.right, bottom: 5, left: gene_chart_margin.left},
    gene_chart_width = 1100 - gene_chart_margin.left - gene_chart_margin.right;

lower_gene_chart_height = 50 - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom,
    gene_chart_height = 300 - gene_chart_margin.top - gene_chart_margin.bottom - lower_gene_chart_height - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom;

function gene_chart(data, exon_data, variant_data) {
    var metric = 'mean';
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
    console.log('total exon length: ', total_exon_length_padded);
    console.log('width: ', gene_chart_width);

    var start_pos = exon_data[0].start;
    var exon_x_scale = d3.scale.linear()
        .domain([0, total_exon_length_padded - padding])
        .range([0, gene_chart_width]);

    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d[metric]; })])
        .range([gene_chart_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(exon_x_scale)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#gene_plot_container').append("svg")
        .attr("width", gene_chart_width + gene_chart_margin.left + gene_chart_margin.right)
        .attr("height", gene_chart_height + gene_chart_margin.top + gene_chart_margin.bottom)
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + gene_chart_margin.left + "," + gene_chart_margin.top + ")");

    var new_data = _.filter(data, function(x) { return x.exon_number >= 0; });

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
            if (d.exon_number > 0) {
                relative_start_pos += (d.exon_number)*padding;
            }
            return exon_x_scale(relative_start_pos);
        })
        .attr("width", 1)
        .attr("y", function(d) { return y(d[metric]); })
        .attr("height", function(d) { return gene_chart_height - y(d[metric]); });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + gene_chart_height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);

    var svg_outer = d3.select('#gene_plot_container').append("svg")
        .attr("width", gene_chart_width + gene_chart_margin_lower.left + gene_chart_margin_lower.right)
        .attr("height", lower_gene_chart_height)
        .append("g")
        .attr('id', 'track')
        .attr("transform", "translate(" + gene_chart_margin_lower.left + "," + 0 + ")");

    var exon_color = "lightsteelblue";
    svg_outer.append("line")
        .attr("y1", lower_gene_chart_height/2)
        .attr("y2", lower_gene_chart_height/2)
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
        .attr("height", lower_gene_chart_height);

//    console.log("Variant data", variant_data);

    var variant_size_scale = d3.scale.log()
        .domain([d3.min(variant_data, function(d) { return d.allele_freq; }), d3.max(variant_data, function(d) { return d.allele_freq; })])
        //Circle/Ellipse
        .range([lower_gene_chart_height/3, 2]);
        //Rectangle
//        .range([lower_gene_chart_height, 2]);

    svg_outer.selectAll("bar")
        .data(variant_data)
        .enter()
        .append("a")
        .attr('class', 'track_variant_link')
        .attr("xlink:href", function(d, i) { return "/variant/" + d.chrom + "-" + d.pos + "-" + d.ref + "-" + d.alt; })
        .attr("data-toggle", "tooltip")
        .attr('filter_status', function(d) {
            return d.filter;
        })
        .attr('category', function(d) {
            return d.category;
        })
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
                if (d.vep_annotations[0]['EXON'] != '') {
                    var variant_exon_info;
                    variant_exon_info = d.vep_annotations[0]['EXON'].split('/');
                    var variant_exon_number;
                    if (d.vep_annotations[0]['STRAND'] == '-1') {
                        variant_exon_number = variant_exon_info[1] - variant_exon_info[0];
                    } else {
                        variant_exon_number = variant_exon_info[0] - 1;
                    }
                    console.log(d);
                    console.log(tx_coord, " ", variant_exon_number);
                    return exon_x_scale(tx_coord + variant_exon_number*padding);
                } else {
                    //TODO: Implement correct intron drawing
                    var variant_intron_info;
                    variant_intron_info = d.vep_annotations[0]['INTRON'].split('/');

                    var variant_intron_number;
                    if (d.vep_annotations[0]['STRAND'] == '-1') {
                        variant_intron_number = variant_intron_info[1] - variant_intron_info[0];
                    } else {
                        variant_intron_number = variant_intron_info[0] - 1;
                    }
                    return 0;
                }
            }
        })
        .attr("cy", lower_gene_chart_height/2)
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
//        .attr("y", function(d, i) { return lower_gene_chart_height/2 - variant_size_scale(d.allele_freq)/2; } )
//        .attr("width", 2)
//        .attr("height", function(d, i) { return variant_size_scale(d.allele_freq); })
//        .attr("rx", 6)
//        .attr("ry", 6);
}

function change_gene_chart_variant_size(variant_data, change_to) {
    var svg_outer = d3.select('#gene_plot_container').select('#track');

    var variant_size_scale;
    if (change_to == 'inverse_af') {
        variant_size_scale = d3.scale.log()
            .domain([d3.min(variant_data, function (d) {
                return d.allele_freq;
            }), d3.max(variant_data, function (d) {
                return d.allele_freq;
            })])
            .range([lower_gene_chart_height / 3, 2]);
    } else {
        variant_size_scale = d3.scale.log()
            .domain([d3.min(variant_data, function (d) {
                return d.allele_freq;
            }), d3.max(variant_data, function (d) {
                return d.allele_freq;
            })])
            .range([2, lower_gene_chart_height / 3]);
    }
    svg_outer.selectAll("a")
        .selectAll("ellipse")
        .transition()
        .duration(500)
        .attr("ry", function(d, i) { return variant_size_scale(d.allele_freq); });
}

function change_gene_chart_metric(data, metric) {
    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d[metric]; })])
        .range([gene_chart_height, 0]);

    var svg = d3.select('#gene_plot_container').select('#inner_graph');

    var new_data = _.filter(data, function(x) { return x.exon_number >= 0; });

//    console.log('Exon data: ', exon_data);
    svg.selectAll("rect")
        .data(new_data)
        .transition()
        .duration(500)
        .attr("y", function(d) { return y(d[metric]); })
        .attr("height", function(d) { return gene_chart_height - y(d[metric]); });

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    svg.select(".y.axis")
        .transition()
        .duration(200)
        .call(yAxis);

}

function update_variants() {
    $('[category]').hide();
    var v = $('.consequence_display_buttons.active').attr('id').replace('consequence_', '').replace('_button', '');
    var f = $('#filtered_checkbox').is(":checked") ? '[filter_status]' : '[filter_status="PASS"]';
    $('[category=' + v + ']' + f).show();
}
