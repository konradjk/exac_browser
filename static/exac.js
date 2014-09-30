quality_chart_margin = {top: 10, right: 30, bottom: 50, left: 50},
    quality_chart_width = 500 - quality_chart_margin.left - quality_chart_margin.right,
    quality_chart_height = 250 - quality_chart_margin.top - quality_chart_margin.bottom;


function draw_histogram_d3(data) {
    console.log(data);
    var x = d3.scale.linear()
        .domain([d3.min(data, function(d) { return d[0]; }), d3.max(data, function(d) { return d[0]; })])
        .range([0, quality_chart_width]);

    var bar_width = 2*x(data[1][0] - data[0][0]);
    var y = d3.scale.linear()
        .domain([d3.min(data, function(d) { return d[1]; }), d3.max(data, function(d) { return d[1]; })])
        .range([quality_chart_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#quality_display_container');
    if (svg.selectAll('rect').length == 0 || svg.selectAll('rect')[0].length == 0) {
        svg = d3.select('#quality_display_container').append("svg")
            .attr("width", quality_chart_width + quality_chart_margin.left + quality_chart_margin.right)
            .attr("height", quality_chart_height + quality_chart_margin.top + quality_chart_margin.bottom)
            .append("g")
            .attr('id', 'inner_graph')
            .attr("transform", "translate(" + quality_chart_margin.left + "," + quality_chart_margin.top + ")");

        var bar = svg.selectAll(".bar")
            .data(data)
            .enter().append("g")
            .attr("class", "bar");

        bar.append("rect")
            .attr("x", function(d) { return x(d[0]); })
            .attr("width", bar_width)
            .attr("height", function(d) { return quality_chart_height - y(d[1]); })
            .attr("y", function(d) { return y(d[1]); });

        svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + quality_chart_height + ")")
            .call(xAxis);

        svg.append("g")
            .attr("class", "y axis")
            .call(yAxis);
    } else {
        svg = d3.select('#quality_display_container').select('svg').select('#inner_graph');
        svg.select(".x.axis")
            .transition()
            .attr("transform", "translate(0," + quality_chart_height + ")")
            .call(xAxis);

        svg.select(".y.axis")
            .transition()
            .call(yAxis);

        svg.selectAll('rect')
            .data(data)
            .transition()
            .duration(500)
            .attr("x", function(d) { return x(d[0]); })
            .attr("width", bar_width)
            .attr("height", function(d) { return quality_chart_height - y(d[1]); })
            .attr("y", function(d) { return y(d[1]); });
    }
}

function draw_region_coverage(raw_data, metric, ref) {
    if (raw_data.length > 1) {
        var data = raw_data;
        var chart_width = _.min([quality_chart_width, data.length*30]);
        var x = d3.scale.linear()
            .domain([0, data.length])
            .range([0, chart_width]);

        var y = d3.scale.linear()
            .domain([0, d3.max(data, function(d) { return d[metric]; })])
            .range([quality_chart_height, 0]);

        var xAxis = d3.svg.axis()
            .scale(x)
            .orient("bottom");

        var yAxis = d3.svg.axis()
            .scale(y)
            .orient("left");

        var svg = d3.select('#region_coverage');

        if (svg.selectAll('rect').length == 0 || svg.selectAll('rect')[0].length == 0) {
            svg = d3.select('#region_coverage').append("svg")
            .attr("width", chart_width  + quality_chart_margin.left + quality_chart_margin.right)
            .attr("height", quality_chart_height + quality_chart_margin.top + quality_chart_margin.bottom)
            .append("g")
            .attr('id', 'inner_graph')
            .attr("transform", "translate(" + quality_chart_margin.left + "," + quality_chart_margin.top + ")");

            var bar = svg.selectAll(".bar")
                .data(data)
                .enter().append("g")
                .attr("class", "bar");

            bar.append("rect")
                .attr("x", function(d, i) { return x(i); })
                .attr("width", chart_width/data.length - 1)
                .attr("height", function(d) { return quality_chart_height - y(d[metric]); })
                .attr("y", function(d) { return y(d[metric]); });

            xAxis = d3.svg.axis()
                .scale(x)
                .tickFormat(function(d) { return ref[d]; })
                .innerTickSize(0)
                .orient("bottom");

            svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + quality_chart_height + ")")
                .call(xAxis);

            svg.append("g")
                .attr("class", "y axis")
                .call(yAxis);
        } else {
            svg = d3.select('#region_coverage').select('svg').select('#inner_graph');
            svg.select(".y.axis")
                .transition()
                .call(yAxis);

            svg.selectAll('rect')
                .data(data)
                .transition()
                .duration(500)
                .attr("x", function(d, i) { return x(i); })
                .attr("width", chart_width/data.length - 1)
                .attr("height", function(d) { return quality_chart_height - y(d[metric]); })
                .attr("y", function(d) { return y(d[metric]); });
        }
    } else {
        var data1 = {};
        $.each(raw_data[0], function(d, i) {
            var num = parseInt(d);
            if (!isNaN(num)) {
                data1[d] = raw_data[0][d];
            }
        });
        var data2 = {};
        data2['mean'] = raw_data[0]['mean'];
        data2['median'] = raw_data[0]['median'];

        var coverages = Object.keys(data1);
        var other_labels = Object.keys(data2);
        var all_labels = coverages.concat([''], other_labels);

        var chart_width = quality_chart_width;
        var total_data_length = coverages.length + other_labels.length + 1;
        var x = d3.scale.linear()
            .domain([0, total_data_length])
            .range([0, chart_width]);

        var y1 = d3.scale.linear()
            .domain([0, d3.max(coverages, function(d) { return data1[d]; })])
            .range([quality_chart_height, 0]);

        var y2 = d3.scale.linear()
            .domain([0, d3.max(other_labels, function(d) { return data2[d]; })])
            .range([quality_chart_height, 0]);

        var xAxis = d3.svg.axis()
            .scale(x)
            .tickFormat(function(d) { return all_labels[d]; })
            .orient("bottom");

        var yAxis1 = d3.svg.axis()
            .scale(y1)
            .orient("left");

        var yAxis2 = d3.svg.axis()
            .scale(y2)
            .orient("right");

        svg = d3.select('#region_coverage').append("svg")
        .attr("width", chart_width + quality_chart_margin.left + quality_chart_margin.right)
        .attr("height", quality_chart_height + quality_chart_margin.top + quality_chart_margin.bottom)
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + quality_chart_margin.left + "," + quality_chart_margin.top + ")");

        var bar = svg.selectAll(".bar")
            .data(coverages)
            .enter().append("g")
            .attr("class", "bar");

        bar.append("rect")
            .attr("x", function(d, i) { return x(i); })
            .attr("width", chart_width/total_data_length)
            .attr("height", function(d) { return quality_chart_height - y1(data1[d]); })
            .attr("y", function(d) { return y1(data1[d]); });

        svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + quality_chart_height + ")")
            .call(xAxis)
            .selectAll("text")
            .attr("transform", "translate(0, 10) rotate(45)");

        var bar = svg.selectAll(".bar").select('g')
            .data(other_labels)
            .enter().append("g")
            .attr("class", "bar");

        bar.append("rect")
            .attr("x", function(d, i) { return x(i + coverages.length + 1); })
            .attr("width", chart_width/total_data_length)
            .attr("height", function(d) { return quality_chart_height - y2(data2[d]); })
            .attr("y", function(d) { return y2(data2[d]); });

        svg.append("g")
            .attr("class", "y axis")
            .call(yAxis1);

        svg.append("g")
            .attr("class", "y axis")
            .attr("transform", "translate(" + chart_width + " ,0)")
            .call(yAxis2);

        d3.select('#region_coverage').append("text").text("Axis");
    }
}

function update_variants() {
    $('[category]').hide();
    var v = $('.consequence_display_buttons.active').attr('id').replace('consequence_', '').replace('_variant_button', '');
    var f = $('#filtered_checkbox').is(":checked") ? '[filter_status]' : '[filter_status="PASS"]';
    $('[category=lof_variant]' + f).show();
    if (v == 'missense') {
        $('[category=missense_variant]' + f).show();
    } else if (v == 'other') {
        $('[category=missense_variant]' + f).show();
        $('[category=other_variant]' + f).show();
    }
}


function get_color(variant) {
    if (variant.category == 'lof_variant') {
        return "darkred";
    } else if (variant.category == 'missense_variant') {
        return "yellow";
    } else {
        return "green";
    }
}