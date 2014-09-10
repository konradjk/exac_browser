google.load('visualization', '1.0', {'packages':['corechart']});

function ChartMarker( options ) {
    this.setValues( options );

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
        .append( this.$inner )
        .css({
            position: 'absolute',
            display: 'none'
        });
};
ChartMarker.prototype = new google.maps.OverlayView;
ChartMarker.prototype.onAdd = function() {
    $( this.getPanes().overlayMouseTarget ).append( this.$div );
};
ChartMarker.prototype.onRemove = function() {
    this.$div.remove();
};

ChartMarker.prototype.draw = function() {
    var marker = this;
    var projection = this.getProjection();
    var position = projection.fromLatLngToDivPixel( this.get('position') );

    this.$div.css({
        left: position.x,
        top: position.y - 30,
        display: 'block'
    })

    this.chart = new google.visualization.PieChart( this.$inner[0] );
    this.chart.draw( this.get('chartData'), this.get('chartOptions') );
};


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


function map_initialize() {
    var center = new google.maps.LatLng(30, 140);

    var map = new google.maps.Map( $('#frequency_map_container')[0], {
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
        }, {
            featureType: "administrative",
            elementType: "geometry.stroke",
            stylers: [
                { visibility: "off" }
            ]
        }, {
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
        var latLng = new google.maps.LatLng( 40.708762, -74.006731 );
        var data = google.visualization.arrayToDataTable([
            [ 'Allele', 'Count' ],
            [ 'Ref', Math.random()*100 ],
            [ 'Alt', Math.random()*100 ]
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

function gene_coverage_chart(data) {
    var margin = {top: 10, right: 30, bottom: 30, left: 50},
        width = 500 - margin.left - margin.right,
        height = 250 - margin.top - margin.bottom;

    var x = d3.scale.linear()
        .domain([0, data.length])
        .range([0, width]);

    console.log('Data: ', data);
    var y = d3.scale.linear()
        .domain([0, d3.max(data)])
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#coverage_container').append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    svg.selectAll("bar")
      .data(data)
      .enter().append("rect")
      .style("fill", "steelblue")
      .attr("x", function(d, i) { return x(i); })
      .attr("width", 1)
      .attr("y", function(d) { return y(d); })
      .attr("height", function(d) { return height - y(d); });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);

}