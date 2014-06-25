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
function map_initialize() {
    var center = new google.maps.LatLng(30, 140);

    var map = new google.maps.Map( $('#map_canvas')[0], {
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
};

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

};