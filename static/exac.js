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
    var latLng = new google.maps.LatLng( 40.708762, -74.006731 );
    var center = new google.maps.LatLng(25, 140);

    var map = new google.maps.Map( $('#map_canvas')[0], {
        zoom: 1,
        center: center,
        mapTypeId: google.maps.MapTypeId.SATELLITE
    });

    var data = google.visualization.arrayToDataTable([
        [ 'Allele', 'Count' ],
        [ 'Ref', 100 ],
        [ 'Alt', 200 ]
    ]);

    var options = {
        fontSize: 8,
        backgroundColor: 'transparent',
        legend: 'none'
    };

    var marker = new ChartMarker({
        map: map,
        position: latLng,
        width: '60px',
        height: '60px',
        chartData: data,
        chartOptions: options
    });
};