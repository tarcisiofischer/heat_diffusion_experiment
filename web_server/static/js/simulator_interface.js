function int_to_color(i)
{
    return "#" + ("0" + i.toString(16)).slice(-2) + "00" + ("0" + (255 - i).toString(16)).slice(-2);
}

function paint_simulation_result(values_2d)
{
    var canvas = document.getElementById('simulation_canvas');
    var context = canvas.getContext("2d");
    var pixel_size = 8;
    var min_value = 0.0;
    var max_value = 25.0;

    for (i = 0; i < values_2d.length; i++) {
        for (j = 0; j < values_2d.length; j++) {
            var normalized_value = (values_2d[i][j] - min_value) / (max_value - min_value);
            context.fillStyle = int_to_color(Math.floor(normalized_value * 255));
            context.fillRect(pixel_size * j, pixel_size * i, pixel_size, pixel_size);
        }
    }
}

var saved_m;
function run_again()
{
    if (saved_m)
        show_animation(saved_m);
}

// Poor man's lock
var show_animation_lock = false;
function show_animation(m)
{
    if (show_animation_lock)
        return;
    show_animation_lock = true;
    var l = 0;
    function myFunction() {
        paint_simulation_result(m[l]);
        if (l < m.length - 1) {
            setTimeout(myFunction, 10);
        } else {
            show_animation_lock = false;
            console.log("Done.");
        }
        l++;
    }
    myFunction();
}

function run_simulation()
{
    $.ajax({
        url: "http://127.0.0.1:5000/run_simulation",
        dataType: "json",
    }).done(function(m) {
        saved_m = m;
        show_animation(m);
    });
}
