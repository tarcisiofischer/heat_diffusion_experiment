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
    // Physical properties inputs
    var k = $("#k_input")[0].value;
    var rho = $("#rho_input")[0].value;
    var cp = $("#cp_input")[0].value;

    // Boundary conditions inputs
    var top_temperature = $("#u_top_input")[0].value;
    var bottom_temperature = $("#u_bottom_input")[0].value;
    var left_temperature = $("#u_left_input")[0].value;
    var right_temperature = $("#u_right_input")[0].value;

    $.ajax({
        url: "/run_simulation",
        dataType: "json",
        data: {
            k: k,
            rho: rho,
            cp: cp,
            top_temperature: top_temperature,
            bottom_temperature: bottom_temperature,
            left_temperature: left_temperature,
            right_temperature: right_temperature,
        },
    }).done(function(m) {
        saved_m = m;
        show_animation(m);
    });
}
