from flask import Flask
from flask.globals import request

app = Flask(__name__)

@app.route("/")
def index():
    return app.send_static_file('index.html')


@app.route('/js/<path:path>')
def send_js(path):
    return app.send_static_file('/'.join(['js', path]))


@app.route("/run_simulation")
def run_simulation():
    solutions = []
    def save_solution(u, t):
        solutions.append(u.tolist())
        print(t)
    import json
    from solver.heat_diffusion_2d import solve

    def get_value(key, default):
        try:
            return float(request.args.get(key, default))
        except:
            return default

    solve(
        n_elements=25,
        total_simulation_time=500.0,
        after_timestep_callback=save_solution,

        # Physical properties
        k=get_value('k', 385.0),
        rho=get_value('rho', 8000.0),
        cp=get_value('cp', 400.0),

        # Boundary conditions
        top_temperature=get_value('top_temperature', 25.0),
        bottom_temperature=get_value('bottom_temperature', 25.0),
        left_temperature=get_value('left_temperature', 0.0),
        right_temperature=get_value('right_temperature', 0.0),
    )
    return json.dumps(solutions)


if __name__ == "__main__":
    app.run()
