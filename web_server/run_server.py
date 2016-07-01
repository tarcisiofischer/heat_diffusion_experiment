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

    try:
        k = float(request.args.get('k', 385.0))
    except:
        k = 385.0
    try:
        rho = float(request.args.get('rho', 8000.0))
    except:
        rho = 8000.0
    try:
        cp = float(request.args.get('cp', 400.0))
    except:
        cp = 400.0

    solve(
        n_elements=25,
        total_simulation_time=500.0,
        after_timestep_callback=save_solution,
        k=k,
        rho=rho,
        cp=cp,
    )
    return json.dumps(solutions)


if __name__ == "__main__":
    app.run()
