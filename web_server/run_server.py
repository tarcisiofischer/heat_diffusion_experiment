from flask import Flask

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
    solve(25, total_simulation_time=500.0, after_timestep_callback=save_solution)
    return json.dumps(solutions)


if __name__ == "__main__":
    app.run()
