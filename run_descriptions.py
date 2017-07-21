class ddg_benchmark_run:
    def __init__(self, prediction_set_name, step_multiplication_factor):
        self.prediction_set_name = prediction_set_name
        self.step_multiplication_factor = step_multiplication_factor

talaris_60k_run = ddg_benchmark_run(
    'zemu_1.2-60000_rscript_validated-t14', # Predicton run name
    5 # Step multiplication factor
)

talaris_control_run = ddg_benchmark_run(
    'zemu_control', # Predicton run name
    1 # Step multiplication factor
)

ref_run = ddg_benchmark_run(
    'zemu_1.2-60000_rscript_validated-ref', # Predicton run name
    5 # Step multiplication factor
)

all_runs = [
    talaris_control_run,
    talaris_60k_run,
    # ref_run,
]
