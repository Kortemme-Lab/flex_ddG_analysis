class ddg_benchmark_run:
    def __init__(self, prediction_set_name, step_multiplication_factor):
        self.prediction_set_name = prediction_set_name
        self.step_multiplication_factor = step_multiplication_factor

talaris_60k_run = ddg_benchmark_run(
    'zemu_1.2-60000_rscript_validated-t14', # Predicton run name
    5 # Step multiplication factor
)

all_runs = [ talaris_60k_run ]
