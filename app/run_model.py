import subprocess


def run_model_1(model, encoding, epi, seq_input):
    subprocess.run(['python', 'epitope_b_cells_predictor/main.py',
                    '--mode', 'predict',
                    '--init', encoding,
                    '--model', model,
                    '--epi', epi,
                    '--test_seq_input', seq_input])


def run_model_2(model, encoding, epi, seq_input):
    subprocess.run(['python', 'epitope_b_cells_predictor/main.py',
                    '--mode', 'predict',
                    '--init', encoding,
                    '--model', model,
                    '--epi', epi,
                    '--test_pdb_list', seq_input])


def run_model_3(model, encoding, epi, seq_input):
    subprocess.run(['python', 'epitope_b_cells_predictor/main.py',
                    '--mode', 'predict',
                    '--init', encoding,
                    '--model', model,
                    '--epi', epi,
                    '--test_pdb_path', seq_input])
