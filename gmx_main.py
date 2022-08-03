import sys, os
import run_alpha
import gmxapi as gmx
#from alphafold.model import config
from absl import flags
from absl import app
import gmx_routines

argv = None
FLAGS = flags.FLAGS
flags_parser = app.parse_flags_with_usage
original_argv = sys.argv if argv is None else argv
args_to_main = flags_parser(original_argv)

flags.DEFINE_boolean('run_moldyn', True, 'Run GROMACS simulation')

# flags.DEFINE_integer('cores', os.cpu_count(),
#                      'The total number of cores allocated for the job.')

def read_flags(inp_params):
    for (key, value) in inp_params.items():
        FLAGS[key].value = value

def get_model_names():
    run_multimer_system = 'multimer' in FLAGS.model_preset
    if run_multimer_system:
        num_predictions_per_model = FLAGS.num_multimer_predictions_per_model
    else:
        num_predictions_per_model = 1
    if FLAGS.run_relax:
        prefix = "relaxed"
    else:
        prefix = "unrelaxed"
    model_runners = []
    model_names = gmx_routines.MODEL_PRESETS[FLAGS.model_preset]
    for model_name in model_names:
        for i in range(num_predictions_per_model):
            model_runners.append(f'{prefix}_{model_name}_pred_{i}.pdb')
    # print(model_runners)
    return model_runners

def run_gmx(alphafold_outdir):
    model_runners = get_model_names()
    for model in model_runners:
        pdb_model = os.path.join(alphafold_outdir, model)
        gmx_routines.gmx_pipeline(pdb_model)
        # a = gmx_routines.gmx_pipeline(pdb_model)

@gmx.function_wrapper(output={'pdb_str': str})
def run_alphafold(required_args: dict, output):
    read_flags(required_args)
    run_alpha.main(" ")
    output.pdb_str = FLAGS.output_dir


if __name__ == '__main__':
    run_alphafold = run_alphafold({
        'fasta_paths': "/home/ssharma/test.fasta",
        'output_dir': "/Users/ssharma/Wrk/gmx_API/api_alpha/alphafold-2.2.0/test_flags-2/output",
        'data_dir': "/home/ssharma/data/dir",
        'uniref90_database_path': "/usr/bin/ref90",
        'run_relax': False,
        'run_moldyn': True
    })
    run_gmx(run_alphafold.output.pdb_str.result())




