from absl import flags
from absl import app
import shutil

FLAGS = flags.FLAGS

flags.DEFINE_string('fasta_paths', None, 'Paths to FASTA files')
flags.DEFINE_string('output_dir', None, 'Paths to output files')
flags.DEFINE_string('data_dir', None, 'Path to directory of supporting data.')
flags.DEFINE_string('exe_python', shutil.which('python'),
                    'Path to the hmmsearch executable.')
flags.DEFINE_string('uniref90_database_path', None, 'UNIREF')
flags.DEFINE_enum('model_preset', 'monomer',
                  ['monomer', 'monomer_casp14', 'monomer_ptm', 'multimer'],
                  'Choose preset model configuration - the monomer model, '
                  'the monomer model with extra ensembling, monomer model with '
                  'pTM head, or multimer model')
flags.DEFINE_integer('num_multimer_predictions_per_model', 5, 'How many '
                     'predictions (each with a different random seed) will be '
                     'generated per model. E.g. if this is 2 and there are 5 '
                     'models then there will be 10 predictions per input. '
                     'Note: this FLAG only applies if model_preset=multimer')
flags.DEFINE_boolean('run_relax', True, 'Whether to run the final relaxation '
                     'step on the predicted models. Turning relax off might '
                     'result in predictions with distracting stereochemical '
                     'violations but might help in case you are having issues '
                     'with the relaxation stage.')
def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')
    print('DATA_DIR  : ', FLAGS.data_dir)
    print('OUTPUT_DIR: ', FLAGS.output_dir)
    print('PYTHON_EXE: ', FLAGS.exe_python)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'fasta_paths',
      'data_dir',
  ])

  app.run(main)


