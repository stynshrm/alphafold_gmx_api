import os
import sys
from pathlib import Path
from absl import flags
from absl import app
import gmxapi as gmx

FLAGS = flags.FLAGS

MODEL_PRESETS = {
    'monomer': (
        'model_1',
        'model_2',
        'model_3',
        'model_4',
        'model_5',
    ),
    'monomer_ptm': (
        'model_1_ptm',
        'model_2_ptm',
        'model_3_ptm',
        'model_4_ptm',
        'model_5_ptm',
    ),
    'multimer': (
        'model_1_multimer_v2',
        'model_2_multimer_v2',
        'model_3_multimer_v2',
        'model_4_multimer_v2',
        'model_5_multimer_v2',
    ),
}
MODEL_PRESETS['monomer_casp14'] = MODEL_PRESETS['monomer']

cores = os.cpu_count()
maxh = 0.01
allocation_size = cores

try:
    from mpi4py import MPI
    rank_number = MPI.COMM_WORLD.Get_rank()
    ensemble_size = MPI.COMM_WORLD.Get_size()
except ImportError:
    rank_number = 0
    ensemble_size = 1
    rank_tag = ''
    MPI = None
else:
    rank_tag = 'rank{}:'.format(rank_number)
try:
    local_cpu_set_size = len(os.sched_getaffinity(0))
except (NotImplementedError, AttributeError):
    threads_per_rank = allocation_size // ensemble_size
else:
    threads_per_rank = min(
        local_cpu_set_size, allocation_size // ensemble_size)

def create_top(inp_pdb):
    # print("Creating topology for ", inp_pdb)
    args = ['pdb2gmx', '-ff', 'amber99sb-ildn', '-water', 'tip3p']
    input_files = {'-f': inp_pdb}
    output_files = {'-p': 'topol.top', '-i': 'posre.itp', '-o': 'conf.gro'}

    make_top = gmx.commandline_operation('gmx', args, input_files, output_files)
    make_top.run()
    # EDITCONF
    args = ['editconf', '-bt', 'dodecahedron', '-d', 2]
    input_files = { '-f': make_top.output.file['-o']}
    output_files = {'-o': 'boxed.gro'}
    editconf = gmx.commandline_operation('gmx', args, input_files, output_files)
    
    top_out_files = [editconf.output.file['-o'].result(), make_top.output.file['-p'].result()]
    return top_out_files

def gen_tpr(mdp_file, gro_file, top_file):
    args = ['grompp', '-maxwarn', 2]
    grompp_input_files = {'-f': mdp_file,
                          '-c': gro_file,
                          '-p': top_file}
    grompp_output_files = {'-o': 'run.tpr'}
    grompp = gmx.commandline_operation('gmx', args, input_files=[
                                       grompp_input_files], output_files=[grompp_output_files])
    tpr_file = grompp.output.file['-o'].result()
    return tpr_file

def solvate(top_outfiles):
    # gmx solvate -cp setup-EM2.gro -cs waterbox.gro -o solv.gro -p topol.top
    args = ['solvate', '-cs', 'spc216']
    input_files = {'-cp': top_outfiles[0],
                    '-p': top_outfiles[1]
                    }
    output_files = {'-o': 'solv.gro'}
    add_wat = gmx.commandline_operation('gmx', args, input_files, output_files)
    # add_wat.run()
    tpr_file = gen_tpr('/Users/ssharma/Wrk/gmx_API/api_alpha/alphafold-2.2.0/test_flags/output_dir/minim.mdp', 
                       add_wat.output.file['-o'], top_outfiles[1])
    # add ions
    # gmx genion -s ions.tpr -p ../top/4ake.top -pname NA -nname CL -neutral -conc 0.15 -o ionized.pdb
    args = ['genion', '-pname', 'NA', '-nname', 'CL', '-neutral', '-conc', '0.15']
    input_files = {'-s': tpr_file, 
                    '-p': top_outfiles[1]}
    output_files = {'-o': 'solv_ions.gro'}
    genion = gmx.commandline_operation('gmx', args, input_files=[input_files], output_files=[output_files], stdin='SOL\n')

    # assert genion.output.file['-o'].result()
    # print(genion.output.file['-o'].result())
    return [genion.output.file['-o'].result(), top_outfiles[1]]

def gmx_minimize(min_input_files):
    tpr_file = gen_tpr('/Users/ssharma/Wrk/gmx_API/api_alpha/alphafold-2.2.0/test_flags/output_dir/minim.mdp',
                       min_input_files[0], min_input_files[1])
    input_list = gmx.read_tpr(tpr_file)
    output_list = {'-c': 'confout.gro'}
    # md = gmx.mdrun(input_list, runtime_args={'-maxh': str(maxh)})
    md = gmx.mdrun(input_list, output_list)
    # return md.output.checkpoint
    md.run()
    md_run_path = os.path.dirname(md.output.trajectory.result())
    return [os.path.join(md_run_path, output_list['-c']), min_input_files[1]]

def gmx_md(md_input_files):
    N = 3
    maxh=0.01
    grompp_input_files = {'-f': '/Users/ssharma/Wrk/gmx_API/api_alpha/alphafold-2.2.0/test_flags/output_dir/grompp.mdp',
                          '-c': md_input_files[0],
                          '-p': md_input_files[1]}
    grompp = gmx.commandline_operation(
        'gmx',
        ['grompp'],
        input_files=[grompp_input_files] * N,
        output_files={'-o': 'run.tpr'})
    tpr_input = grompp.output.file['-o'].result()
    input_list = gmx.read_tpr(tpr_input)
    md = gmx.mdrun(input_list, runtime_args={'-maxh': str(maxh)})
    md.run()

def gmx_pipeline(inp_model):
    topol_prot = create_top(inp_model)
    topol_solv = solvate(topol_prot)
    minimized_prot = gmx_minimize(topol_solv)
    if FLAGS.run_moldyn:
        gmx_md(minimized_prot)
        print("="*10) 

