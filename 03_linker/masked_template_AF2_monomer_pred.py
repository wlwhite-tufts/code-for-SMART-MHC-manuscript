#modified slightly from scripts made by nrbennett
import os
import mock
import numpy as np
import sys
import datetime

from typing import Any, Mapping, Optional, Sequence, Tuple
import collections
from collections import OrderedDict

from timeit import default_timer as timer
import argparse

import io
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import PDBParser
from Bio.PDB.mmcifio import MMCIFIO

sys.path.insert( 0, '<path/to/your>/alphafold' ) #alphafold installation location

import jax
import jax.numpy as jnp

from alphafold.common import residue_constants
from alphafold.common import protein
from alphafold.common import confidence
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.data import mmcif_parsing
from alphafold.model import data
from alphafold.model import config
from alphafold.model import model
from alphafold.data.tools import hhsearch

sys.path.append( '<path/to/your>/silent_tools' ) #siletn tools installation location
import silent_tools

# Run with Nate's ampere environment
from pyrosetta import *
from rosetta import *
init( '-in:file:silent_struct_type binary' )

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-fasta", required=True, type=str,
                        help="fasta file to predict")
    parser.add_argument("-outname", required=True, type=str,
                        help="base name for score and silent outputs")
    parser.add_argument("-interface1", required=True, type=str,
                        help="range of residues to use for side 1 of the interface for calculating interface pae")
    parser.add_argument("-interface2", required=True, type=str,
                        help="range of residues to use for side 2 of the interface for calculating interface pae")
    parser.add_argument("-template_pdb",  type=str, required=True,
                        help="Path to the template pdb file")
    parser.add_argument("-recycle", dest="recycle", type=int, default=3,
                        help="Number of recycle iterations to perform")
    parser.add_argument("-exclude_template", type=str, default='',
                        help="Which residues of the template to mask out. Comma-separated idxs or ranges. ex/ 1-3,5,8-30")
    parser.add_argument("-batch_size", type=int, default=8,
                        help="How many structures to read into AF2 at once")

    args = parser.parse_args()
    return args

args = get_args()

model_name = "model_1_ptm"

model_config = config.model_config(model_name)
model_config.data.eval.num_ensemble = 1

model_config.data.common.num_recycle = args.recycle
model_config.model.num_recycle = args.recycle

model_config.model.embeddings_and_evoformer.initial_guess = True

model_config.model.global_config.mixed_precision = False

model_config.data.common.max_extra_msa = 5
model_config.data.eval.max_msa_clusters = 5

model_params = data.get_model_haiku_params(model_name=model_name, data_dir="<path/to/your/alphafold_model_parameters_dir/>")
model_runner = model.RunModel(model_config, model_params)

def get_seq_from_pdb( pdb_fn, slash_for_chainbreaks ):
  to1letter = {
    "ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C',
    "GLN":'Q', "GLU":'E', "GLY":'G', "HIS":'H', "ILE":'I',
    "LEU":'L', "LYS":'K', "MET":'M', "PHE":'F', "PRO":'P',
    "SER":'S', "THR":'T', "TRP":'W', "TYR":'Y', "VAL":'V' }

  seq = ''
  with open(pdb_fn) as fp:
    for line in fp:
      if line.startswith("TER"):
        if not slash_for_chainbreaks: continue
        seq += "/"
      if not line.startswith("ATOM"):
        continue
      if line[12:16].strip() != "CA":
        continue
      resName = line[17:20]
      #
      seq += to1letter[resName]
  return seq

def af2_get_atom_positions( pdbfilename ) -> Tuple[np.ndarray, np.ndarray]:
  """Gets atom positions and mask from a list of Biopython Residues."""

  with open(pdbfilename, 'r') as pdb_file:
    lines = pdb_file.readlines()

  # indices of residues observed in the structure
  idx_s = [int(l[22:26]) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
  num_res = len(idx_s)

  all_positions = np.zeros([num_res, residue_constants.atom_type_num, 3])
  all_positions_mask = np.zeros([num_res, residue_constants.atom_type_num],
                                dtype=np.int64)

  residues = collections.defaultdict(list)
  # 4 BB + up to 10 SC atoms
  xyz = np.full((len(idx_s), 14, 3), np.nan, dtype=np.float32)
  for l in lines:
    if l[:4] != "ATOM":
        continue
    resNo, atom, aa = int(l[22:26]), l[12:16], l[17:20]

    residues[ resNo ].append( ( atom.strip(), aa, [float(l[30:38]), float(l[38:46]), float(l[46:54])] ) )

  for resNo in residues:

    pos = np.zeros([residue_constants.atom_type_num, 3], dtype=np.float32)
    mask = np.zeros([residue_constants.atom_type_num], dtype=np.float32)

    for atom in residues[ resNo ]:
      atom_name = atom[0]
      x, y, z = atom[2]
      if atom_name in residue_constants.atom_order.keys():
        pos[residue_constants.atom_order[atom_name]] = [x, y, z]
        mask[residue_constants.atom_order[atom_name]] = 1.0
      elif atom_name.upper() == 'SE' and res.get_resname() == 'MSE':
        # Put the coordinates of the selenium atom in the sulphur column.
        pos[residue_constants.atom_order['SD']] = [x, y, z]
        mask[residue_constants.atom_order['SD']] = 1.0

    idx = idx_s.index(resNo)
    all_positions[idx] = pos
    all_positions_mask[idx] = mask
  # _check_residue_distances(
  #     all_positions, all_positions_mask, max_ca_ca_distance) # AF2 checks this but if we want to allow massive truncations we don't want to check this

  return all_positions, all_positions_mask

def af2_all_atom_from_struct( template_pdb, template_seq, template_mask, query_seq, just_target=False ):
  all_atom_positions, all_atom_mask = af2_get_atom_positions( template_pdb )

  all_atom_positions = np.split(all_atom_positions, all_atom_positions.shape[0])

  templates_all_atom_positions = []

  # Initially fill will all zero values
  for _ in query_seq:
    templates_all_atom_positions.append(
        jnp.zeros((residue_constants.atom_type_num, 3)))

  for idx in range( -1, -min(len(query_seq),len(template_seq)) - 1, -1 ): #matches the C-term of query with C-term of template

    if not template_mask[idx]: continue

    templates_all_atom_positions[ idx ] = all_atom_positions[ idx ][0] # assign target indices to template coordinates

  return jnp.array(templates_all_atom_positions)

def template_from_struct( template_pdb, template_seq, query_seq, mask ):

  ret_all_atom_positions, ret_all_atom_mask = af2_get_atom_positions( template_pdb )

  all_atom_positions = np.split(ret_all_atom_positions, ret_all_atom_positions.shape[0])
  all_atom_masks = np.split(ret_all_atom_mask, ret_all_atom_mask.shape[0])
  
  output_templates_sequence = []
  output_confidence_scores = []
  templates_all_atom_positions = []
  templates_all_atom_masks = []

  # Initially fill will all zero values
  for _ in query_seq:
    templates_all_atom_positions.append(
        np.zeros((residue_constants.atom_type_num, 3)))
    templates_all_atom_masks.append(np.zeros(residue_constants.atom_type_num))
    output_templates_sequence.append('-')
    output_confidence_scores.append(-1)
  
  confidence_scores = []
  for _ in query_seq: confidence_scores.append( 9 )

  for idx in range( -1, -min(len(query_seq),len(template_seq)) - 1, -1 ): #matches the C-term of query with C-term of template

    if not mask[ idx ]: continue

    templates_all_atom_positions[ idx ] = all_atom_positions[ idx ][0] # assign target indices to template coordinates
    templates_all_atom_masks[ idx ] = all_atom_masks[ idx ][0]
    output_templates_sequence[ idx ] = template_seq[ idx ]
    output_confidence_scores[ idx ] = confidence_scores[ idx ] # 0-9 where higher is more confident

  output_templates_sequence = ''.join(output_templates_sequence)

  templates_aatype = residue_constants.sequence_to_onehot(
      output_templates_sequence, residue_constants.HHBLITS_AA_TO_ID)

  template_feat_dict = {'template_all_atom_positions': np.array(templates_all_atom_positions)[None],
       'template_all_atom_masks': np.array(templates_all_atom_masks)[None],
       'template_sequence': [output_templates_sequence.encode()],
       'template_aatype': np.array(templates_aatype)[None],
       'template_confidence_scores': np.array(output_confidence_scores)[None],
       'template_domain_names': ['none'.encode()],
       'template_release_date': ["none".encode()]}

  return template_feat_dict, ret_all_atom_positions, ret_all_atom_mask

def get_final_dict(score_dict, string_dict):
    final_dict = OrderedDict()
    keys_score = [] if score_dict is None else list(score_dict)
    keys_string = [] if string_dict is None else list(string_dict)

    all_keys = keys_score + keys_string

    argsort = sorted(range(len(all_keys)), key=lambda x: all_keys[x])

    for idx in argsort:
        key = all_keys[idx]

        if ( idx < len(keys_score) ):
            final_dict[key] = "%8.3f"%(score_dict[key])
        else:
            final_dict[key] = string_dict[key]

    return final_dict

def add2scorefile(tag, scorefilename, write_header=False, score_dict=None):
    with open(scorefilename, "a") as f:
        add_to_score_file_open(tag, f, write_header, score_dict)

def add_to_score_file_open(tag, f, write_header=False, score_dict=None, string_dict=None):
    final_dict = get_final_dict( score_dict, string_dict )
    if ( write_header ):
        f.write("SCORE:     %s description\n"%(" ".join(final_dict.keys())))
    scores_string = " ".join(final_dict.values())
    f.write("SCORE:     %s        %s\n"%(scores_string, tag))

def register_output( outtag, start_time, interface1, interface2, prediction_result, scorefilename ):

  plddt_array = prediction_result['plddt']
  plddt = np.mean( plddt_array )
  plddt_1 = np.mean( plddt_array[interface1] )
  plddt_2 = np.mean( plddt_array[interface2] )

  pae = prediction_result['predicted_aligned_error']

  pae_interaction1 = np.mean( pae[interface1,:][:,interface2] )
  pae_interaction2 = np.mean( pae[interface2,:][:,interface1] )
  pae_1 = np.mean( pae[interface1,:][:,interface1] )
  pae_2 = np.mean( pae[interface2,:][:,interface2] )

  pae_interaction_total = ( pae_interaction1 + pae_interaction2 ) / 2
  time = timer() - start_time

  score_dict = {
          "plddt_total" : plddt,
          "plddt1" : plddt_1,
          "plddt2" : plddt_2,
          "pae_interaction1" : pae_interaction1,
          "pae_interaction2" : pae_interaction2,
          "pae1" : pae_1,
          "pae2" : pae_2,
          "pae_interaction" : pae_interaction_total,
          "time" : time
  }

  write_header=False
  if not os.path.isfile(scorefilename): write_header=True
  add2scorefile(outtag, scorefilename, write_header=write_header, score_dict=score_dict)
  
  print(score_dict)
  print( f"Tag: {outtag} reported success in {time} seconds" )

def combine_batches(tags, feature_dict_dict, initial_guess_dict):
    print( tags )
    allkeys = feature_dict_dict[tags[0]].keys()
    processed_feature_dict = {}
    for key in allkeys:
        processed_feature_dict[key] = jnp.stack([feature_dict_dict[tag][key] for tag in tags], axis=0)
    processed_initial_guess_dict = jnp.stack([initial_guess_dict[tag] for tag in tags], axis=0)
    return processed_feature_dict, processed_initial_guess_dict

def unpack_batches( tags, start, feature_dict_dict, prediction_result, sfd_out, scorefilename ):

    # First unpack the structures
    all_struct_module = prediction_result['structure_module']
    proteins = {} 
    for i in range(all_struct_module['final_atom_positions'].shape[0]):
        key = tags[i]

        proteins[key] = protein.Protein(
           aatype=feature_dict_dict[key]['aatype'][0],
           atom_positions=all_struct_module['final_atom_positions'][i,...],
           atom_mask=all_struct_module['final_atom_mask'][i,...],
           residue_index=feature_dict_dict[key]['residue_index'][0] + 1,
           b_factors=np.zeros_like(all_struct_module['final_atom_mask'][i,...]) )

    # Then unpack and compute confidence metrics
    confidence_metrics = {}
    for i in range(prediction_result['predicted_lddt']['logits'].shape[0]):
        key = tags[i]

        curr_metrics = {}
        curr_metrics['plddt'] = confidence.compute_plddt(
                prediction_result['predicted_lddt']['logits'][i,...])
        if 'predicted_aligned_error' in prediction_result:
            curr_metrics.update(confidence.compute_predicted_aligned_error(
                prediction_result['predicted_aligned_error']['logits'][i,...],
                prediction_result['predicted_aligned_error']['breaks'][i,...]))

        confidence_metrics[key] = curr_metrics

    for tag in tags:

        unrelaxed_pdb_lines = protein.to_pdb(proteins[tag])
        
        outtag = f'{tag}_af2pred'
        unrelaxed_pdb_path = f'{outtag}.pdb'
        with open(unrelaxed_pdb_path, 'w') as f: f.write(unrelaxed_pdb_lines)
        
        add2silent( outtag, unrelaxed_pdb_path, sfd_out ) 
        os.remove( unrelaxed_pdb_path )

        interface1 = np.zeros(confidence_metrics[tag]['plddt'].shape)
        for resn in range_str_to_list(args.interface1):
          interface1[resn - 1] = 1

        interface2 = np.zeros(confidence_metrics[tag]['plddt'].shape)
        for resn in range_str_to_list(args.interface2):
          interface2[resn - 1] = 1

        register_output( outtag, start, interface1.astype(bool), interface2.astype(bool), confidence_metrics[tag], scorefilename)

def range_str_to_list(str_in):
    ranges = str_in.split(',')
    list_out = []
    for r in ranges:
        if '-' in r:
            endpts = r.split('-')
            assert len(endpts) == 2

            list_out += list(range(int(endpts[0]),int(endpts[1])+1)) #inclusive ranges

        else:
            list_out += [int(r)]

    return list_out

def add2silent( tag, pdb, sfd_out ):
    pose = pose_from_file( pdb )

    struct = sfd_out.create_SilentStructOP()
    struct.fill_struct( pose, tag )
    sfd_out.add_structure( struct )
    sfd_out.write_silent_struct( struct, args.outname+".silent" )

def predict_structure(tags, feature_dict_dict, initial_guess_dict, sfd_out, scorefilename, random_seed=0):  
  """Predicts structure using AlphaFold for the given sequence."""

  start = timer()
  print(f"running {model_name}")
  model_runner.params = model_params
  
  processed_feature_dict, processed_initial_guess_dict = combine_batches(tags, feature_dict_dict, initial_guess_dict)

  prediction_result = jax.vmap(model_runner.apply, in_axes=(None,None,0,0))(model_runner.params,
          jax.random.PRNGKey(0), processed_feature_dict, processed_initial_guess_dict)

  unpack_batches( tags, start, feature_dict_dict, prediction_result, sfd_out, scorefilename ) 

  print( f"{len(tags)} predictions made in {timer() - start} seconds" )

def generate_feature_dict( query_seq, template_pdb, template_mask_str ):
  template_sequence = get_seq_from_pdb(template_pdb, slash_for_chainbreaks=False)

  template_mask = np.ones([len(template_sequence)]) #start by including everything, then exclude indicated residues
  if len(template_mask_str):
    for res_range in template_mask_str.split(','):
      if '-' in res_range:
        start = int(res_range.split('-')[0])
        end = int(res_range.split('-')[1])
        template_mask[start:end] = 0
      else:
        template_mask[int(res_range)] = 0

  initial_guess = af2_all_atom_from_struct( template_pdb, template_sequence, template_mask, query_seq, just_target=False )

  template_dict, all_atom_positions, all_atom_masks = template_from_struct(template_pdb, template_sequence, query_seq, template_mask)
  
  # Gather features
  feature_dict = {
      **pipeline.make_sequence_features(sequence=query_seq,
                                        description="none",
                                        num_res=len(query_seq)),
      **pipeline.make_msa_features(msas=[[query_seq]],
                                   deletion_matrices=[[[0]*len(query_seq)]]),
      **template_dict
  }

  return feature_dict, initial_guess 

def tag_buffer2features(input_seqs, tag_buffer, template_pdb, template_mask_str):
  
  feature_dict_dict = {}
  initial_guess_dict = {}

  for tag in tag_buffer:
    
    feature_dict, initial_guess = generate_feature_dict(input_seqs[tag], template_pdb, template_mask_str)
    feature_dict_dict[tag] = model_runner.process_features(feature_dict, random_seed=0)
    initial_guess_dict[tag] = initial_guess

  return feature_dict_dict, initial_guess_dict

# Checkpointing Functions

def record_checkpoint( tag_buffer, checkpoint_filename ):
    with open( checkpoint_filename, 'a' ) as f:
        for tag in tag_buffer:
            f.write( tag )
            f.write( '\n' )

def determine_finished_structs( checkpoint_filename ):
    done_set = set()
    if not os.path.isfile( checkpoint_filename ): return done_set

    with open( checkpoint_filename, 'r' ) as f:
        for line in f:
            done_set.add( line.strip() )

    return done_set

def fasta_dict(fasta):
  with open(fasta,'r') as f:
    lines = f.readlines()
    
    seqs = {}
    for l in range(0,len(lines)):
      if lines[l].startswith('>'):
        if l > 1:
          seqs[name] = curr_seq
        
        name = lines[l].replace('>','').replace('\n','')
        curr_seq = ''
      else:
        curr_seq += lines[l].replace('\n','')
            
    seqs[name] = curr_seq
        
  return seqs

# End Checkpointing Functions

################## Begin Main Function ##################

sfd_out = core.io.silent.SilentFileData( args.outname+".silent", False, False, "binary", core.io.silent.SilentFileOptions())

input_dict = fasta_dict(args.fasta)

alltags = list(input_dict.keys())

scorefilename = args.outname+".sc"

for i in range(0,len(alltags),args.batch_size):

    subtags = alltags[i:min(i+args.batch_size,len(alltags))]
    subdict = {}
    for tag in subtags:
        subdict[tag] = input_dict[tag]
    feature_dict_dict, initial_guess_dict = tag_buffer2features(subdict, subtags, args.template_pdb, args.exclude_template)
    predict_structure(subtags, feature_dict_dict, initial_guess_dict, sfd_out, scorefilename)

print('done predicting')
