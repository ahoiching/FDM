import numpy as np
import cvxpy as cp

from tqdm import tqdm
from scipy.sparse import coo_matrix
from scipy.interpolate import interp1d

import cobra
import pandas as pd

import os, sys

import gurobipy
import warnings

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

current_objective_reactions=[]
current_objective_coefficients=[]
  
################################################################################
################################################################################
#----------------------Summary for flux decomposition--------------------------#
def Get_yield_flux_functional_gene_decomposition_from_file(file_path):
  FBA_model=cobra.io.load_json_model(file_path)

  (flux_modes, flux_modes_curate, flux_decomposition, functional_decomposition,
   gene_expressions_decomposition)=Get_yield_flux_functional_gene_decomposition(FBA_model)

  return (flux_modes, flux_modes_curate, flux_decomposition, functional_decomposition,
          gene_expressions_decomposition)

def Get_yield_flux_functional_gene_decomposition(FBA_model):
  ### 1. Load the FBA model to CVXpy
  (stoi_mat,
   reactions_dict,
   metabolites_ids,
   reactions_ids,
   reactions_lb,
   reactions_ub,
   linear_coefficients)=Cobra_cvxpy_conversion(FBA_model)

  gene_reactions_df=Get_gene_reaction_matrix(FBA_model)

  ### 2. Get the flux modes
  flux_modes=Run_flux_decomposition(stoi_mat,
                                  reactions_ids,
                                  linear_coefficients,
                                  reactions_lb, 
                                  reactions_ub)

  ### 3. Define the coupling matrix
  try:
    if ((flux_modes['__FBA_SOLUTIONS__']['AKGDH']>1e-6) and 
        (flux_modes['__FBA_SOLUTIONS__']['EX_o2_e']<-1e-6)):
        print('\n..\nThe condition is identified as aerobic E.coli..\n..\n..\n')
        coupling_matrix=np.diag(np.ones(len(flux_modes.columns)))
        coupling_matrix[flux_modes.columns=='ATPM']=(
          flux_modes.loc['AKGDH']/flux_modes['ATPM']['AKGDH'])
    elif ((flux_modes['__FBA_SOLUTIONS__']['EX_o2_e']==0) and 
         (flux_modes['ATPM']['EX_etoh_e']>0)):
        print('\n..\nThe condition is identified as anaerobic E.coli..\n..\n..\n')
        coupling_matrix=np.diag(np.ones(len(flux_modes.columns)))
        #coupling_matrix[flux_modes.columns=='ATPM']=(
        #flux_modes.loc['EX_etoh_e']/flux_modes['ATPM']['EX_etoh_e'])
        coupling_values=pd.concat([
                                    flux_modes.loc['EX_etoh_e']/flux_modes['ATPM']['EX_etoh_e'],
                                    flux_modes.loc['EX_ac_e']/flux_modes['ATPM']['EX_ac_e'],
                                    flux_modes.loc['EX_for_e']/flux_modes['ATPM']['EX_for_e']
                                  ],axis=1).min(axis=1)
        coupling_matrix[flux_modes.columns=='ATPM']=coupling_values
    else:
      warnings.warn("We cannot find an appropriate coupling matrix, either AKGDH (aerobic) or EX_etoh_e (anaerobic) is not carrying any fluxes. You may need to find your own coupling matrix to fix this")
      coupling_matrix=np.diag(np.ones(len(flux_modes.columns)))
  except:
    warnings.warn("We cannot find an appropriate coupling matrix. You may need to find your own coupling matrix to fix this")
    coupling_matrix=np.diag(np.ones(len(flux_modes.columns)))

  ### 4. Coupling the base fluxes (energy compensation), while decoupling the flux modes,
  ###    and get the flux decomposition
  (flux_modes_cmat,
   flux_decomposition_cmat)=Couple_flux_decomposition(flux_modes, coupling_matrix)

  ### 5. Get the functional decomposition
  functional_decomposition=Flux_to_functional_decomposition(flux_decomposition_cmat)

  ### 6. Get the gene expressions decomposition
  gene_expressions_decomposition = Reaction_to_gene_decomposition(flux_modes,
                                                                  functional_decomposition, 
                                                                  gene_reactions_df)

  return (flux_modes, flux_modes_cmat, flux_decomposition_cmat, functional_decomposition,
          gene_expressions_decomposition)

################################################################################
#------------------------Detailed function definitions-------------------------#
# The loading model function, input: path of the COBRApy FBA model in .json format.
def Load_json_model(file_path):
  FBA_model=cobra.io.load_json_model(file_path)
  (stoi_mat,
   reactions_dict,
   metabolites_ids,
   reactions_ids,
   reactions_lb,
   reactions_ub,
   linear_coefficients)=Cobra_cvxpy_conversion(FBA_model)

  return (FBA_model,
          stoi_mat,
          reactions_dict,
          metabolites_ids,
          reactions_ids,
          reactions_lb,
          reactions_ub,
          linear_coefficients)

def Cobra_cvxpy_conversion(FBA_model):
  metabolites_ids=[]
  metabolites_dict={}
  i=0
  for m in tqdm(FBA_model.metabolites):
    metabolites_ids.append(m.id)
    metabolites_dict[m]=i
    i+=1

  print('\n.\n.\nMetabolites loaded...')

  reactions_ids=[]
  reactions_dict={}
  #reactions_stoi=[]
  stoi_row=[]
  stoi_col=[]
  stoi_data=[]

  reactions_lb=[]
  reactions_ub=[]

  linear_coefficients=[]

  i=0
  for r in tqdm(FBA_model.reactions):
    reactions_ids.append(r.id)
    reactions_dict[r]=i
    reactions_lb.append(r.lower_bound)
    reactions_ub.append(r.upper_bound)
    linear_coefficients.append(r.objective_coefficient)
    #stoi_dict={}
    for m in r.metabolites:
      #stoi_dict[metabolites_dict[m]]=r.metabolites[m]
      stoi_col.append(i)
      stoi_row.append(metabolites_dict[m])
      stoi_data.append(r.metabolites[m])
    #reactions_stoi.append(stoi_dict)
    i+=1

  stoi_mat=coo_matrix((stoi_data,(stoi_row, stoi_col)))
  print('\n.\n.\nReactions loaded...')

  metabolites_ids=np.array(metabolites_ids)
  reactions_ids=np.array(reactions_ids)
  reactions_lb=np.array(reactions_lb)
  reactions_ub=np.array(reactions_ub)
  linear_coefficients=np.array(linear_coefficients)

  return (stoi_mat,
          reactions_dict,
          metabolites_ids,
          reactions_ids,
          reactions_lb,
          reactions_ub,
          linear_coefficients)

# The simulation function
# The parameters in this function would be:
# 1. The stoichiometry matrix.
# 2. List of reactions id.
# 3. Linear objective coefficient (usually it indicates which reaction flux need to be maximize).
# 4. Upper bounds and 5. lower bounds.
# The output of the function is a dictionary of the fluxes of each reaction.
def Simulate_ub_lb(stoi_mat,
                   reactions_ids,
                   linear_coefficients,
                   reactions_lb, 
                   reactions_ub,
                   verbose=False):
  x=cp.Variable(len(reactions_ids))
  prob = cp.Problem(cp.Maximize(linear_coefficients.T@x),
                  [stoi_mat @ x == 0, x<=reactions_ub, x>=reactions_lb])
  prob.solve(solver=cp.GUROBI, verbose=verbose,max_iters=1000)

  new_reactions_ub=np.array(reactions_ub)
  new_reactions_lb=np.array(reactions_lb)

  new_reactions_ub[linear_coefficients!=0]=x.value[linear_coefficients!=0]
  new_reactions_lb[linear_coefficients!=0]=x.value[linear_coefficients!=0]

  x_L2=cp.Variable(len(reactions_ids))
  fixed_index=(new_reactions_lb==new_reactions_ub)
  unfixed_index=(new_reactions_lb!=new_reactions_ub)
  intracelluar_index=np.array([(('EX_' not in reactions_ids[i]) or 
                                ('tpp' not in reactions_ids[i])) 
                                for i in range(len(reactions_ids))])

  prob2 = cp.Problem(cp.Minimize(
        cp.sum_squares(x_L2[unfixed_index & intracelluar_index])
        #cp.norm2(x_L2[unfixed_index & intracelluar_index])
        ),
                      [stoi_mat @ x_L2 == 0, 
                      new_reactions_ub[fixed_index]==x_L2[fixed_index],
                      new_reactions_lb<=x_L2,
                      new_reactions_ub>=x_L2])

  prob2.solve(solver=cp.GUROBI, 
                verbose=verbose
                ,BarConvTol=1e-12
                ,BarQCPConvTol=1e-12
                ,FeasibilityTol=1e-9
                ,OptimalityTol=1e-9)

  result_dict=dict(zip(reactions_ids,x_L2.value))

  return result_dict

def Run_flux_decomposition(stoi_mat,
                           reactions_ids,
                           linear_coefficients,
                           reactions_lb, 
                           reactions_ub):
  result_dict=Simulate_ub_lb(stoi_mat,
                             reactions_ids, 
                             linear_coefficients,
                             reactions_lb, 
                             reactions_ub)
  
  #eliminate the small fluxes
  abs_fluxes=np.abs(list(result_dict.values()))
  reactions_lb_fd=np.array(reactions_lb)
  reactions_lb_fd[abs_fluxes<1e-10]=0

  reactions_ub_fd=np.array(reactions_ub)
  reactions_ub_fd[abs_fluxes<1e-10]=0

  demand_fluxes=[]
  #Residuals
  der_reactions_id=[]
  #Performs derivatives on every reaction that is fixed.
  #If you only want to differentiate biomass vs. energy, you can only take the
  #derivative for the energy reactions (EX_ac_e and ATPM, for instance). This
  #only applies to when the biomass composition is fixed.
  sel_bool=((reactions_ub-reactions_lb<1e-8) &
                (abs(reactions_ub)>1e-8))
  
  print('Reactions that are constraint to carry fixed non-zero fluxes:')
  print('Lower bound\t\tUpper bound \t\tReaction ID ')
  for i in range(np.sum(sel_bool)):
    print(np.round(reactions_lb[sel_bool][i],10),
          ' \t\t',np.round(reactions_ub[sel_bool][i],10),
          '\t\t', reactions_ids[sel_bool][i])
    #Change base fluxes to demand fluxes.
    demand_fluxes.append(reactions_lb[sel_bool][i])
    der_reactions_id.append(reactions_ids[sel_bool][i])

  demand_fluxes=np.array(demand_fluxes)
  der_reactions_id=np.array(der_reactions_id)
   
  derivatives=demand_fluxes*1e-3
  demand_fluxes_mat=(np.array([demand_fluxes for i in range(len(demand_fluxes))])
                    +np.diag(derivatives))

  all_solution_series=[]
  for i in tqdm(range(len(der_reactions_id))):
    der_rxn=der_reactions_id[i]
    der=derivatives[i]
    reactions_lb_der=np.array(reactions_lb_fd)
    reactions_ub_der=np.array(reactions_ub_fd)
    #if der_rxn!='__FLUX_INTERCEPTS__':
    rxn_sel_bool=(reactions_ids==der_rxn)
    reactions_lb_der[rxn_sel_bool]+=der
    reactions_ub_der[rxn_sel_bool]+=der
    with HiddenPrints():
      result_dict_der=Simulate_ub_lb(stoi_mat,
                                    reactions_ids, 
                                    linear_coefficients,
                                    reactions_lb_der, reactions_ub_der)
      
    my_result_der=pd.DataFrame.from_dict(result_dict_der, 
                        orient='index')
    my_result_der.columns=[der_rxn]

    all_solution_series.append(my_result_der)

  all_solutions_pd=pd.concat(all_solution_series, axis=1)
  all_solutions_pd[np.abs(all_solutions_pd)<1e-8]=0

  fluxes_mat=all_solutions_pd.T.values

  solution=np.linalg.solve(demand_fluxes_mat, fluxes_mat)

  flux_modes = pd.concat([pd.DataFrame(demand_fluxes, index=all_solutions_pd.columns,
                          columns=['__DEMAND_FLUXES__']).T,
                          pd.DataFrame(solution, index=all_solutions_pd.columns,
                          columns=all_solutions_pd.index).T])
  
  result_dict['__DEMAND_FLUXES__']=0
  result_series=pd.Series(result_dict)
  result_series.name='__FBA_SOLUTIONS__'
  flux_modes['__FBA_SOLUTIONS__']=result_series
    
  return flux_modes

# The function "run_flux_decomposition_file" takes the 
# file path of an FBA model and return the yield table as a Pandas DataFrame
def Run_flux_decomposition_file(json_FBA_path):
  #json_FBA_path="iML1515_model_with_sinks.json"
  print('Loading FBA model:', json_FBA_path)
  (FBA_model,stoi_mat,reactions_dict,
   metabolites_ids,reactions_ids,
   reactions_lb,reactions_ub,
   linear_coefficients)=Load_json_model(json_FBA_path)

  flux_modes=Run_flux_decomposition(stoi_mat,
                                    reactions_ids,
                                    linear_coefficients,
                                    reactions_lb,
                                    reactions_ub)
  return flux_modes

def Run_flux_decomposition_model(FBA_model):
  (stoi_mat,reactions_dict,
   metabolites_ids,reactions_ids,
   reactions_lb,reactions_ub,
   linear_coefficients) = Cobra_cvxpy_conversion(FBA_model)

  flux_modes=Run_flux_decomposition(stoi_mat,
                                  reactions_ids,
                                  linear_coefficients,
                                  reactions_lb,
                                  reactions_ub)
  return flux_modes

def Get_flux_decomposition(flux_modes):
  flux_decomposition=flux_modes[1::]*flux_modes.loc['__DEMAND_FLUXES__']
  flux_decomposition=flux_decomposition[flux_decomposition.columns[0:-1]]
  return flux_decomposition

def Couple_flux_decomposition(flux_modes, coupling_matrix):
  #Rename "coupling_matrix" as "coupling_matrix"
  #Rename "flux_modes" as "flux_modes"
  inv_coupling_matrix=np.linalg.inv(coupling_matrix)
  flux_modes_cmat=np.matmul(flux_modes, inv_coupling_matrix)
  flux_modes_cmat.columns=flux_modes.columns
  flux_modes_cmat.loc['__DEMAND_FLUXES__']=np.matmul(coupling_matrix,
                                                flux_modes.loc['__DEMAND_FLUXES__'])
  flux_modes_cmat['__FBA_SOLUTIONS__']=flux_modes['__FBA_SOLUTIONS__']
  flux_decomposition_cmat=Get_flux_decomposition(flux_modes_cmat)
  flux_decomposition_cmat['__RESIDUALS__']=flux_modes['__FBA_SOLUTIONS__'
                                                    ]-flux_decomposition_cmat.sum(axis=1)

  return (flux_modes_cmat, flux_decomposition_cmat)

def Flux_to_functional_decomposition(flux_decomposition_cmat):
  reaction_directions=(flux_decomposition_cmat.sum(axis=1)/
                      np.abs(flux_decomposition_cmat.sum(axis=1)))

  Fmix = (1-(np.abs(flux_decomposition_cmat.sum(axis=1))/
            np.abs(flux_decomposition_cmat).sum(axis=1))).fillna(1)

  flux_decomposition_pos=flux_decomposition_cmat.multiply(reaction_directions, axis=0)
  flux_decomposition_pos[flux_decomposition_pos<0]=0

  functional_decomposition=flux_decomposition_pos.divide(
      flux_decomposition_pos.sum(axis=1), axis=0).fillna(0).multiply(
          1-Fmix, axis=0)

  functional_decomposition['__MIX_FUNCTIONS__']=Fmix

  return functional_decomposition

def Get_gene_reaction_matrix(FBA_model):
  gene_ids=[]
  gene_dict={}
  i=0
  for g in tqdm(FBA_model.genes):
    gene_ids.append(g.id)
    gene_dict[g]=i
    i+=1

  print('\n.\n.\nGetting gene reactions matrix from FBA model...')
  print('\n.\n.\nGenes loaded...')

  reactions_ids=[]
  reactions_dict={}
  gene_reaction_row=[]
  gene_reaction_col=[]
  gene_reaction_data=[]

  i=0
  j=0
  k=0
  for r in tqdm(FBA_model.reactions):
    if len(r.genes)!=0:
      reactions_ids.append(r.id)
      reactions_dict[r]=i
      for g in r.genes:
        gene_reaction_col.append(i)
        gene_reaction_row.append(gene_dict[g])
        gene_reaction_data.append(1)
      i+=1

  gene_reaction_mat=coo_matrix((gene_reaction_data,(gene_reaction_row, gene_reaction_col)))
  print('\n.\n.\nReactions loaded...')

  gene_ids=np.array(gene_ids)
  reactions_ids=np.array(reactions_ids)

  gene_reactions_df=pd.DataFrame(gene_reaction_mat.toarray()
                              ,index=gene_ids, columns=reactions_ids)

  return gene_reactions_df

def Reaction_to_gene_decomposition(flux_modes, functional_decomposition, gene_reactions_df):
  gene_reactions_fractions=gene_reactions_df.divide(gene_reactions_df.sum(axis=1), axis=0)

  functional_decomposition_flux_dist=functional_decomposition.multiply(
      np.abs(flux_modes['__FBA_SOLUTIONS__'][1::]), axis=0)

  gene_expressions_decomposition_flux_dist=(functional_decomposition_flux_dist.reindex(
                                    gene_reactions_fractions.columns
                                    ).T
                                  ).dot(gene_reactions_fractions.T).T

  gene_expressions_decomposition_flux_dist.loc[gene_expressions_decomposition_flux_dist.sum(axis=1)==0,
                                              '__RESIDUALS__']=1
                                              
  gene_expressions_decomposition=gene_expressions_decomposition_flux_dist.divide(
      gene_expressions_decomposition_flux_dist.sum(axis=1), axis=0)

  return gene_expressions_decomposition
