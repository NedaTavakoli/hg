import argparse
import time
from gurobipy import *
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import os 
import scipy.sparse as sp


parser = argparse.ArgumentParser(description='ILP solution for variant picking')

parser.add_argument('--hap', type = str, default='chr22_haplotypes_frq_all_strings.txt',
                    help=' A text file containg the list of haplotypes per SNP location for a single chromosome')

parser.add_argument('--pos', type = str, default='chr22_snp_positions_numbered.txt',
                    help=' A text file containg the list of SNP positions for a single chromosome')

parser.add_argument('--alleles', type = str, default='chr22_num_allels_per_positions_numbered.txt',
                    help=' A text file containg the list of number of alleles per SNP positions for a single chromosome') 

parser.add_argument('--samples', type = str, default='chr22_sample_dup_numbered.txt',
                    help=' A text file containg the list of samples per SNP location for a single chromosome')                                      

parser.add_argument('--chr', type = int, default=22,
                    help='Chromosome number (Default is 22)')                    

parser.add_argument('--k', type=int,  default=5008,
										help='Number of haplotypes per chromosome (default is 5008)')

parser.add_argument('--n', type=int,  default=1059517,
										help='Number of SNP positions per chromosome (Default is 1059517 for chr22)')

parser.add_argument('--N', type=int,  default=2123135,
										help='Number of decision variables for ILP (number of edges in variation graph, default is 2123135 for chr22)')

parser.add_argument('--alpha', type=int,  default=1000,
										help='Read Length (Default 1000)')

parser.add_argument('--delta', type=int, default=100,
                    help=' Error threshold (Default 100')

#****************************************************************************************************************
args, unparsed = parser.parse_known_args()

alpha = args.alpha
delta = args.delta
k = args.k
N = args.N
n = args.n
chr = args.chr
haplotypes = args.hap
positions = args.pos
alleles = args.alleles
samples = args.samples

#****************************************************************************************************************

# Reading list of haplotypes per position
# Read all lines and save as a list
#with open('chr1_haplotypes_frq_all_strings.txt', "r") as file1:
with open(haplotypes, "r") as file1:
    FileasList = file1.readlines()

#****************************************************************************************************************
# Convert haplotype names to numbers
# Numbers start at 1 
haplotype_name = {}  # this dictionary is used to save information of sample and numbers assign to it
#with open("chr22_sample_dup_numbered.txt") as file2:
with open(samples, "r") as file2:
 for line in file2:
 
    (key, value) = line.split()
 
    haplotype_name[int(key)] = value
#****************************************************************************************************************    
# Pos is a dictionary of rank for each positin and statrs from 1
POS = {}  # this dictionary is used to save information of sample and numbers assign to it
#with open("chr22_snp_positions_numbered.txt") as file3:
with open(positions, "r") as file3:

 for line in file3:
 
    (key, value) = line.split()
 
    POS[int(key)] = int(value)

#****************************************************************************************************************

N_alleles={}  # starts from 1
#with open("chr22_num_allels_per_positions_numbered.txt") as file4:
with open(alleles, "r") as file4:


  for line in file4:

    (key, value) = line.split()
    N_alleles[int(key)] = (int(value) -1)

#****************************************************************************************************************  
# defining corresponding x for alt alleles per position
# Starts from 1
x_pos_index = {new_list: [] for new_list in range(1,n+1)}   # starts from 1

t_prev = 0
for i in range(1,n+1):  # i from 0 to n
  t = N_alleles[i]
  for j in range(1,t+1): # repeat j times , j from 1 to t
    x_pos_index[i] += [j+ t_prev]
  t_prev += t

#****************************************************************************************************************
# Find all the positions located in range alpha from the current position
pos_in_range = {new_list: [] for new_list in range(1,n+1)}   # starts from 1

for p in range(1,n+1):  # SNP position, p starts from 1 to n
  j = p

  while(j>=0 and j<=n and (abs(POS[p]- POS[j])) < alpha):
    pos_in_range[p] +=[j]
    j = j +1
 
#****************************************************************************************************************

def ILP_haplotype_resolved( k, n, N, delta, alpha, LogToConsole=True):

    
    # starting time
    start = time.time()
    
    # Create a new mode
    model = gp.Model()

    model.params.LogToConsole = LogToConsole

    # Create variables
    x = model.addMVar(lb=0.0, ub=GRB.INFINITY,shape=N, vtype=GRB.BINARY, name="x")

    # set objective function
    obj = gp.LinExpr(0)
    for i in range(N):
      obj += x[i]

    # maximize c.x
    model.setObjective(obj, GRB.MAXIMIZE)


    # Add constraints
    #model.addConstr(A @ x <= B, name="c")

    # cols[0]: POS. 
    # cols[1]:N_allels, 
    # cols[2]:REF, 
    # cols[3]:ALT
    # cols[4]:List of haplotypes for REF
    # cols[5]:List of haplotypes for ALT1
    # cols[6]:List of haplotypes for ALT2
    # cols[7]:List of haplotypes for ALT3
    #********************************
    # Length of haplotypes for REF: len(cols[4].split(","))
    # Length of haplotypes for ALT1: len(cols[5].split(","))
    # Length of haplotypes for ALT2: len(cols[6].split(","))
    # Length of haplotypes for AL3: len(cols[6].split(","))
    #********************************
    # First haplotype for ALT2:  cols[6].split(",")[0]
    
    for p in range(1,n+1):  # SNP position, p starts from 1 to n
      #print("position:",p)
      for p_range in pos_in_range[p]:       # all positions in range, including itself
        cols = [x for x in FileasList[(p_range)-1].split()]
        n_x = int(cols[1]) - 1 #  Number of edges, to preserve the first haplotype we use minus 1

        # List of haplotypes   
        REF_l_hap = cols[4].split(",")
        ALT1_l_hap = cols[5].split(",")

        for h in range(1,k+1):  # For list of haplotypes, h starts from 1 to k
            
            hap_name = haplotype_name[h]
            lhs = gp.LinExpr(0)

            if hap_name in ALT1_l_hap:

                # Find x corrsponding to position
                l = x_pos_index[p_range][0]
                lhs +=x[l-1]

            try:
                ALT2_l_hap = cols[6].split(",")
                
                
        
                if hap_name in ALT2_l_hap:
                    
                    l = x_pos_index[p_range][1]
                    lhs +=x[l-1]  # because x starts from 0, but items starts from 1 

            except IndexError:
                pass 

            try:
                ALT3_l_hap = cols[7].split(",")
                
                if hap_name in ALT3_l_hap:
                    # Find x corrsponding to positions
                    l = x_pos_index[p_range][2]
                    lhs +=x[l-1]

            except IndexError:
                pass 

            model.addConstr(lhs , GRB.LESS_EQUAL, 1.0 * delta)
                
      
    # Optimize model
    model.optimize()

    total_profit = model.objVal
    items_selected_to_removed = [i for i in range(N) if x[i].x > 0.5]
    items_selected_to_retain= [i for i in range(N) if x[i].x < 0.5]

    
    # end time
    end = time.time()
    totol_time_ILP = end - start

    return items_selected_to_removed, items_selected_to_retain, total_profit, totol_time_ILP 
  #****************************************************************************************************************

items_selected_to_removed, items_selected_to_retain, total_profit,total_time_ILP = ILP_haplotype_resolved(k, n, N, delta, alpha, LogToConsole=True)
print("\n")
print("total_profit", total_profit)
print("total time", total_time_ILP)
#print("list of removed edges", items_selected_to_removed)
#print("list of retained edges", items_selected_to_retain)
print("count of variants retained ", len(items_selected_to_retain))
R = N-len(items_selected_to_retain)
print("size_of_graph_reduction", R)
print("ratio_of_graph_reduction", (R/N)*100)





