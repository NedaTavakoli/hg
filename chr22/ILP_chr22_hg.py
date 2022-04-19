import argparse
import time
from gurobipy import *
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import os 
import scipy.sparse as sp


parser = argparse.ArgumentParser(description='ILP solution for variant picking')

parser.add_argument('--hap', type = str, default='col4_chr22_haplotypes_frq_all_2.txt',
                    help='List of haplotypes per SNP location for a single chromosome')

parser.add_argument('--k', type=int,  default=2504,
										help='Number of haplotypes per chromosome')

parser.add_argument('--n', type=int,  default=1064080,
										help='Number of SNP positions per chromosome')

parser.add_argument('--N', type=int,  default=2128158,
										help='Number of decision variables for ILP (from col4_chr22_haplotypes_frq_all_2.txt)')

parser.add_argument('--alpha', type=int,  default=150,
										help='Read Length')

parser.add_argument('--delta', type=int, default=15,
                    help=' Error threshold')


args, unparsed = parser.parse_known_args()

# k = 2504
N = 2128158
# n = 1064079
# alpha = 150
# delta = 15

alpha = args.alpha
delta = args.delta
k = args.k
N = args.N
n = args.n

################################

#  Compute C

def compute_C(n):
  #N = int(N/2)
  C = np.ones((n))
  return C

C = compute_C(n)

################################

with open('chr22_snp_positions.txt', "r") as file2:
    positions = file2.readlines()

    # Remove all the \n at the end of each line
    for i in range(len(positions)):
       positions[i] =  positions[i].rstrip(os.linesep)
    
################################

# Reading list of haplotypes per position
# Read all lines and save as a list
with open('col4_chr22_haplotypes_frq_all_2.txt', "r") as file1:
    FileasList = file1.readlines()

# Remove all the \n at the end of each line
for i in range(len(FileasList)):
  FileasList[i] = FileasList[i].rstrip(os.linesep)
  # Split values of the file on comma
  FileasList[i] = FileasList[i].split(',')

  ################################

def ILP_haplotype_resolved(c, k, n, N, delta, alpha, LogToConsole=True, TimeLimit=60):

    
    # starting time
    start = time.time()
    
    # Create a new mode
    #model = gp.Model(env=e)
    model = gp.Model()
    #model = Model()
    model.params.LogToConsole = LogToConsole
    model.params.TimeLimit = TimeLimit # seconds

    # Create variables
    x = model.addMVar(lb=0.0, ub=GRB.INFINITY,shape=n, vtype=GRB.BINARY, name="x")
    model.setObjective(c @ x, GRB.MAXIMIZE)
 

    # Add constraints
    #model.addConstr(A @ x <= B, name="c")

    lhs = gp.LinExpr(0)
    p = 0
    while(p < n):
        #print("p", p)
        if (FileasList[p] != ['']): # ignore empty line, I noticed that line 81 for chr22 is empty of SNPs (position 16053107 )
            lhs = gp.LinExpr(0)
            for index, value in enumerate(FileasList[p]):
                lhs = gp.LinExpr(0)
                if(value != " " and value !=","): # for line  1371 
                    lhs +=x[p]
                    j = p
                    while(j>=0 and (int(value)- int(positions[j])) < alpha):
                    # A[j, p] =1
                        #print("position", p,"index", index, "value", int(value))
                        lhs +=x[j]
                        j = j -1
                    
                    model.addConstr(lhs , GRB.LESS_EQUAL, 1.0 * delta);   
        p = p + 1

      
   
    # Optimize model
    model.optimize()

    total_profit = model.objVal
    items_selected_to_removed = [i for i in range(n) if x[i].x > 0.5]
    items_selected_to_retain= [i for i in range(n) if x[i].x < 0.5]

    
    # end time
    end = time.time()
    totol_time_ILP = end - start

    return items_selected_to_removed, items_selected_to_retain, total_profit, totol_time_ILP 
    #return total_profit, totol_time_ILP   
  ################################
items_selected_to_removed, items_selected_to_retain, total_profit,total_time_ILP = ILP_haplotype_resolved(C, k, n, N, delta, alpha, LogToConsole=True, TimeLimit=60)
print("\n")
print("total_profit", total_profit)
print("total time", total_time_ILP)
#print("list of removed edges", items_selected_to_removed)
#print("list of retained edges", items_selected_to_retain)
print("count of variants retained ", len(items_selected_to_retain))
