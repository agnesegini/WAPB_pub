#!/usr/bin/python3
## SageMath version 9.5

##Modules
from sage.crypto.boolean_function import BooleanFunction
import  itertools 
prodit=itertools.product
from multiprocessing import Pool,cpu_count
from functools import partial
#value for parlallelisation
cpu=cpu_count()

from datetime import datetime
def now():
   return '{:%Y-%m-%d-%H:%M:%S}'.format(datetime.now())


#returns an iterator over the v length subsequences of range(n) (lex ordering)
def ithp(n,v):
  return itertools.combinations(range(n), v)

#loading functions
load("supp_prop.sage")# indication functions, supports, verification WPB,WAPB and SWAPB
load("nlk.sage")#computing NLk both sequential and parallel 
load("AI.sage")#compute algebraic immunity using the uppbound n/2 

load("cmr.sage")#building CMR functions
load("LM-fun.sage")#building LM functions
load("TL-fun.sage")#building TL functions

load("GM22b/constructions_GM22b.sage")#WAPB constuctions from the paper GM22b
load("GM22b/hybrid16.sage")#building h16 hybrid function
load("GM22b/data_collect.sage")#statistics/ printing tables GM22b

load("GM22c/construction_GM22c.sage")#Bounds  constuctions from the paper GM22c
load("GM22c/bounds_GM22c.sage")#Bounds  constuctions from the paper GM22c
load("GM22c/stat_NL.sage")#collecting data NL for GM22c
load("GM22c/my_plot_3.sage")#plotting base function for experiments GM22c

load("GM23/stat_all.sage")
load("GM23/uptosymm.sage")
load("GM23/get_all.sage")
load("GM23/porcelain.sage")

