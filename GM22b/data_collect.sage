""" 
Functions to collect data and build tables used in [GM22b]
"""


def stat_fun_GM22b(g):
  """Given a Boolean function g prints a string including information about: 
      - number of variables (n:)
      - if g is SWAPB (SWAPB:)
      - the algebraic degree of g (degree:)
      - the algebraic immunity of g (AI:)
      - the nonlinearity of g (NL:)
      - all weightwise nonlinearity NLk for k=2..n-1 (NL_{k,n}:)
    
   Args:
      g (BooleanFunction): a Boolean function
  """      
  n=g.nvariables()
  print('n: ', n, ', SWAPB: ', is_SWAPB(g) , ', degree:',  g.algebraic_normal_form().degree(),', AI:', g.algebraic_immunity(),', NL: ', g.nonlinearity(),end=' ')
  for k in range(2,n-1):
    print( ", NL_{%d,%d}: " %(k,n), NLk_par(k,g), end=' ')
  print('\n') 
  

#print table gln in file name (must be a string)
def stat_gln_2_file_GM22b(n,name):
  """Creates a file including the NLk's of all the SWAPB functions g_{\ell,n} functions as
     defined in Definition 11 [GM22b]. 

  Args:
      n (int): 
      name (string): name output file
  """
  f=open(name, "w")
  for ell in range(2,n+1,2):
      print(ell,n, '\n')
      f.write(str(ell)+','+str(n)+ ':')
      g=construction_2_GM22b(n-ell,CMR(ell))
      #print  is_SWAPB(g) , g.algebraic_normal_form().degree(), g.algebraic_immunity(), g.nonlinearity(),
      for k in range(2,n-1):
        rr=NLk_par(k,g)
        print(rr, end=' , ')
        f.write(str(rr)+' , ')
      f.write('\n')
      print        
  f.close()
 
def table_1_GM22b(n):
   """Prints a tex-format table including the main cryptographic criteria of th  g_{\ell,n} functions as
     defined in Definition 11 [GM22b]. 
     For instance, results Table 1 and 2 [GM22b] can be reproduced by calling  table_1_GM22b(8)
     and table_1_file_GM22b(16), respectively. """
   
   for ell in range(2,n+1):
      print('$\const{%d}{%d}$' %(ell,n))
      g=construction_2_GM22b(n-ell,CMR(ell))
      L=[g.algebraic_normal_form().degree(), g.algebraic_immunity(), g.nonlinearity()] 
      for k in range(2,n-1):
        L+=[NLk_par(k,g)]
      print(*L, sep='~&~')
      print('\\\\')    
