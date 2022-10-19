##Collect data and build tables    

#print  profile of g
def stat_fun(g):
      n=g.nvariables()
      print('n: ', n, ', SWAPB: ', is_SWAPB(g) , ', degree',  g.algebraic_normal_form().degree(),', AI', g.algebraic_immunity(),', NL', g.nonlinearity(),end=' ')
      for k in range(2,n-1):
        print( ", NL_{%d,%d}: " %(k,n), NLk_par(k,g), end=' ')
      print('\n') 
  

#print table gln in file name (must be a string)
def stat_gln_2_file(n,name):
  f=open(name, "w")
  for ell in range(2,n+1,2):
      print(ell,n, '\n')
      f.write(str(ell)+','+str(n)+ ':')
      g=construction_2(n-ell,CMR(ell))
      #print  is_SWAPB(g) , g.algebraic_normal_form().degree(), g.algebraic_immunity(), g.nonlinearity(),
      for k in range(2,n-1):
        rr=NLk_par(k,g)
        print(rr, end=' , ')
        f.write(str(rr)+' , ')
      f.write('\n')
      print        
  f.close()

#print table gln  
def table_1(n):
   for ell in range(2,n+1):
      print('$\const{%d}{%d}$' %(ell,n))
      g=construction_2(n-ell,CMR(ell))
      L=[g.algebraic_normal_form().degree(), g.algebraic_immunity(), g.nonlinearity()] 
      for k in range(2,n-1):
        L+=[NLk_par(k,g)]
      print(*L, sep='~&~')
      print('\\\\')    
