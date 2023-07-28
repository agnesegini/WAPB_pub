#load("my_plot_3.sage")

def get_NL_from_file(infile):
  with open(infile, 'r') as file:
    for line in file:
      if line[0:3]=='NL:': 
       s=line[5:-2]
       break

  s0=s.split(", ")
  S=[int(a) for a in s0 ]
  return S
  
def tot_NL(m):
  if m==4: n=32768
  if m==3: n=128
  V=vector(ZZ,n)
  import os
  D=os.listdir()
  for file in D:
    if file[-1]=='t':
          V+=vector(ZZ,get_NL_from_file(file))      
  return V
  
"""  
def table(L,  cols=14):
  V=vector(ZZ,L)
  l=len(L)
  nsamp=sum(V)
  print("sample size: ",nsamp)
  Y=V/float(nsamp)*100
  a,b,c='','',''
  count=0  
  print("|c|"+"c|"*cols) 
  for x in range(l):

   if Y[x]!=0:
      a+=' ~$%d$~ &' %x
      b+=' ~$%.3f' %Y[x] +'$~ &'
      c+=' ~$%d$~ &' %V[x]
      count+=1 
      if count%cols==0:
         print("\\hline")  
         print("~$x$~ &"+ a[:-1]+'\\\\'+" \\hline")
         print("~$p_{\\DistNn{n}}(x)\\%$~ &"+b[:-1]+'\\\\'+" \\hline")
         print("~$\\#$~ &"+c[:-1]+'\\\\'+ "\\hline")
         a,b,c='','',''  
  u=cols-(count%cols)
  print("\\hline")  
  print("~$x$~ &"+ a[:-1]+'&'*u+'\\\\'+" \\hline")
  print("$p_{\\DistNn{n}}(x)\\%$~ &"+b[:-1]+'&'*u+'\\\\'+" \\hline")
  print("~$\\#$~ &"+c[:-1]+'&'*u+'\\\\'+ "\\hline")
  #return a,b,c
  
def propD(L):
  V=vector(ZZ,L)
  l=len(L)

  nsamp=sum(V)
  print("sample size: ",nsamp)
  Y=V/float(nsamp)
  expect=sum([i*Y[i] for i in range(l)])
  mode=list(L).index(max(L))
  print(expect, mode)

"""  
