
def get_all_from_file(infile):
  with open(infile, 'r') as file:
    for line in file:
      if line[0:3]=='AI:': 
         ai=line[5:-2]
      if line[0:3]=='NL:': 
         nl=line[5:-2]
      if line[0:4]=='DEG:': 
         deg=line[6:-2]
         break
  ai0=ai.split(", ")
  AI=[int(a) for a in ai0 ]
  nl0=nl.split(", ")
  NL=[int(a) for a in nl0 ]
  deg0=deg.split(", ")
  DEG=[int(a) for a in deg0 ]
  return AI,NL,DEG
  
def get_some_from_file(infile,what): #what must be list with no duplicates
  toget={'0': []}
  with open(infile, 'r') as file:
    for line in file:
      for w in what:
         if line[0:len(w)]==w: 
            getting=line[len(w)+2:-2]
            gotten=getting.split(", ")
            toget[w]=[int(a) for a in gotten ]
            what.remove(w)
      if not len(what): break
  return toget
  
def tot_all(m):
  if m==4: nlmax=32768
  if m==3: nlmax=128
  n=2^m
  Vnl=vector(ZZ,nlmax)
  Va=vector(ZZ,n//2+1)
  Vd=vector(ZZ,n)
  import os
  D=os.listdir()
  for file in D:
    if file[-4:]=='.txt':
          AI,NL,DEG=get_all_from_file(file)
          Vnl+=vector(ZZ,NL)
          Va+=vector(ZZ,AI)
          Vd+=vector(ZZ,DEG)        
  return Va,Vnl,Vd
  
def tot_one(what):
  As=[]
  import os
  D=os.listdir()
  for file in D:
    if file[-4:]=='.txt':
          A=get_some_from_file(file,[copy(what)])
          As+=[A[what]]       
  return sum(map( lambda a : vector(ZZ,a), As ))
 
  
  
def table_dist(L,  cols=14, dist='~'):
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
         print("~$p_{"+dist+"}(x)\\%$~ &"+b[:-1]+'\\\\'+" \\hline")
         print("~$\\#$~ &"+c[:-1]+'\\\\'+ "\\hline")
         a,b,c='','',''  
  u=cols-(count%cols)
  print("\\hline")  
  print("~$x$~ &"+ a[:-1]+'&'*u+'\\\\'+" \\hline")
  print("$p_{"+dist+"}(x)\\%$~ &"+b[:-1]+'&'*u+'\\\\'+" \\hline")
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
  L0=[]
  o=0
  for i in range(l):
     o+=L[i]
     #L0.append(o)
     if nsamp//2<o: 
      print(i,o, nsamp//2, i-1, o-L[i])
      median=i-1
      break      
  print("Expectation: ", expect,", mode: ", mode,", median: ", median)
  
