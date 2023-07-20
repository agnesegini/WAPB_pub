orange=(1,0.5,0)  
green=(0,0.8,0) 
blue=(0,0,0.8)



import matplotlib.pyplot as plt
def plot_dist(L,col=green,title='',fig_name='temp.eps'):
  V=vector(ZZ,L)
  l=len(L)
  nsamp=sum(V)
  Y=V/float(nsamp)*100
  plt.bar(range(l),Y,color=col)
  plt.title(title)
  #plt.xlim([30000, 40000])
  plt.ylabel('%')
  plt.xlabel('x')
  plt.savefig( fig_name,format='eps')
  plt.close()  
  

def plot_dist_NL(L,col=green,title='',fig_name='temp.eps'):
  V=vector(ZZ,L)
  l=len(L)
  nsamp=sum(V)
  Y=V/float(nsamp)*100
  plt.bar(range(l),Y,color=col)
  plt.title(title)
  mi=0
  while L[mi]==0: mi+=1
  plt.xlim([mi-l//100, l])
  plt.ylabel('%')
  plt.xlabel('x')
  plt.savefig( fig_name,format='eps')
  plt.close()    
  
  
    
