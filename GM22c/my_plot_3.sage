"""
Plotting results of statistics.

"""

import matplotlib.pyplot as plt

"""Global color name definition"""
orange=(1,0.5,0)  
green=(0,0.8,0) 
blue=(0,0,0.8)




def plot_dist(L,col=green,title='',fig_name='temp'):
  """Generate barplot from list.

  Args:
      L (list): list of values, namely L[i]/sum(L)% is the y-axis value of i.
      col (tuple, optional): color of bars. Defaults to green.
      title (str, optional): title of the graph. Defaults to ''.
      fig_name (str, optional): name of the outputfile (.eps automatically included). Defaults to 'temp'.
  """  
  V=vector(ZZ,L)
  l=len(L)
  nsamp=sum(V)
  Y=V/float(nsamp)*100
  plt.bar(range(l),Y,color=col)
  plt.title(title)
  #plt.xlim([30000, 40000])
  plt.ylabel('%')
  plt.xlabel('x')
  plt.savefig( fig_name+'.eps',format='eps')
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
  
  
    
