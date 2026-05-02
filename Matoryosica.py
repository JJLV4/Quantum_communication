from numpy._core.fromnumeric import mean
import numpy as np
import pandas
from matplotlib import pyplot as plt
import networkx as nx
import random
import time
from IPython.display import clear_output, display
import os
from scipy.optimize import brentq
#for the perm function
import itertools
#qutip install
import qutip as qt
import numpy as np
from qutip_qip.operations import cnot # CNOTをqutip_qipからインポート
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter # ぼかし用のライブラリ



class RepeaterSegment:
    def __init__(self,id,n_els):
  

#ここからマトリョーシカ専用クラス
        self.ARCstate = False
        self.ARCvalue = []

    def value_assign(self,data):
        value = [] 
        n= random.randint(0,100000 - 1)
        value.append(data[n][0])
        value.append(data[n][1])
        return value


    def step(self,data,flag):
        if flag == 1:
            self.ARCstate == False
        else:
            pass

        if self.ARCstate == False:
           self.ARCvalue = self.value_assign(data)
           self.ARCstate = True
        else:
            pass   

        return self.ARCvalue   
        






    def reset(self):
        self.qstR_states = [False] * (self.n_els - 1)
        self.qstL_states = [False] * (self.n_els - 1)
        self.el_states = [False] * self.n_els
        self.ec_states = [False] * (self.n_els - 1)
        self.is_complete = False
        self.current_bucket_idx = 0
        self.global_swap_done = False         


def storage_time(simdata,memory,compare,fail_id,step):
    new_time=[[] for i in range(fail_id)]
    fail_count = 0


    if step == 0:#1週目
      max =[]
      max = max(row[0] for row in simdata)
      for i in range(len(simdata)):
        compare.append(max - simdata[i][0])
        if compare[i] > memory:
            fail_id.append(i)
            compare[i] = 0 
    else:#2週目以降
      max_next = []
      for i in range(len(fail_id)):
        new_time.append(simdata[fail_id[i]])  
      max_next = max(row[0] for row in new_time)
      for i in range(len(simdata)):
          if i == fail_id[fail_count]:
              compare[i] = max_next - new_time[i]
              fail_count+=1
          else:
              compare[i] += max_next
      
      fail_id = [] #2回目以降新しいfail_idを組み込む為のコード

      for i in range(len(simdata)):
        if compare[i] > memory:
           fail_id.append(i)
           compare[i]=0
           
                      





    return compare,fail_id


def Distillation(compare,memory,simdata,fail_id):
    original_time=[]
    fail_count = 0
    for i in range(len(simdata)):
        original_time.append(-memory*np.log((2*simdata[i][1]/((0.95)**2))-1))
        if i != fail_id[fail_count]:
            simdata[i][1] = (1 - 0.05) * (1 - 0.05) * (1 - (1/2) * (1 - np.exp(-(original_time[i]+compare[i]) / memory)))
        else:
            pass#ここでsimdataを初期化するか迷ったが、simulation関数で置き換えているのでやらなくて良い！
                

    
    return simdata

   
          

def run_simulation(segments,data,simdata,memory,compare,fail_id,step):
    #simdataに対してデコヒーレンスを考えて上げる。デコヒーレンスに関してはcompareの部分をmax以外でかけて上げる。このとき、failした物はデコヒーレンスをかけたりかけなかったりするが、それに関してはデコヒーレンスの関数を作ってから考える！！

    count = 0
    fail_count = 0

       
    for seg in segments:
       if step == 0:
           flag=0
           simdata.append(seg.step(data,flag)) 
       else:    
        if fail_id[fail_count] == count:
            flag = 1
            fail_count += 1
        else:
            flag =0        
        simdata[count] = seg.step(data,flag)     
       count += 1
       


    compare,fail_id = storage_time(simdata,memory,compare,fail_id,step)


    simdata = Distillation(compare,memory,simdata,fail_id)



    

    if not fail_id:
        all_done = True
    else:
        all_done = False    

    
    return all_done,simdata,compare,fail_id     

def paramset_ARC():
    MaxARCs_input = int(input('How many max ARCs (Enter number, e.g., 100): '))
    vals=MaxARCs_input
    return vals

def paramset_ELs():
    MaxELs_input = int(input('How many max ELs (Enter number, e.g., 100): '))
    vals=MaxELs_input
    return vals

def paramset_memory():
    memory_input = int(input('How long memory (Enter number, e.g., 100): '))
    vals=memory_input
    return vals

def main_loop():
    
    imp data #仮の変数として置いておく！Google Driveから後で持ってくる

    attempts_input = input('How many attempts (Enter number, e.g., 100): ')
    if not attempts_input.isdigit():
        attempts = 100
    else:
        attempts = int(attempts_input)
    #---define parameter---
    n_ARC =  paramset_ARC()
    n_ELs =  paramset_ELs()
    memory = paramset_memory()
    simdata = [[] for i in range(len(data))]
    Fidelity = []
    #----------------------

    for num_len in range(n_ARC):
        for attempt in range(attempts):
            step = 0
            compare = []
            fail_id =[]
            segments = [RepeaterSegment(i, n_ELs) for i in range(num_len) ]
            while True:
                all_complete,simdata,compare,fail_id = run_simulation(segments,data,simdata,memory,compare,fail_id,step)
                step += 1
                if all_complete:
                    time = max(row[0] for row in simdata)
                    for i in range(len(simdata)):
                        Fidelity.append(simdata[i][1])
                    #next task is incorporate bellswaping function    


                    break


