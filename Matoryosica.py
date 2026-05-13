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
            self.ARCstate = False
        else:
            pass

        if self.ARCstate == False:
           self.ARCvalue = self.value_assign(data)
           #print(f"{self.ARCvalue}debugだお")
           self.ARCstate = True
        else:
            pass   

        return self.ARCvalue   
        
def storage_time(simdata,memory,compare,fail_id,step,total_time,pileup_data):
    
    fail_count = 0


    if step == 0:#1週目
      max_list =[]
      max_list = max(row[0] for row in simdata)
      for i in range(len(simdata)):
        compare.append(max_list - simdata[i][0])
        pileup_data.append(simdata[i][0])
        if compare[i] > memory:
            fail_id.append(i)
            #compare[i] = 0 #this code is for the safety reason to apply Fidelity caluculation

      total_time =max_list
      print(f"{step}回目、積み立て{pileup_data}") #デバック     
            
    else:#2週目以降
      
      

      max_list = []
      max_list = max(pileup_data)  #identify max value of formerdata
   

      for i in range(len(simdata)):
          #print(f"{fail_id[fail_count]} and now value is {i}") #this code is for debug but may couse erro because of out of list
          if not fail_id or len(fail_id) < fail_count+1:
              
              pass
              
          elif i == fail_id[fail_count]:
            pileup_data[i] = simdata[i][0]+pileup_data[i] + 1 #sumarize total time for new data
            fail_count += 1
            #print("the confirm of the change of pileup_data")
          else:
              #print("somethingwrong")
              pass  
          
      max_resuccess_total_time = max(pileup_data)#identify max value of semarized total time

      if max_resuccess_total_time == max_list:#if there is no lager value than formerdata,there is only max value equal to the former max
        fail_count = 0
        for i in range(len(simdata)):
              if not fail_id or len(fail_id) < fail_count+1:
                pass
              elif i == fail_id[fail_count]:
                  compare[i] = max_list - pileup_data[i] 
              else:
                  #print("somethingwrong")
                  pass #this is effected former results    
      #tatal_time stays formertime

      elif  max_resuccess_total_time > max_list:
        for i in range(len(simdata)):
            compare[i] = max_resuccess_total_time - pileup_data[i] 

      else:
          print("something wrong")        
          
      total_time = max_resuccess_total_time


      fail_id = [] #2回目以降新しいfail_idを組み込む為のコード

      print(f"{step}回目、積み立て{pileup_data}") #デバック 

      for i in range(len(simdata)):
        if compare[i] > memory:
           fail_id.append(i)
           #compare[i]=0
           
                      





    return compare,fail_id,total_time,pileup_data


def Distillation(compare,memory,simdata,fail_id):
    original_time=[]
    fail_count = 0
    for i in range(len(simdata)):
        original_time.append(-memory*np.log((2*simdata[i][1]/((0.95)**2))-1))
        if not fail_id:
            simdata[i][1] = (1 - 0.05) * (1 - 0.05) * (1 - (1/2) * (1 - np.exp(-(original_time[i]) / memory)))
        else:    
            if i != fail_id[fail_count]:
                simdata[i][1] = (1 - 0.05) * (1 - 0.05) * (1 - (1/2) * (1 - np.exp(-(original_time[i]+compare[i]) / memory)))
            else:
                pass#ここでsimdataを初期化するか迷ったが、simulation関数で置き換えているのでやらなくて良い！
                

    
    return simdata

   
          

def run_simulation(segments,data,simdata,memory,compare,fail_id,step,totaltime,pileup_data):
    #simdataに対してデコヒーレンスを考えて上げる。デコヒーレンスに関してはcompareの部分をmax以外でかけて上げる。このとき、failした物はデコヒーレンスをかけたりかけなかったりするが、それに関してはデコヒーレンスの関数を作ってから考える！！

    count = 0
    fail_count = 0

    #print(f"debug={step}={len(segments)}")

    for seg in segments:
       if step == 0:
           #print("debug")
           flag=0
           simdata.append(seg.step(data,flag)) 
       else:
        if  not fail_id or len(fail_id) < fail_count+1:#consider if there is no fail
            #print("pass")
            flag =0
            pass

        else:
            #print(f"debugid{fail_id[fail_count]}and{count}")        
            if fail_id[fail_count] == count:
                flag = 1
                fail_count += 1
                simdata[count] = seg.step(data,flag) 
            else:
                flag =0        
        
           
       count += 1

    print(simdata)   
       


    compare,fail_id,totaltime,pileup_data = storage_time(simdata,memory,compare,fail_id,step,totaltime,pileup_data)


    simdata = Distillation(compare,memory,simdata,fail_id)



    

    if not fail_id:
        all_done = True
    else:
        all_done = False    

    
    return all_done,simdata,compare,fail_id,totaltime,pileup_data     

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

def bellswaping_Hop(Fidelity):
    step = 0
    for i in range(len(Fidelity)-1):
        if step == 0:
            Ftotal = Fidelity[i]*Fidelity[i+1]+((1-Fidelity[i])*(1-Fidelity[i+1]))/3
            step +=1
        else:
            Ftotal = Ftotal*Fidelity[i+1] + ((1-Ftotal)*(1-Fidelity[i+1]))/3

    return Ftotal            

def plotter_lneth(y_Fidelity,y_time,x_data):
    
  plt.figure(figsize=(10, 5))
  plt.plot(x_data, y_Fidelity, color='blue', marker='o', linestyle='None', label='Fidelity')
  plt.title('Quantum Entanglement Fidelity', fontsize=14, fontweight='bold') # 題名
  plt.xlabel('distance', fontsize=12) # 横軸ラベル
  plt.ylabel('Fidelity', fontsize=12)    # 縦軸ラベル
  plt.grid(True, linestyle='--', alpha=0.7) # グリッド線
  plt.legend(loc='best') # 凡例
  plt.show()

  plt.figure(figsize=(10, 5))
  plt.plot(x_data,y_time,color='blue',marker='o',linestyle = 'None',label='entangurument time')
  plt.title('Quantum Entanglement Time', fontsize=14, fontweight='bold') # 題名
  plt.xlabel('distance', fontsize=12) # 横軸ラベル
  plt.ylabel('Time', fontsize=12)    # 縦軸ラベル
  plt.grid(True, linestyle='--', alpha=0.7) # グリッド線
  plt.legend(loc='best') # 凡例
  plt.show()

def difference(mean_time,data2):
    epsilon = abs((data2-mean_time[1]) / data2) * 100

    print(f"誤差率: {epsilon:.4f} %")

def main_loop():
    
    
    #データベースのインポート
    # ユーザーのホームディレクトリ（C:/Users/ユーザー名）を自動取得
    home = os.path.expanduser("~")
    
    # ホームディレクトリ以下の相対パスを指定
    # 例：デスクトップの「research」フォルダにある場合
    relative_path = "Desktop\研究データ\simulation_database.npy"
    
    # パスを結合
    full_path = os.path.join(home, relative_path)

    try:
        data = np.load(full_path)
        print(f"✅ ローカルCドライブからロード完了: {full_path}")
    except FileNotFoundError:
        print(f"❌ ファイルが見つかりません。パスを確認してください: {full_path}")
        return
    
    #誤差率簡易版の為のmean_timeのインポート
    relative_path2 = "Desktop\研究データ\mean_time_segment2.npy"
    
    # パスを結合
    full_path2 = os.path.join(home, relative_path2)

    try:
        data2 = np.load(full_path2)
        print(f"✅ ローカルCドライブからロード完了: {full_path2}")
    except FileNotFoundError:
        print(f"❌ ファイルが見つかりません。パスを確認してください: {full_path2}")
        return
    
    
    
    


    attempts_input = input('How many attempts (Enter number, e.g., 100): ')
    if not attempts_input.isdigit():
        attempts = 100
    else:
        attempts = int(attempts_input)
    #---define parameter---
    n_ARC =  paramset_ARC()
    n_ELs =  paramset_ELs()
    memory = paramset_memory()
    Fidelity = []
    Total_Fidelity = []
    
    mean_Fidelity = []
    mean_time = []
    x_data = []
    
    theta = []

    
    relative_path3 = "Desktop\研究データ\95tau.npy"
    full_path3 = os.path.join(home, relative_path3)

    try:
        tau = np.load(full_path3)
        print(f"✅ ローカルCドライブからロード完了: {full_path3}")
    except FileNotFoundError:
        print(f"❌ ファイルが見つかりません。パスを確認してください: {full_path3}")
        return





          


    #----------------------

    for num_len in range(n_ARC):#データの取得
        ftheta = []
        time = []
        for attempt in range(attempts):
            step = 0
            compare = [] #consider as waiting time
            fail_id =[]
            segments = [RepeaterSegment(i, n_ELs) for i in range(num_len + 1) ]
            totaltime = []
            # pileup_data = [[] for i in range(num_len)]
            # simdata = [[] for i in range(num_len)]
            pileup_data = [] #this is singl data so this consider only time
            simdata = []
            print("succed")
            while True:
                all_complete,simdata,compare,fail_id,totaltime,pileup_data = run_simulation(segments,data,simdata,memory,compare,fail_id,step,totaltime,pileup_data)
                step += 1
                if all_complete:
                    time.append(totaltime)
                    for i in range(len(simdata)):
                        Fidelity.append(simdata[i][1])
                    if num_len+1 >= 2: 
                        Total_Fidelity.append(bellswaping_Hop(Fidelity))
                    else:
                        Total_Fidelity.append(simdata[i][1])        
                    
                    break

            
       
            dif = tau[num_len] - totaltime
            if dif >= 0:
                ftheta.append(1)
            elif dif < 0:
                ftheta.append(0)        

        theta.append(np.mean(ftheta))
            
        
        x_data.append((num_len+1)*20*n_ELs)


        mean_Fidelity.append(np.mean(Total_Fidelity))
        mean_time.append(np.mean(time))
        
       
        

    



    plotter_lneth(mean_Fidelity,mean_time,x_data)

    #誤差率簡易版
    difference(mean_time,data2)

    for i in range(num_len+1):
         print(f"θ{theta[i]}")



 
if __name__ == "__main__":
    main_loop()
    
    

