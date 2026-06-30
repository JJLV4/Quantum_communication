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
import mpmath as mp
from mpmath import mp

#from google.colab import drive
#drive.mount('/content/drive') # Excelを使用する場合はコメントアウトを外す2


# --- 5. パラメータ読み込み関数群 ---
def load_params_dict(arch, arc_method, param_set=None):
        """
        「Parameters.xlsx」という名前のファイルを確認し、param_setで指定されたパラメータのセットをシミュレーションで使用するために辞書に読み込みます。
        パラメータセットが指定されていない状態で関数を呼び出すと、実装されているパラメータセットの一覧が返されます。
        """
        print(param_set)
        ##### Load presets from excel files
        if arch == "askarani":
           #df_data = pandas.read_excel(r'/content/drive/MyDrive/temp/Parameters_askarani.xlsx') Colab用
            df_data = pandas.read_excel(r"C:\Users\joe10\OneDrive\temp\Parameters_askarani.xlsx")

            if param_set == 1:
                column_key = "askarani near term"
            elif param_set == 2:
                column_key = "askarani long term"
            else:
                raise ValueError("The parameter set you indicated does not exist.")

            param_list = df_data[column_key].to_list()
            "df_data[colum_key]でExcelでdf_data = pandas.read_excel(r'/content/drive/MyDrive/temp/Parameters_askarani.xlsx'読み込んだデータから特定の列に選択して、.to_list()でデータの列を数値で上から順にリスト化する。"
            dictkeys = ["eta_NV", "t_NV", "t_13C", "t_CNOT", "eta_13C", "eta_QFC_1588", "Gamma_t", "eta_BSM", "eta_DET", "Gamma_f", \
                        "eta_AFC", "t_AFC", "l", "eta_EPPS", "R_EPPS", "eta_BUFF", "eta_MAP", "t_BUFF_spin", "eta_shift", "eta_pol", \
                        "eta_QFC_637", "alpha", "t_QR"]
            "dictkeysはdf_data[column_key].to_list()の値と変数の名前を対応関係にするための準備"
        elif arch == "ca":
            #df_data = pandas.read_excel(r'/content/drive/MyDrive/temp/Parameters_ca_trap.xlsx')#colab用
            df_data = pandas.read_excel(r"C:\Users\joe10\OneDrive\temp\Parameters_ca_trap.xlsx")
            
            if param_set == 1:
                column_key = "lquom"
            elif param_set == 2:
                column_key = "askarani_nearterm"
            else:
                raise ValueError("The parameter set you indicated does not exist.")

            param_list = df_data[column_key].to_list()

            # ★★★ ここにデバッグコードを追加！ ★★★
            print("\n" + "="*30)
            print("【デバッグ開始】Excelから読み込んだ生データ")
            print(f"データの個数: {len(param_list)}")
            print("-" * 20)

            # 全データをインデックス付きで表示
            for i, value in enumerate(param_list):
                # NaN (空白) はスキップせずに表示してみる
                print(f"Index {i}: {value}")

            print("="*30 + "\n")
            # ★★★ ここまで ★★★

            dictkeys = ["t_Ca43", "t_CNOT", "eta_BSM", "eta_DET", "Gamma_f", "eta_AFC", "t_AFC", "l", "eta_EPPS", "R_EPPS", "eta_BUFF", "t_buff", "eta_shift", "alpha" \
                        ,"eta_aom", "eta_interface", "t_QR"]
#etaは論文ではη
        loaded_dict = dict(zip(dictkeys, param_list))
        "dict(...): この関数は、zip が作ったペアの集まりを受け取り、各ペアの最初の要素をキー、2番目の要素を値とする辞書を作成します。"
        "zip(dictkeys, param_list): この関数は、先ほど準備したdictkeys リストと param_list リストの要素を、先頭から順番にペアにしていきます。"

#ここで言うmethodは動機方式によるものであり、長さを増やす方式ではないことに注意
        if arc_method == "A":
            # loaded_dict["n"] = np.arange(1,11,1)　  <-これが長さを増やす方式Bこれは論文のResultであるA,Bの長さの増やし方に関係しており、
            #今回はAを使いn=1でNを増やして行くやり方を採用している。Bでnを増やすやり方は、コメントアウトに追記していく
            loaded_dict["n"] = int(input('Enter the number of ELs per segment (n): '))
            # loaded_dict["big_N"] = 1　  <-長さを増やす方式B
            loaded_dict["big_N"] = int(input('Enter the number of ARCs per segment (N): '))
            loaded_dict["eps"] = 0.05
            loaded_dict["eta_loss"] = 0.9
            loaded_dict["separate"] = int(input('分割数を入力してください: '))

        elif arc_method == "B":
                      # loaded_dict["n"] = np.arange(1,11,1)　  <-これが長さを増やす方式Bこれは論文のResultであるA,Bの長さの増やし方に関係しており、
            #今回はAを使いn=1でNを増やして行くやり方を採用している。Bでnを増やすやり方は、コメントアウトに追記していく
            initial_n_els = int(input('Enter the number of ELs per segment (n): '))
            loaded_dict["big_N"] = int(input('Enter the number of ARCs per segment (N): '))
            loaded_dict["eps"] = 0.05
            loaded_dict["eta_loss"] = 0.9
            loaded_dict["separate"] = int(input('分割数を入力してください: '))
        ##### Load ARC variables depending on chosen arc_method
    #arcはarchitecher
                                    # Generic fiber loss coefficient

            ##### Method B makes some adjustments to ARC parameters, perform that adjustment here.
            # big_Nの計算には最初のnの値を使用
            loaded_dict["big_N"] = loaded_dict["big_N"] + loaded_dict["big_N"] * (initial_n_els - 1)
            print(f'===========================================> Check now we have big N equal to {loaded_dict["big_N"]} ')
            # loaded_dict["n"] = np.arange(2,7,1)
            loaded_dict["n"] = loaded_dict["separate"]
            #つまり、変数を増やさずに、nでxiを定義している。また今回はn=1の場合で始めARCを設定していたため、
            #分割した後はそれぞれのARCのn=1が分割されるためn=2となる。
            loaded_dict["l"] = loaded_dict["l"]/loaded_dict["n"]

        else:
            raise NotImplementedError

        return loaded_dict

def calculate_probabilities_from_params(demo_choice, param_dict, architecture, arc_method="A"):
    
    l = []
    mp.dps = 100

    # 修正点: 1でも2でも計算ロジックに入るようにする (あるいはデモモードを明確に分ける)
    # ここでは「0以外なら計算する」というロジックにします
    if demo_choice == 1:

        # --- パラメータから確率を計算するロジック ---

        # 1. 必要なパラメータの取り出し
        eta_bsm = param_dict.get("eta_BSM", 0.9)
        alpha = param_dict.get("alpha", 0.2)

        param_dict["t_AFC"] = 0.0001
        t_AFC = param_dict.get("t_AFC")


        # l = param_dict.get("l", 10)
        # print(f"the lenth of l1 = {l}")
        # l = param_dict.get("t_AFC")*200000 #lはt_AFCに依存すると考える。
        # print(f"the lenth of l2 = {l}")

        
        # 1. 10から800までの普通の配列を作る
        l = np.arange(10, 201, 10)

        # 2. mp.mpf を配列対応にパワーアップさせた「特製関数」を作る
        v_mpf = np.vectorize(mp.mpf)

        # 3. まとめて一撃で変換！（l_mp は超高精度数が詰まった配列になります）
        l_mp = v_mpf(l)


        eta_det = param_dict.get("eta_DET", 0.9)
        gammaf = param_dict.get("Gamma_f", 1)
        eta_EPPS = param_dict.get("eta_EPPS")
        R_EPPS = param_dict.get("R_EPPS")

        # 2. p_single_link (EL生成確率) の計算
        required_measurements = 2
        omega_epps = eta_EPPS*R_EPPS
        eps = param_dict.get("eps")
        n = param_dict.get("n")
        big_N = param_dict.get("big_N")
        tARC = t_AFC * (n-1)
        ttrans = param_dict.get("t_QR", 0.0)+param_dict.get("t_CNOT", 0.0)+tARC





        # 3. eta_afc (QST成功確率)
        prob_afc = param_dict.get("eta_AFC", 0.9)
        
        l = np.array(l)
        v_exp = np.vectorize(mp.exp)

        p_link_pure = (1-(1-v_exp(-0.046*l_mp)*eta_bsm*(eta_det)**required_measurements)**gammaf) * (prob_afc)**2
        #p_link_pure = (1-(1-mp.exp(-0.046*l_mp)*eta_bsm*(eta_det)**required_measurements)**gammaf) * (prob_afc)**2 #コンピュータの計算限界
        #p_link_pure = 1
        #p_link_pure = R_EPPS*eta_EPPS*(1-(1-np.exp(-alpha*l)*eta_bsm*(eta_det)**required_measurements)**gammaf) * (prob_afc)**2
        #p_link_pure = (1-(1-eta_EPPS*np.exp(-alpha*l)*eta_bsm*(eta_det)**required_measurements)**gammaf) * (prob_afc)**2
        #p_link_pure = (1-(1-np.exp(-alpha*l)*eta_bsm*(eta_det)**required_measurements)**(gammaf*eta_EPPS*R_EPPS*t_AFC)) * (prob_afc)**2

        # 4. eta_qr (QR成功確率)
        if architecture == "ca":
            eta_loss = param_dict.get("eta_loss", 0.9)
            eta_shift = param_dict.get("eta_shift", 1.0)
            eta_aom = param_dict.get("eta_aom", 0.95)

            eta_interface = eta_aom * eta_bsm * (eta_det)**required_measurements
            #eta_interface = param_dict.get("eta_interface", 0.3)
            #eta_qr = eta_loss * eta_shift * eta_aom * eta_interface
            #eta_qr = 1
            eta_qr = eta_loss * eta_shift * eta_aom * 0.3

        else:
            eta_qr = 0.9

        # 5. eta_conn (内部EC成功確率)
        if architecture == "ca":
            eta_conn = eta_bsm * (param_dict.get("eta_shift", 1.0)) * (eta_det)**required_measurements
            #eta_conn =1
        else:
            eta_conn = eta_bsm
        prob_arcgen_n = (eta_qr)**2* (p_link_pure)** n * (eta_conn)**(n-1)
        #tau = (1/omega_epps * np.log(1-(1-eps)**(1/big_N))/np.log(1-prob_arcgen_n)) + ttrans #計算限界により除外

        # ★デバッグ用: 計算された確率を表示する
        print(f"eta_interface:{eta_interface}")
        print(f"Calculated EL Prob: {p_link_pure}")
        print(f"Calculated AFC Prob: {prob_afc:.5f}")
        print(f"Calculated QR Prob: {eta_qr:.5f}")
        print(f"Calculated Conn Prob: {eta_conn:.5f}")
        #print(tau)


        return [p_link_pure, prob_afc, eta_qr, eta_conn,l]

    else:
        # デモ用確率 (demo_choice が 0 とかの場合のみここに来るようにする)
        print("Using Fixed Demo Probabilities (0.6)")
        return [0.6, 0.95, 0.95, 0.9]


def plotter(xdata, ydata, xlabel, dlabels):
        """
       シーケンスのプロットを処理する汎用関数です。入力の ydata はリストにすることもでき、
       その場合関数は複数の曲線をプロットします。Labels の入力は ydata の各セットごとに1件必要です。
       常にリストである必要があります。
        """

        if ydata.ndim == 1:
            plt.plot(xdata, ydata,'x', label = dlabels[0])

        else:
            for dataset in range(ydata.shape[0]):#nの分割数を複数で見る場合の動作
           # ydata.shape: 配列の形状（(行数, 列数)）をタプルで返します。ydata.shape[0]: タプルの最初の要素、つまり行数（＝データ系列の数）を取得します。
                plt.plot(xdata[dataset], ydata[dataset], 'x', label = dlabels[dataset])


        if dlabels == "EDRplot":
          plt.ylabel("EDR")
        elif dlabels == "EDR95plot":
          plt.ylabel("EDR95")
        else:
          plt.ylabel("prob whith cahnging l")

        plt.yscale("log")
        plt.grid()
        plt.legend()#label を指定した凡例をグラフ上に表示します(dlabelsの表示)
        plt.xlabel(xlabel)
        plt.show()

def option_select(options, input_message):
    print(input_message)
    for index, item in enumerate(options): print(f'{index+1}) {item}')
    user_input = input('Your choice (Enter number): ')
    if not user_input.isdigit() or int(user_input) not in range(1, len(options) + 1):
        return options[0], 1
    return options[int(user_input) - 1], int(user_input)


def main_loop():
   
    total_prob = []


    options_list = ["askarani", "ca"]
    architecture, _ = option_select(options_list, "Arch:")
    options_list = ["A", "B"]
    method_choice, _ = option_select(options_list, "同期方式:")
    if architecture =="A":
       options_list = ["askarani near term", "askarani long term"]
    else:
       options_list = ["LQUOM parameter", "askarani near term"]

    param_choice, param_index = option_select(options_list, "Params:")



    options_list =["実際のパラメータ","動作確認用のパラメータ"]
    # 「_」で名前を受け流し、番号だけを demo_choice に入れる
    _, demo_choice = option_select(options_list, "現実かDemo:")
    print(f"\nInit: {architecture}, {method_choice}, {param_choice}")
    param_dict = load_params_dict(architecture, method_choice, param_set=param_index)
    probs = calculate_probabilities_from_params(demo_choice,param_dict, architecture, method_choice)
    prob_el, prob_afc, prob_qr, prob_ec,l = probs





    sim_params = {
        "num_segments": int(param_dict["big_N"]),
        "n_ELs": int(param_dict["n"]),
        "prob_EL_gen": prob_el,
        "prob_AFC": prob_afc,
        "prob_Global_Swap": (prob_qr)**2,
        "prob_Internal_EC": prob_ec,
        "prob_AFC_Hold": 1.00,
        "R_EPPS": param_dict["R_EPPS"],
        "t_AFC": param_dict["t_AFC"],
        "separate": param_dict["separate"],
        "eta_EPPS": param_dict["eta_EPPS"]
    }

    t_latency_total = param_dict.get("t_CNOT") + param_dict.get("t_QR")

    for i in range(len(l)):
        N_tt = 200/l[i]#the maximum lenth is described as this place and param dfine place
        total_prob.append((1-(1-prob_el[i]*(prob_qr**2))**(sim_params["eta_EPPS"]*sim_params["R_EPPS"]*(1-t_latency_total-((1.47*l[i])/299792.458))))**N_tt)
        print((1-(1-prob_el[i]*(prob_qr**2))**(sim_params["eta_EPPS"]*sim_params["R_EPPS"]*(1-t_latency_total-((1.47*l[i])/299792.458))))**N_tt)
        print(total_prob[i])



    total_prob = np.array(total_prob)

    #plotter(l,total_prob,"distance","The change of the prob with l")

    plt.figure()
    #plt.ylim(1e-10000, 1e-19)
    #plt.yscale("log")
    plt.plot(l, total_prob, 'x-', label="Total Probability") # 点を 'x'、線で繋ぐ
    plt.plot(l,prob_el,'x', label="EL Probability")
    plt.xlabel("distance")
    plt.ylabel("The change of the prob with l")
    plt.title("Probability vs Distance")
    plt.legend()
    plt.grid(True)
    plt.show()



    

 

if __name__ == '__main__':
    main_loop()  