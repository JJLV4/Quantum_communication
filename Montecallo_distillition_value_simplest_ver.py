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

#from google.colab import drive
#drive.mount('/content/drive') # Excelを使用する場合はコメントアウトを外す2

# --- 1. 確率判定関数 ---
def check_success(prob):
    return random.random() < prob

# --- 2. セグメント管理クラス ---
class RepeaterSegment:
    def __init__(self, id, n_els):
        self.id = id
        self.n_els = int(n_els)
        self.max_buckets = max(1, self.n_els - 1)

        self.el_states = [False] * self.n_els
        self.ec_states = [False] * (self.n_els - 1)

        self.qstR_states = [False] * (self.n_els - 1)
        self.qstL_states = [False] * (self.n_els - 1)

        self.is_complete = False
        self.stored_in_afc = False
        self.current_bucket_idx = 0
        self.global_swap_done = False

        self.now_helald = 0
        self.helald_distance = 0
        self.helald_time = False

        self.storage_count =0
        self.step_count = 0

        self.qstR_fails = [False] * (self.n_els - 1)
        self.qstL_fails = [False] * (self.n_els - 1)

        self.tQR = 0
        self.tCNOT =0

        self.golobal_count = 0

        self.RP = False


    def reset_process(self):
        self.helald_time = True
        self.now_helald = 0 #この時点でヘラルディング信号は1進んでいると考えて良いかも
        failed_indices = [i for i, state in enumerate(self.ec_states) if not state]
        #if self.n_els ==1:
        self.helald_time =False
        self.RP =True
        #self.reset()#ここにリセットがあるから駄目説
        #else:
          #self.helald_distance = self.n_els




    # --- バケツの進行状況を数える関数 (クラス内のメソッドにしました) ---
    def count_active_buckets(self, states_list):
        """連続して埋まっているバケツの数を数える"""
        count = 0
        for state in states_list:
            if state:
                count += 1
            else:
                break
        return count

    def step(self, params,method):

        #print(f"Debag Qr{params["prob_Global_Swap"]}")
        #print(f"Debag AFC{params["prob_AFC"]}") # Fixed: Changed 'prob_afc' to 'prob_AFC'
        #print(f"Debag HOLD{params["prob_AFC_Hold"]}")
        #print(f"Debag EL{params["prob_EL_gen"]}")
        #print(f"Debag EC{params["prob_Internal_EC"]}")
        #print(f"debag{method}")


        if self.global_swap_done:
           self.storage_count+=1
           if self.storage_count == 10001 and method == "A":
              self.storage_count =0
              self.global_swap_done=False
              self.reset_process()

           elif self.storage_count == (params["separate"]*10000)+1  and method == "B": # Fix: prams -> params
              self.storage_count = 0
              self.global_swap_done=False
              self.reset_process()
           else:
            return

        if self.helald_time == True:

            self.now_helald+=1
            #print("debug")
            if not self.global_swap_done and method =="A":
                if self.el_states[0] and self.qstR_fails==True:
                        if len(self.qstL_states) > 0: # Check if qstL_states is not empty
                            # 現在どこまでバケツが進んでいるか確認
                            current_idx = self.count_active_buckets(self.qstL_states)

                            # まだゴールまで達していない場合
                            if current_idx < len(self.qstL_states):

                                # 次に進むための条件: 「その場所のECが成功していること(ヘラルディング信号)」
                                # current_idx 番目のECが成功していれば、バケツを1つ進められる
                                if self.ec_states[current_idx]:

                                    # メモリへの転送/保持確率 (eta_AFC または prob_AFC_Hold)
                                    if check_success(params["prob_AFC"]): # Fixed: Changed 'params.get("eta_AFC", 0.99)' to 'params["prob_AFC"]'
                                        # 成功！ バケツを埋める
                                        self.qstL_states[current_idx] = True
                                        self.current_bucket_idx += 1
                                    else:
                                      self.qstR_fails == False


                            else:
                                # バケツが最後まで埋まった状態（完了待機）
                                pass
                        elif self.n_els == 1 and self.is_complete: # Special case for single EL segments
                            pass # No buckets to fill, completion is already 'is_complete'



                    #右
                if self.el_states[self.n_els - 1] and self.qstL_fails==True:
                    if len(self.qstR_states) > 0: # Check if qstR_states is not empty
                        current_idx = self.count_active_buckets(self.qstR_states)

                        # まだゴールまで達していない場合
                        if current_idx < len(self.qstR_states):

                            # 次に進むための条件: 「その場所のECが成功していること(ヘラルディング信号)」
                            # current_idx 番目のECが成功していれば、バケツを1つ進められる
                            if self.ec_states[self.n_els - 2 - current_idx]:

                                # メモリへの転送/保持確率 (eta_AFC または prob_AFC_Hold)
                                if check_success(params["prob_AFC"]): # Fixed: Changed 'params.get("eta_AFC", 0.99)' to 'params["prob_AFC"]'
                                    # 成功！ バケツを埋める
                                    self.qstR_states[current_idx] = True
                                    self.current_bucket_idx += 1
                                else:
                                    self.qstR_fails==False



                        else:
                            # バケツが最後まで埋まった状態（完了待機）
                            pass
                    elif self.n_els == 1 and self.is_complete: # Special case for single EL segments
                        pass # No buckets to fill, completion is already 'is_complete'



            if self.now_helald == self.helald_distance:
                self.reset()
                self.helald_time = False


        else:

            if self.RP == True:
              self.reset()
              self.RP = False


                  # Check if the segment has at least one bucket (n_els - 1 > 0)
                # If n_els == 1, then (n_els - 1) is 0, and qstR_states, qstL_states are empty lists.
                # In this case, accessing [-1] would cause an IndexError.
            if (self.n_els - 1) > 0 and self.qstR_states[-1] == True and self.qstL_states[-1] == True and not self.global_swap_done:
                # Original completion condition for segments with multiple ELs (and thus, buckets)

                if check_success(params["prob_Global_Swap"]):
                    self.global_swap_done = True


                else:

                    self.reset_process()

                self.golobal_count +=1

            if self.n_els == 1 and self.is_complete and not self.global_swap_done: # seg.n_els == 1, meaning qstR_states and qstL_states are empty
                # For a segment with only 1 EL, its completion is based on its 'is_complete' status
                # (which means its single EL is established).

                    if check_success(params["prob_Global_Swap"]):
                        self.global_swap_done = True

                    else:
                        self.reset_process()

                    self.golobal_count +=1




            # --- Phase 1: 内部リンク構築 ---
            if not self.el_states[-1] and not self.el_states[0]:
                #print(f"試行回数{(params["R_EPPS"])*params["t_AFC"]}") # Removed dictionary access here
                if not self.is_complete:
                      for i in range(self.n_els):
                          if not self.el_states[i]:
                            if method == "A":
                              for s in range(int(params["eta_EPPS"]*params["R_EPPS"]*params["t_AFC"])):
                                    if check_success(params["prob_EL_gen"]):
                                        self.el_states[i] = True
                                        break
                            else:

                                expected_trials = (params["eta_EPPS"] * params["R_EPPS"] * params["t_AFC"]) / params["separate"]

                                if expected_trials < 1.0:
                                    # 期待値が1未満の場合（例: 0.5回）
                                    # 「expected_trialsの確率」で「1回試行する権利」を得る
                                    if check_success(expected_trials):
                                        if check_success(params["prob_EL_gen"]):
                                            self.el_states[i] = True
                                else:
                                  # 1. まず期待値（小数）を計算
                                  expected_value = (params["eta_EPPS"] * params["R_EPPS"] * params["t_AFC"]) / params["separate"]

                                  # 2. 整数部分（確実に行う回数）と小数部分（追加チャンスの確率）に分ける
                                  base_trials = int(expected_value)        # 例: 1.5 -> 1
                                  fractional_part = expected_value - base_trials  # 例: 1.5 - 1 = 0.5

                                  # 3. 今回のステップで行う総試行回数を決定
                                  num_trials = base_trials
                                  if check_success(fractional_part): # 残りの確率で +1回 チャンスを得る
                                      num_trials += 1

                                  # 4. 決定した回数分だけループを回す
                                  for s in range(num_trials):
                                      if check_success(params["prob_EL_gen"]):
                                          self.el_states[i] = True
                                          break

                      for i in range(self.n_els - 1):#ECを後で光らせたかったらこのfor文を最初に持って行く
                          if self.el_states[i] and self.el_states[i+1] and not self.ec_states[i]:

                                if check_success(params["prob_Internal_EC"]):
                                  self.ec_states[i] = True





                      if all(self.el_states) and all(self.ec_states):
                          self.is_complete = True

                      else:
                        self.reset_process()




            else:
                # --- Phase 2: バケツリレー ---
                # 1. まず左端のELが成功しているかチェック
                if method =="A":
                  if self.el_states[0]:
                      if len(self.qstL_states) > 0: # Check if qstL_states is not empty
                          # 現在どこまでバケツが進んでいるか確認
                          current_idx = self.count_active_buckets(self.qstL_states)

                          # まだゴールまで達していない場合
                          if current_idx < len(self.qstL_states):

                              # 次に進むための条件: 「その場所のECが成功していること(ヘラルディング信号)」
                              # current_idx 番目のECが成功していれば、バケツを1つ進められる
                              if self.ec_states[current_idx]:

                                  # メモリへの転送/保持確率 (eta_AFC または prob_AFC_Hold)
                                  if check_success(params["prob_AFC"]): # Fixed: Changed 'params.get("eta_AFC", 0.99)' to 'params["prob_AFC"]'
                                      # 成功！ バケツを埋める
                                      self.qstL_states[current_idx] = True

                                  else:
                                      # メモリ転送失敗 -> リセット
                                      self.qstL_fails == True
                                      self.reset_process()
                              else:
                                  # ECがまだ成功していない（信号が来ていない）
                                  # 「待てない」厳格な設定ならここでリセット
                                  self.reset_process()#current_idxでヘラルディング信号の伝達は進めているから、最大の数から(current_idx+1)を引くのがよい

                          else:
                              # バケツが最後まで埋まった状態（完了待機）
                              pass
                      elif self.n_els == 1 and self.is_complete: # Special case for single EL segments
                          pass # No buckets to fill, completion is already 'is_complete'
                      else:
                          self.reset_process() # No ELs or no buckets to progress, reset.


                  #右
                  if self.el_states[self.n_els - 1]:
                      if len(self.qstR_states) > 0: # Check if qstR_states is not empty
                          current_idx = self.count_active_buckets(self.qstR_states)

                          # まだゴールまで達していない場合
                          if current_idx < len(self.qstR_states):

                              # 次に進むための条件: 「その場所のECが成功していること(ヘラルディング信号)」
                              # current_idx 番目のECが成功していれば、バケツを1つ進められる
                              if self.ec_states[self.n_els - 2 - current_idx]:

                                  # メモリへの転送/保持確率 (eta_AFC または prob_AFC_Hold)
                                  if check_success(params["prob_AFC"]): # Fixed: Changed 'params.get("eta_AFC", 0.99)' to 'params["prob_AFC"]'
                                      # 成功！ バケツを埋める
                                      self.qstR_states[current_idx] = True
                                  else:
                                      # メモリ転送失敗 -> リセット
                                      self.qstR_fails == True
                                      self.reset_process()

                              else:
                                  # ECがまだ成功していない（信号が来ていない）
                                  # 「待てない」厳格な設定ならここでリセット
                                  self.reset_process()#current_idxでヘラルディング信号の伝達は進めているから、最大の数から(current_idx+1)を引くのがよい

                          else:
                              # バケツが最後まで埋まった状態（完了待機）
                              pass
                      elif self.n_els == 1 and self.is_complete: # Special case for single EL segments
                          pass # No buckets to fill, completion is already 'is_complete'
                      else:
                          self.reset_process() # No ELs or no buckets to progress, reset.

                  if not self.el_states[self.n_els - 1] or not self.el_states[0]:
                      self.reset_process()

                elif method =="B":#バケツリレーのステップ処理から
                  #左のバケツリレー処理

                  if self.el_states[0]:
                    if len(self.qstL_states) > 0: # Check if qstL_states is not empty
                        # 現在どこまでバケツが進んでいるか確認
                        current_idx = self.count_active_buckets(self.qstL_states)

                        # まだゴールまで達していない場合
                        if current_idx < len(self.qstL_states):

                            # 次に進むための条件: 「その場所のECが成功していること(ヘラルディング信号)」
                            # current_idx 番目のECが成功していれば、バケツを1つ進められる
                            if self.ec_states[current_idx]:
                               # 成功！ バケツを埋める
                               self.qstL_states[current_idx] = True
                            else:
                               self.reset_process()
                        else:
                            # バケツが最後まで埋まった状態（完了待機）
                            pass
                    elif self.n_els == 1 and self.is_complete: # Special case for single EL segments
                        pass # No buckets to fill, completion is already 'is_complete'
                    else:
                        self.reset_process() #current_idxでヘラルディング信号の伝達は進めているから、最大の数から(current_idx+1)を引くのがよい



                  #右のバケツリレー処理
                  if self.el_states[self.n_els - 1]:

                    if len(self.qstR_states) > 0: # Check if qstR_states is not empty
                      current_idx = self.count_active_buckets(self.qstR_states)

                      if current_idx < len(self.qstR_states):
                        if self.ec_states[self.n_els - 2 - current_idx]:
                          self.qstR_states[current_idx] = True
                        else:
                          self.reset_process()
                    elif self.n_els == 1 and self.is_complete: # Special case for single EL segments
                        pass # No buckets to fill, completion is already 'is_complete'
                    else:
                        self.reset_process()#current_idxでヘラルディング信号の伝達は進めているから、最大の数から(current_idx+1)を引くのがよい





    def reset(self):
        self.qstR_states = [False] * (self.n_els - 1)
        self.qstL_states = [False] * (self.n_els - 1)
        self.el_states = [False] * self.n_els
        self.ec_states = [False] * (self.n_els - 1)
        self.is_complete = False
        self.current_bucket_idx = 0
        self.global_swap_done = False

# --- 3. シミュレーション実行関数 ---
def run_simulation(segments, params,method):
    for seg in segments:
          seg.step(params,method)

    all_done = True
    for i in range(len(segments)):
        seg = segments[i]



        if not seg.global_swap_done:
            all_done = False
            break
        #else:
            #print(f"how many els{seg.n_els}")
    return all_done

# --- 4. 可視化関数 (バケツリレー可視化強化版) ---
def draw_network(segments, params, step, method):
    G = nx.Graph()
    pos = {}
    node_colors = []
    edges_to_draw = []

    # 色定義
    C_QR_MAIN = "#00FF00"; C_QR_WAIT = "#CCCCCC"
    C_BKT_CURRENT = "#00FFFF" # 現在地 (明るい水色)
    C_BKT_PASSED  = "#008888" # 通過済み (暗い水色)
    C_BKT_OFF     = "#EEEEEE" # 空

    C_EL_ON = "#FFFF00"; C_EL_OFF = "#333333"; C_PATH_DONE = "#FF4500"
    C_EC_NODE = "#FFA500"; C_EC_OFF = "#DDDDDD"

    x = 0; y = 0

    # スタートQR
    G.add_node("QR_Start"); pos["QR_Start"] = (x, y)
    if segments[0].global_swap_done:
        node_colors.append(C_QR_MAIN);
        x += 1.2
    else:
        node_colors.append(C_QR_WAIT);
        x += 1.2

    for i, seg in enumerate(segments):
        prev_node = "QR_Start" if i==0 else f"QR_{i}"

        # ★修正ポイント1: 変数を先に初期化しておく
        l_buckets = []
        r_buckets = []
        last_node = prev_node # Method BのときはこれがそのままELへの接続元になる

        # === 左QST (バケツ) ===
        if method == "A":
            if len(seg.qstL_states) > 0: # Only draw buckets if they exist
                for b in reversed(range(len(seg.qstL_states))):
                    b_id = f"QST_{i}_L_B{b}"
                    G.add_node(b_id); pos[b_id] = (x, y)

                    if seg.global_swap_done:
                        col = C_BKT_PASSED
                    elif seg.qstL_states[b]:
                        is_head = True
                        if b < len(seg.qstL_states) - 1 and seg.qstL_states[b+1]:
                            is_head = False
                        col = C_BKT_CURRENT if is_head else C_BKT_PASSED
                    else:
                        col = C_BKT_OFF

                    node_colors.append(col)
                    l_buckets.append(b_id)
                    x += 0.5

                edges_to_draw.append((prev_node, l_buckets[0], "#000000"))
                for k in range(len(l_buckets)-1):
                    edges_to_draw.append((l_buckets[k], l_buckets[k+1], "#555555"))
                last_node = l_buckets[-1]
            else:
                last_node = prev_node

            x += 0.5

        else:
             pass


        # === 内部構造 (ELとEC) ===
        for j in range(seg.n_els):
            target = None
            if j < seg.n_els - 1:
                id_ec = f"EC_{i}_{j}"
                G.add_node(id_ec); pos[id_ec] = (x + 1.0, y)
                is_ec = seg.ec_states[j] or seg.global_swap_done
                node_colors.append(C_EC_NODE if is_ec else C_EC_OFF)
                target = id_ec; x += 2.0
            else:
                x += 1.0

                if method == "A" and len(seg.qstR_states) > 0:
                    target = f"QST_{i}_R_B0"
                else:
                    target = f"QR_{i+1}"

            col = C_EL_OFF
            if seg.global_swap_done: col = C_PATH_DONE
            elif seg.el_states[j]: col = C_EL_ON

            if "QR" in target and target not in G.nodes():
                G.add_node(target); pos[target] = (x, y)
                node_colors.append(C_QR_MAIN if seg.global_swap_done else C_QR_WAIT)

            edges_to_draw.append((last_node, target, col))
            last_node = target

        x += 0.5

        # === 右QST (バケツ) ===
        if method == "A":
            start_x_r = x
            if len(seg.qstR_states) > 0:
                for b in range(len(seg.qstR_states)):
                    b_id = f"QST_{i}_R_B{b}"
                    G.add_node(b_id); pos[b_id] = (start_x_r, y)

                    if seg.global_swap_done:
                        col = C_BKT_PASSED
                    elif seg.qstR_states[b]:
                        col = C_BKT_PASSED
                    else:
                        col = C_BKT_OFF

                    node_colors.append(col)
                    r_buckets.append(b_id)
                    start_x_r += 0.5

                for k in range(len(r_buckets)-1):
                    edges_to_draw.append((r_buckets[k], r_buckets[k+1], "#555555"))

                x = start_x_r + 0.5

        # === 右QR ===
        id_qr = f"QR_{i+1}"
        if id_qr not in G.nodes():
            G.add_node(id_qr); pos[id_qr] = (x, y)
            node_colors.append(C_QR_MAIN if seg.global_swap_done else C_QR_WAIT)

        if len(r_buckets) > 0:
            edges_to_draw.append((r_buckets[-1], id_qr, "#000000"))

        x += 1.5

    plt.figure(figsize=(18, 5))
    plt.title(f"Method {method} Trace View - Step: {step}")

    sizes = []
    for n in G.nodes():
        if "QR" in str(n) or "Start" in str(n): sizes.append(800)
        elif "EC" in str(n): sizes.append(400)
        else: sizes.append(150)

    final_edge_colors = []
    for u, v, c in edges_to_draw:
        G.add_edge(u, v)
        final_edge_colors.append(c)

    nx.draw(G, pos, node_color=node_colors, node_size=sizes,
            edge_color=final_edge_colors, width=2, with_labels=False)

    labels = {n: n.split('_')[0] for n in G.nodes() if "QR" in str(n) or "Start" in str(n)}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight='bold')

    plt.axis('off'); plt.tight_layout(); plt.show()

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

    # 修正点: 1でも2でも計算ロジックに入るようにする (あるいはデモモードを明確に分ける)
    # ここでは「0以外なら計算する」というロジックにします
    if demo_choice == 1:

        # --- パラメータから確率を計算するロジック ---

        # 1. 必要なパラメータの取り出し
        eta_bsm = param_dict.get("eta_BSM", 0.9)
        alpha = param_dict.get("alpha", 0.2)

        param_dict["t_AFC"] = 0.0001
        t_AFC = param_dict.get("t_AFC")


        l = param_dict.get("l", 10)
        print(f"the lenth of l1 = {l}")
        l = param_dict.get("t_AFC")*200000 #lはt_AFCに依存すると考える。
        print(f"the lenth of l2 = {l}")

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

        p_link_pure = (1-(1-np.exp(-alpha*l)*eta_bsm*(eta_det)**required_measurements)**gammaf) * (prob_afc)**2
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
        tau = (1/omega_epps * np.log(1-(1-eps)**(1/big_N))/np.log(1-prob_arcgen_n)) + ttrans

        # ★デバッグ用: 計算された確率を表示する
        print(f"eta_interface:{eta_interface}")
        print(f"Calculated EL Prob: {p_link_pure:.5f}")
        print(f"Calculated AFC Prob: {prob_afc:.5f}")
        print(f"Calculated QR Prob: {eta_qr:.5f}")
        print(f"Calculated Conn Prob: {eta_conn:.5f}")
        print(tau)


        return [p_link_pure, prob_afc, eta_qr, eta_conn]

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
                plt.plot(xdata, ydata[dataset], 'x', label = dlabels[dataset])


        if dlabels == "EDRplot":
          plt.ylabel("EDR")
        elif dlabels == "EDR95plot":
          plt.ylabel("EDR95")
        else:
          plt.ylabel("STD")

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


def Distilation_caluculation(e,F_total,mode="normal"):
  x1 = (1 - 0.05) * (1 - 0.05) * ((1 - e))**F_total[0]
  y1 = (1 - 0.05) * (1 - 0.05) * ((1 - e))**F_total[1]
  z1 = (1 - 0.05) * (1 - 0.05) * ((1 - e))**F_total[2]
  q1 = (1 - 0.05) * (1 - 0.05) * ((1 - e))**F_total[3]

  # x2, y2, z2, q2 の計算
  x2 = (1 - x1) / 3
  y2 = (1 - y1) / 3
  z2 = (1 - z1) / 3
  q2 = (1 - q1) / 3

  # p1, p2, p3, p4 の計算
  p1 = (y1*x2 + x1*y2)*(q1*z2 + z1*q2) + 4*x2*y2*z2*q2
  p2 = 2*x2*y2*(q1*z2 + z1*q2) + 2*z2*q2*(y1*x2 + x1*y2)
  p3 = (x1*y1 + x2*y2)*(z1*q1 + z2*q2) + 4*x2*y2*z2*q2

  # 元のコードの記述通り (x1*y1 + x1*y2)
  p4 = 2*x2*y2*(z1*q1 + z2*q2) + 2*z2*q2*(x1*y1 + x1*y2)

  Distilation_value = p3 / (p1 + p2 + p3 + p4)

  if mode =="brentq":
    return Distilation_value - q1
  else:
    return Distilation_value





# --- 6. メインループ ---
def main_loop():

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
    prob_el, prob_afc, prob_qr, prob_ec = probs





    sim_params = {
        "num_segments": int(param_dict["big_N"]),
        "n_ELs": int(param_dict["n"]),
        "prob_EL_gen": prob_el,
        "prob_AFC": prob_afc,
        "prob_Global_Swap": (prob_qr)**2,
        "prob_Internal_EC": prob_ec,
        "prob_AFC_Hold": 1.00,
        "l" : int(param_dict["l"]),
        "R_EPPS": param_dict["R_EPPS"],
        "t_AFC": param_dict["t_AFC"],
        "separate": param_dict["separate"],
        "eta_EPPS": param_dict["eta_EPPS"]
    }




    print(f"\n--- Simulation Parameters ---")
    print(f"Config: {sim_params['num_segments']} Segments, {sim_params['n_ELs']} ELs")
    time.sleep(2)

    segments = [RepeaterSegment(i, sim_params["n_ELs"]) for i in range(sim_params["num_segments"]) ]
    step = 0
    success_count = 0



    options_list = ["アニメーション (1回のみ)", "統計解析 (多数回試行)"]
    mode_name, mode_idx = option_select(options_list, "Select Mode:")

    # ==========================================
    # モード1: アニメーション (詳細表示)
    # ==========================================
    if mode_idx == 1:
        # セグメント初期化
        segments = [RepeaterSegment(i, sim_params["n_ELs"]) for i in range(sim_params["num_segments"]) ]
        step = 0
        success_count = 0

        print(">>> Starting Animation Mode <<<")
        time.sleep(2)

        while True:
            step += 1
            all_complete = run_simulation(segments, sim_params,method_choice)

            #clear_output(wait=True)
            draw_network(segments, sim_params, step, method_choice)

            if all_complete:
                print("\n>>> NETWORK ESTABLISHED! <<<")
                break
            else:
                print(f"Step {step}: Simulating...")


            time.sleep(0.1)

        print(f"Success at Step: {step}")
        print(f"Estimated Time: {step  * param_dict.get("t_AFC"):.6f} sec")
        print(param_dict.get("t_AFC"))
    # ==========================================
    # モード2: 統計解析 (高速周回・グラフなし)
    # ==========================================
    elif mode_idx == 2:
          options_list = ["A", "B"]
          Sndmethod_choice, _ = option_select(options_list, "長さの方式:")
          # 試行回数の入力
          forth = 0
          prob = []
          mean_list = []
          tau_list = []
          x_list = []
          y_list = []
          x_list2 = []
          y_list2 = []
          y_list3 = []
          y_err_lower = []  # エラーバー下側
          y_err_upper = []  # エラーバー上側
          F_list = [] #フィデリティーのカウント
          F_listR = []#round ごとのF
          F_total = []#stepを加算したフィデリティーのカウント
          step_R = []#ラウンド毎のステップ数

          if Sndmethod_choice == "A":
            MaxELs_input = int(input('How many max ELs (Enter number, e.g., 100): '))
            vals=MaxELs_input

          else:
            MaxARCs_input = int(input('How many max ARCs (Enter number, e.g., 100): '))
            vals=MaxARCs_input

          attempts_input = input('How many attempts (Enter number, e.g., 100): ')
          if not attempts_input.isdigit():
              attempts = 100
          else:
              attempts = int(attempts_input)

          print(f"\n>>> Starting Statistics Mode ({attempts} trials) <<<")
          print("Running...", end="")

          # 結果を保存するリスト

          start_time = time.time()

          # --- 試行ループ ---
          for num_len in range(vals):
                execution_times = [] # かかった時間 (ステップ数 * 単位時間)
                step_counts = []     # かかったステップ数
                probabilty_check = [] #確率でチェック
                Fiderity_times = []

                if Sndmethod_choice == "A": # n (EL数) を増やす
                    sim_params["n_ELs"] = num_len + 1
                  # N (セグメント数) は固定
                else:
                    sim_params["num_segments"] = num_len +1

                for forth in range(4):    #4 denotes first round of distillation.


                  for attempt in range(attempts):
                      # 【重要】毎回セグメントを新品に作り直す (リセット)

                      # print(f"the number of n:{sim_params["n_ELs"]:.5f}") # Removed excessive print
                      # print(f"N:{sim_params["num_segments"]:.5f}") # Removed excessive print

                      segments = [RepeaterSegment(i, sim_params["n_ELs"]) for i in range(sim_params["num_segments"]) ]
                      step = 0



                        # 1回のシミュレーション

                      while True:
                          step += 1
                          all_complete = run_simulation(segments, sim_params,method_choice)


                          if all_complete:
                              # ---------------------------------------------------
                              # 1. 光を出していた時間 (Generation Time)
                              # ---------------------------------------------------

                              Fiderity_count = sum(seg.storage_count for seg in segments)+num_len #num_lenはELのメモリの分のカウントこれは、Nを増やす方式にしか対応していないことに注意


                              if method_choice == "A":
                                  t_generation = (step) * param_dict.get("t_AFC")#＋1はデッドタイム
                              else:
                                  t_generation = ((step) * param_dict.get("t_AFC")) / param_dict.get("separate")#＋1はデッドタイム





                              t_latency_total =param_dict.get("t_CNOT") + param_dict.get("t_QR")

                              # Sum the golobal_count from all segments
                              total_golobal_count = sum(seg.golobal_count for seg in segments)

                              t_elapsed = t_generation + t_latency_total #* total_golobal_count
                              # print(total_golobal_count) # Removed excessive print






                              Fiderity_times.append(Fiderity_count)

                              execution_times.append(t_elapsed)

                              if method_choice == "A":
                                ttrans = param_dict.get("t_QR", 0.0) + param_dict.get("t_CNOT", 0.0) + param_dict.get("t_AFC")*(sim_params["n_ELs"]-1)
                                eta_qst_total = prob_qr * (param_dict.get("eta_AFC") ** (sim_params["n_ELs"] - 1))

                              else :
                                ttrans = param_dict.get("t_QR", 0.0) + param_dict.get("t_CNOT", 0.0) + param_dict.get("t_AFC")
                                eta_qst_total = prob_qr

  # その修正した効率を使って tau を計算
  # 分母の log の中身: 1 - (eta_qst_total^2 * ...)
                              tau = ((1/(sim_params["R_EPPS"]*sim_params["eta_EPPS"])) * np.log(1-(1-param_dict.get("eps"))**(1/sim_params["num_segments"])) / np.log(1 - (eta_qst_total**2) * (prob_el**sim_params["n_ELs"]) * (prob_ec**(sim_params["n_ELs"]-1)))) * sim_params["separate"] + ttrans
                              #print(tau)
                              #tau = (1/0.95 + ((1+sim_params["n_ELs"])/0.05))* (param_dict.get("eta_AFC"))
                              #prob_el, prob_afc, prob_qr, prob_ec, tau = probs

                              # print(t_elapsed) # Removed excessive print

                              if tau-t_elapsed >= 0:
                                probabilty_check.append(1)
                              else:
                                probabilty_check.append(0)

                              step_counts.append(step)
                              break

                          # 無限ループ防止 (適当な上限)
                          if step > 1000000000:
                              step_counts.append(step) # 失敗扱い
                              break

                      # 進捗表示 (10%ごと)
                      if (attempt + 1) % (attempts // 10 + 1) == 0:
                          print(".", end="")





                  #------Distilliationの要素構築------#
                  mean_step = np.mean(step_counts)
                  step_R.append(mean_step)
                  mean_Fiderity = np.mean(Fiderity_times)#intendをここにしないと,attemptsと、for(4)に入らない
                  F_listR.append(mean_Fiderity)


                #if forth == 0:
                  #pass
                #else:
                  #for r in range(forth-1):
                    #F_listR[forth] += step_R[r+1]
                #F_total.append(F_listR[forth])


                              # F_listR[0] ～ [2] までを処理したいので range(3)
                for i in range(3):
                    # step_R の [i+1] から [3] までを合計して足す
                    # スライスは「最後の数字を含まない」ので、3まで入れたければ 4 と書く
                    F_listR[i] += sum(step_R[i+1 : 4])

                # F_total への追加も、ループを使わず一発で追加できます
                F_total.extend(F_listR)








          # --- 統計計算 (numpyを使用) ---
                mean_prob = np.mean (probabilty_check)

                tau_list.append(tau)

                mean_Fiderity = np.mean(Fiderity_times)
                mean_time = np.mean(execution_times)
                time_95_percentile = np.percentile(execution_times, 95)
                edr_95 = 1.0 / time_95_percentile
                std_dev_time = np.std(execution_times) # 標準偏差

                mean_step = np.mean(step_counts)
                std_dev_step = np.std(step_counts)

                time_high = mean_time + std_dev_time
                time_low = mean_time - std_dev_time

                print(mean_time)

                # 平均レート (1 / 平均時間)
                #edr = 1.0 / mean_time #平均時間で割らずにそれぞれのタイムで割ってedr求める


                #edr_pre = 1.0 / np.array(execution_times)
                #edr = np.mean(edr_pre)
                t_err = np.std(execution_times,ddof=1)/np.sqrt(attempts)


                #N =len(execution_times)
                #t_sem =  t_std / np.sqrt(N)#これは標準誤差に限る
                #edr_err = (edr **2) * t_sem
                #edr_err = (edr **2) * t_std * 1.96




                #edr_min = time_low / mean_time
                #edr_max = time_high / mean_time
                #edr_min = (1-std_dev_time/mean_time)/mean_time
                #edr_max = (1+std_dev_time/mean_time)/mean_time
                #if time_low > (mean_time * 0.1):
                    # 平均の10%よりは時間かかってるなら、普通に計算
                    #edr_max = 1.0 / time_low
                #else:
                    # 時間が0付近、またはマイナスになる場合は、エラーバーを「平均の2倍」程度で止める
                    # (これ以上伸ばすとグラフが見えなくなるため)
                    #edr_max = edr * 2.0

                # 4. エラーバーの長さ (中心からの距離)
                #err_down = edr - edr_min
                #err_up = edr_max - edr
                #edr=edr/10 #デバック
                mean_list.append(mean_time)
                x_list.append(sim_params["num_segments"]*sim_params["n_ELs"]*sim_params["l"])
                y_list.append(mean_time)
                y_list3.append(edr_95)
                x_list2.append(sim_params["num_segments"]*sim_params["n_ELs"]*sim_params["l"])
                y_list2.append(std_dev_step)
                #y_err_lower.append(err_down)
                #y_err_upper.append(err_up)
                prob.append( mean_prob)
                F_list.append(mean_Fiderity)





          end_time = time.time()
          print(" Done!")
          x_data = np.array(x_list)
          y_data = np.array(y_list)
          y_err = t_err
          dlab = ["EDRplot"] # ラベル





          # プロッター呼び出し (y_dataは配列にする)


          plt.errorbar(x_data, y_data, y_err, fmt='o', capsize=5,ecolor='red', color='blue', label='EDR with Time-STD Error')
          plt.plot(x_data, tau_list, color='black', marker='o', linestyle='None', label='LQUOM Analytical')

          plt.grid()
          plt.legend()
          plt.show()


          x_data2 = np.array(x_list2)
          y_data2 = np.array(y_list2)

          y_data3 = np.array(y_list3)

          dlab2 = ["STDplot"] # ラベル
          plt.figure() # <--- これで「新しい白紙」を用意する！
          plotter(x_data2, y_data2, xlabel="ARC-R distance (km)", dlabels=dlab2)
          dlab3 = ["EDR95plot"]
          plt.figure() # <--- これで「新しい白紙」を用意する！
          plotter(x_data, y_data3, xlabel="ARC-R distance (km)", dlabels=dlab3)


          #ヒストグラム
          print("\nDisplaying Histogram for the last distance setting...")
          plt.figure(figsize=(10, 5))
          plt.hist(step_counts, bins=20, color='skyblue', edgecolor='black', alpha=0.7)

          # 平均値のライン
          plt.axvline(mean_step, color='red', linestyle='dashed', linewidth=1.5, label=f'Mean: {mean_step:.1f}')

          plt.title(f"Distribution of Steps to Success (Last Trial)")
          plt.xlabel("Steps")
          plt.ylabel("Frequency")
          plt.legend()
          plt.grid(axis='y', alpha=0.5)
          plt.show()
          #Debag
          if method_choice == "A":
            print(f"自分の値{mean_list[0]}")
            print(f"理想値{(1/(prob_el*prob_qr**2))*param_dict['t_AFC']}")
          else:
            print(f"自分の値{mean_list[0]}")
            print(f"理想値{(1/(prob_ec*(prob_el**2)*(prob_qr**2)))*(param_dict['t_AFC']/2)}")
            print(mean_step)
            print((1/(prob_ec*(prob_el**2)*(prob_qr**2))))

          if method_choice == "A":
            print(f"自分の値{mean_list[0]}")
            print(f"理想値{(1/(prob_el*prob_qr**2))*param_dict['t_AFC']}")
            print(param_dict['t_AFC'])
            print(sim_params["t_AFC"])

          num = 0
          ydata4 = np.array(prob)
          plotter(x_data, ydata4, xlabel="ARC-R distance (km)", dlabels=dlab3)

          while True :

            print(f"{num+1}番目のθ関数{prob[num]}")


            num+=1
            if num==len(prob):
              break
          num = 0

          while True :


            if Sndmethod_choice == "B":
              print(f"{num+1}の時のトータルフィデリティーの平均値{F_list[num]} ここ0になるのなぜか一応確認すべき（恐らくインテンドを変えたことによるミス）")
            else:
              print(f"{num+1}の時のトータルフィデリティーの平均値{(num+1)}")

            num+=1
            if num==len(prob):
              break

    #-----------Fiderityの計算式(N=1,n=1の時のみ対応)------------#
          for confirmation in range(len(F_total)):
            #print(f"The individual Fiderity value of {confirmation+1} value {(1-0.05)*(1-0.05)*((1-0.0001))**F_total[confirmation]}") #正しいバージョン
            print(f"The individual Fiderity value of {confirmation+1} value {(1-0.05)*(1-0.05)*((1-0.00001))**F_total[confirmation]}")#理想的なバージョン

          max_value=max(F_total)
          min_value = min(F_total)


          try:
            mathematic_value=brentq(lambda e: Distilation_caluculation(e,F_total,mode="brentq"),1e-9,0.1)

            check_diff = Distilation_caluculation(mathematic_value, F_total)

            print(f"【結果】 F' = F4 となる損益分岐点が見つかりました！")
            print(f"  ε (epsilon) = {mathematic_value:.10f}")
            print(f"  (検算時の誤差: {check_diff:.2e})")

          except ValueError:
            print("【エラー】 指定した範囲 (1e-9 〜 0.1) では F' = F4 となる点が見つかりませんでした。")
            print("条件が厳しすぎるか、F_total の設定値を見直す必要があります。")

          F4_threshold = (1 - 0.05) * (1 - 0.05) * ((1 - mathematic_value))**min_value #ここで、数値解析を行い閾値を定義
          Fiderity = (1-0.05)*(1-0.05)*((1-mathematic_value))**max_value
          Fiderity2 = (1-0.05)*(1-0.05)*((1-mathematic_value))**min_value



          print(f"基準となる初期フィデリティ (F4相当): {F4_threshold:.6f}")
          print("-" * 50)

          # ---------------------------------------------------------
          # 2. ループ探索
          # ---------------------------------------------------------
          # 0.00001 から 0.0001 まで 0.00001 刻みで変化させる
          epsilon_values = np.arange(0.00001, 0.00011, 0.00001)

          results = []
          crossed_threshold = False

          print(f"{'e (epsilon)':<12} | {'Distilation Value':<20} | {'判定 (vs F4)'}")
          print("-" * 50)

          for e in epsilon_values:
              # --- 計算ロジック (提供された式) ---
              # x1, y1, z1, q1 の計算 (理想の値)
              Distilation_value=Distilation_caluculation(e,F_total)
              # 結果の保存
              results.append(Distilation_value)

          realistic = 0.00001
          Distilation_value=Distilation_caluculation(realistic,F_total)

          #Distilation_value = ((Fiderity)**2 + (1/9) *(1-Fiderity)**2)/((Fiderity)**2+(2/3)*Fiderity*(1-Fiderity)+(5/9)*(1-Fiderity)**2) #これは、簡略化したディスティレーションの値
          if Distilation_value > Fiderity2 :
            print(f"distillation is succeed!!!----> Distilated value is {Distilation_value}")
            print(f"In this case the decohirence value is assumed as {realistic}")
            print("These value is gained by the state of |φ+> so we have to apply σy gate")
          else:
            print(f"distillation is failed...----> Distilated value is {Distilation_value}")
            print(f"In this case the decohirence value is assumed as {realistic}")
            print("These value is gained by the state of |φ+> so we have to apply σy gate")

          results = np.array(results)
          threshold_line = np.full_like(epsilon_values, F4_threshold)
          print(f"The threshold value is ---->{mathematic_value}")
          print(f"so we must obtain memorytime more than {0.0001/mathematic_value}s")
          print("but I confused why this below plot is nonliner")

          #プロッター関数を使った場合これで、一つのグラフにまとめたいというか、コードの意味の深層理解をしたい所
          # plt.figure() # <--- これで「新しい白紙」を用意する！
          # plotter(epsilon_values, results, xlabel="Decoherence Value", dlabels="Distillation Result (F')")
          # plotter(epsilon_values, threshold_line, xlabel="Decoherence Value", dlabels="Threshold ")

          plt.plot(epsilon_values, results, label="Distillation Result (F')")

          # 2本目：閾値のライン（比較しやすいように linestyle="--" で点線にします）
          plt.plot(epsilon_values, threshold_line, label="Threshold (F4)", linestyle="--", color="red")

          # ---------------------------------------------------------
          # グラフの見た目を整える（ラベルや凡例）
          # ---------------------------------------------------------
          plt.xlabel("Decoherence Value (epsilon)")
          plt.ylabel("Fidelity")
          plt.title("Distillation vs Threshold")

          plt.legend()      # これを書くことで、右上に「どの線が何か」の凡例が出ます
          plt.grid(True)    # 縦横に薄い補助線（グリッド）を入れて見やすくします

          # 最後にグラフを表示！
          plt.show()





if __name__ == '__main__':
    main_loop()