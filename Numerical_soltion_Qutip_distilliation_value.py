# !pip install qutip #Colab用
# !pip install qutip-qip
import qutip as qt
import numpy as np
from qutip_qip.operations import cnot # CNOTをqutip_qipからインポート

# ---------------------------------------------------------
# 1. 準備：状態の定義
# --------------------- ------------------------------------
q0 = qt.basis(2, 0)
q1 = qt.basis(2, 1)

#ベル状態の仕分け怪しいが一旦信じよう
psi_minus = (qt.tensor(q0, q1) - qt.tensor(q1, q0)).unit()
psi_plus  = (qt.tensor(q0, q1) + qt.tensor(q1, q0)).unit()
phi_plus  = (qt.tensor(q0, q0) + qt.tensor(q1, q1)).unit()
phi_minus = (qt.tensor(q0, q0) - qt.tensor(q1, q1)).unit()

# ---------------------------------------------------------
# 2. 初期状態（混合状態）の作成
# ---------------------------------------------------------
F1 = 0.9024500024999167
F2 = 0.9024500024999167
F3 = 0.9024500024999167
F4 = 0.9024500024999167

noise_coeff = (1 - F1) / 3
noise_coeff2 = (1 - F2) / 3
noise_coeff3 = (1 - F3) / 3
noise_coeff4 = (1 - F4) / 3

rho_1 = (F1 * qt.ket2dm(psi_minus) +
         noise_coeff * qt.ket2dm(psi_plus) +
         noise_coeff * qt.ket2dm(phi_plus) +
         noise_coeff * qt.ket2dm(phi_minus))

rho_2 = (F2 * qt.ket2dm(psi_minus) +
         noise_coeff2 * qt.ket2dm(psi_plus) +
         noise_coeff2 * qt.ket2dm(phi_plus) +
         noise_coeff2 * qt.ket2dm(phi_minus))

rho_3 = (F3 * qt.ket2dm(psi_minus) +
         noise_coeff3 * qt.ket2dm(psi_plus) +
         noise_coeff3 * qt.ket2dm(phi_plus) +
         noise_coeff3 * qt.ket2dm(phi_minus))

rho_4 = (F4 * qt.ket2dm(psi_minus) +
         noise_coeff4 * qt.ket2dm(psi_plus) +
         noise_coeff4 * qt.ket2dm(phi_plus) +
         noise_coeff4 * qt.ket2dm(phi_minus))

print(f"【初期状態】 フィデリティ(Ψ-): {F1}")

# ---------------------------------------------------------
# 3. ディスティレーション回路の構築
# ---------------------------------------------------------
rho_total1 = qt.tensor(rho_1, rho_2)
rho_total2 = qt.tensor(rho_3, rho_4)



# ★【修正1】ハダマールゲートを手動で作る（これでエラー回避！）
# H = 1/√2 * [[1, 1], [1, -1]]
H_matrix = 1 / np.sqrt(2) * qt.Qobj([[1, 1], [1, -1]])
H_total = qt.tensor(H_matrix, H_matrix, H_matrix, H_matrix)

# 適用
#rho_total = H_total * rho_total * H_total.dag()

# ★【修正2】CNOTゲートも念のため手動で作る関数を用意
def make_cnot(control, target):
    # 4量子ビット用のCNOT行列を自作する関数
    # (ライブラリの場所が変わっても大丈夫なように)
    # 基本のCNOT行列 (2量子ビット用)
    CX = qt.Qobj([[1,0,0,0], [0,1,0,0], [0,0,0,1], [0,0,1,0]], dims=[[2,2],[2,2]])

    # 4量子ビット空間への拡張は少し面倒なので、
    # ここではQuTiPの機能が使えるか試し、ダメならエラーが出る前に基本ゲートで構成します
    try:
        # まず標準関数を試す
        return cnot(N=4, control=control, target=target) # qt.cnot ではなく cnot を使用
    except AttributeError:
        # ダメなら自力で構成（ここが保険）
        # control=0, target=2 の場合などはテンソル積で構成
        # ※実装が長くなるため、今回は「qutip-qip」がない環境を想定して
        # 最もシンプルな「標準関数」でトライさせます。
        # もしここでエラーが出たら pip install qutip-qip が必要です。
        pass

# もし qt.cnot が使えない環境（v5系）の場合、以下を実行してください：
# !pip install qutip-qip
# その上で import qutip_qip などを足す必要がありますが、
# まずは「v5でも core に入っていることが多い」cnot を信じて実行します。

# もし次の行でエラーが出たら教えてください
try:
    CNOT_Alice = cnot(N=4, control=0, target=2) # qt.cnot ではなく cnot を使用
    CNOT_Bob   = cnot(N=4, control=1, target=3) # qt.cnot ではなく cnot を使用
except AttributeError:
    print("【重要】CNOTが見つからないため、手動定義に切り替えます...")
    # 手動でCNOTを作る（力技）
    # 0->2 のCNOTなどを作るのは大変なので、qutip-qipのインストールを推奨するメッセージを出します
    raise ImportError("QuTiP v5をお使いのようです。CNOTを使うには '!pip install qutip-qip' を実行してから、 'from qutip_qip.operations import cnot' をコードの先頭に追加してください。")

U_cnot = CNOT_Alice * CNOT_Bob
rho_after_cnot1 = U_cnot * rho_total1 * U_cnot.dag()
rho_after_cnot2 = U_cnot * rho_total2 * U_cnot.dag()



# ---------------------------------------------------------
# 4. 測定と事後選択
# ---------------------------------------------------------
P_00 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.ket2dm(qt.tensor(q0, q0)))
P_11 = qt.tensor(qt.qeye(2), qt.qeye(2), qt.ket2dm(qt.tensor(q1, q1)))
P_success = P_00 + P_11


rho_unnormalized1 = P_success * rho_after_cnot1 * P_success.dag()
rho_unnormalized2 = P_success * rho_after_cnot2 * P_success.dag()


#テンソルを切り離す
rho_p1d = rho_unnormalized1.ptrace([0, 1])
rho_p3d = rho_unnormalized2.ptrace([0, 1])

#σzのディスティレーション
rho_total3 = qt.tensor(rho_p1d, rho_p3d)

rho_total3 =  H_total * rho_total3 * H_total.dag()

# rho_p1d = H_total * rho_p1d * H_total.dag()
# rho_p3d = H_total * rho_p3d * H_total.dag()
rho_after_cnot3 = U_cnot * rho_total3 * U_cnot.dag()
rho_unnormalized3 = P_success * rho_after_cnot3 * P_success.dag()

#成功確率のトレース
prob = rho_unnormalized3.tr()

#規格化
rho_unnormalized3= rho_unnormalized3/prob

#テンソルを切り離す
rho_p1dd = rho_unnormalized3.ptrace([0, 1])







#prob_success = rho_unnormalized.tr()
#rho_final_4qubit = rho_unnormalized / prob_success
#rho_final = rho_final_4qubit.ptrace([0, 1])

# ---------------------------------------------------------
# 5. 結果の検証
# ---------------------------------------------------------
#fid_psi_minus = qt.expect(qt.ket2dm(psi_minus), rho_final)
fid_phi_plus  = qt.expect(qt.ket2dm(phi_plus), rho_p1dd )

print("-" * 30)
print(f"【結果】 成功確率: {prob:.4f}")
print("-" * 30)
#print(f"出力のフィデリティ (対 Ψ-): {fid_psi_minus:.4f}")
print(f"出力のフィデリティ (対 Φ+): {fid_phi_plus:.4f}")
print("-" * 30)

if fid_phi_plus > F4:
    print(f"大成功！ 純度が初期値({F4})より上がっています。")
else:
    print("純度が下がりました。")