reset
set nokey                 # 凡例の非表示
set term gif animate      # 出力をgifアニメに設定
set output "map.gif"  # 出力ファイル名の設定

#-------------------------------------------------------------------------------
# 変数の設定
#-------------------------------------------------------------------------------
n0 = 0    # ループ変数の初期値
n1 = 99 # ループ変数の最大値
dn = 1    # ループ変数の増加間隔

#-------------------------------------------------------------------------------
# ループの開始
#-------------------------------------------------------------------------------
load "map_animation.plt" 