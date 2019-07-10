# frasyr23
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R for data poor fisheries stocks (type 2 and type3)
- ２系・３系用のABC計算パッケージ（試験運用中）です

# インストールと呼び出し
```
# install.pakcages("devtools") # <-- devtoolsをインストールしていない人
devtools::install_github("ichimomo/frasyr23") # frasyrのインストール
library(frasyr23) # frasyrの呼び出し
```

# 主な関数と使い方
- 2系の計算
   - calc_abc2 ABCの計算
   - plot_abc2 結果のプロット

- 3系の計算
   - calc_abc3 ABCの計算
   - plot_abc3 結果のプロット

```
help(calc_abc2) # helpを見ると引数の説明などが見れます

# データの作成
catch <- c(10,5,10,3,2,1,3)
cpue <- c(11,10,2,3,2,5,2)
example_data <- data.frame(year=1:7,cpue=cpue,catch=catch)

# 2系
example_abc2 <- calc_abc2(example_data)
example_abc2$ABC
# $ABC=1.894191
graph_abc2 <- plot_abc2(example_abc2)

# 3系
## dataにCPUEが入っていても無視します
example_abc3 <- calc_abc3(example_data)
example_abc3$ABC
# [1] 2.658756
graph_abc3 <- plot_abc3(example_abc3)

```
