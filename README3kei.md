# frasyr23
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R for data poor fisheries stocks (type 2 and type3)
- ２系(３系)用のABC計算パッケージです
- ３系（計算方法検討中）の関数について説明します

# インストールと呼び出し
```
# install.pakcages("devtools") # <-- devtoolsをインストールしていない人
devtools::install_github("ichimomo/frasyr23") # frasyrのインストール

# 過去の安定版を指定してインストールする場合
# @以下にリリースバージョンを指定します
devtools::install_github("ichimomo/frasyr23@v1.00")

library(frasyr23) # frasyrの呼び出し
library(tidyverse) # こちらのパッケージを使うので呼び出しておく		  
```
- うまくインストールできない場合
- frasyr23と一緒に多くのパッケージが同時にインストールされます．そのパッケージのどれか1つでもうまくインストールできないと，frasry23もインストールできません．対処法としては．．
- 問題があってインストールできないと言われたパッケージを手動でインストールしてみる（install.packages("パッケージ名")
- 古いパッケージが残っていてそれを削除できないためにインストールできない場合もあるみたい．以下のサイトを参考に古いパッケージのファイルを消し，それを手動でインストールしてから再トライ http://www.thirtyfive.info/entry/2017/07/28/R%E3%81%AEplyr%E3%83%91%E3%83%83%E3%82%B1%E3%83%BC%E3%82%B8%E3%81%8C%E8%AA%AD%E3%81%BF%E8%BE%BC%E3%82%81%E3%81%AA%E3%81%84%E5%95%8F%E9%A1%8C%E3%81%AE%E5%AF%BE%E5%87%A6

# 3系の関数
- 3系の計算
- calc_abc3 ABCの計算
- plot_abc3 結果のプロット

# Rコード例
```
# 例データ
catch <- c(15,20,13,14,11,10,5,10,3,2,1,3)
data_example <- data.frame(year=2001:2012,catch=catch)

# 3系
## dataにCPUEが入っていても無視します
abc3_ex <- calc_abc3(data_example)
# ABCが決定できる魚種で、かつ漁期が暦の年に一致する場合
graph3_ex <- plot_abc3(abc3_ex)
# ABCが決定できる魚種で、かつ漁期が暦の年に一致しない場合
graph3_ex <- plot_abc3(abc3_ex,fishseason=1)
# ABCが決定できない魚種で、かつ漁期が暦の年に一致しない場合
graph3_ex <- plot_abc3(abc3_ex,fishseason=1,detABC=1)

```
# 実データの解析例とグラフ

```	  	   	
# アカガレイデータの呼び出し
data(data_aka)
# 3系
abc3_aka <- calc_abc3(data_aka)
graph3_aka <- plot_abc3(abc3_aka)
# グラフをセーブする場合
# ggsave(graph3_aka[[2]],file="aka3.png")
```
![](tools/aka3.png)