# frasyr23
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R for data poor fisheries stocks (type 2 and type3)
- ２系・３系用のABC計算パッケージ（試験運用中）です

# インストールと呼び出し
```
# install.pakcages("devtools") # <-- devtoolsをインストールしていない人
devtools::install_github("ichimomo/frasyr23") # frasyrのインストール
library(frasyr23) # frasyrの呼び出し
library(tidyverse) # こちらのパッケージを使うので呼び出しておく		  
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

# 例データ
catch <- c(15,20,13,14,11,10,5,10,3,2,1,3)
cpue <- c(10,9,8,4,8,11,10,2,3,2,5,2)
data_example <- data.frame(year=2001:2012,cpue=cpue,catch=catch)

# 2系
abc2_ex <- calc_abc2(data_example)
graph2_ex <- plot_abc2(abc2_ex)
# AAVのちがいを見る	   
abc2_ex_AAV1 <- calc_abc2(data_example,AAV=1)	     

# 3系
## dataにCPUEが入っていても無視します
abc3_ex <- calc_abc3(data_example)
graph3_ex <- plot_abc3(abc3_ex)

# 2系（アカガレイデータ
data(data_aka)
abc2_aka <- calc_abc2(data_aka)
graph2_aka <- plot_abc2(abc2_aka)
# グラフをセーブする場合
# ggsave(graph2_aka[[2]],file="aka2.png")
  	   	
abc3_aka <- calc_abc3(data_aka)
graph3_aka <- plot_abc3(abc3_aka)
# グラフをセーブする場合
# ggsave(graph3_aka[[2]],file="aka3.png")
	  
```
![](tools/aka2.png)	
![](tools/aka3.png)	
