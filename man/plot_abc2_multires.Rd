% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_23.r
\name{plot_abc2_multires}
\alias{plot_abc2_multires}
\title{複数の2系ABC計算結果を同時プロットするための関数}
\usage{
plot_abc2_multires(
  res.list,
  stock.name = NULL,
  fishseason = 0,
  detABC = 0,
  abc4 = FALSE,
  cpueunit = "",
  fillarea = FALSE,
  RP = TRUE,
  leftalign = FALSE,
  proposal = TRUE,
  hcrdist = FALSE,
  BThcr = FALSE
)
}
\arguments{
\item{res.list}{calc_abc2の返り値。比較結果が多すぎるとみづらくなるので一度に比較するのは５つまで。}

\item{fishseason}{X軸のラベルを変更（0なら年、1なら漁期年）}

\item{abc4}{北海道東部の跨り資源で使用する図を描画（TRUEなら使用、デフォルトはFALSE））}

\item{cpueunit}{資源量指標値の縦軸見出しに追記したい指標値の単位（例えば"（トン/網）"のように指定する）}

\item{RP}{資源量指標値/年のプロットでReference Point（目標・限界管理基準線）を載せる・載せない（デフォルトはTRUE、FALSEでは直近年の資源量指標値をポイントでハイライトする）}

\item{leftalign}{資源量指標値の時系列の長さが漁獲量に比べて短い時、データが無い範囲の空間を削除する（TRUEなら使用、デフォルトはFALSE）}
}
\description{
複数の2系ABC計算結果を同時プロットするための関数
}
