library(frasyr23)

context("type23 check")

test_that("type23 check",{
    # example dataを使ったテスト
    catch <- c(10,5,10,3,2,1,3)
    cpue <- c(11,10,2,3,2,5,2)
    example_data <- data.frame(year=1:7,cpue=cpue,catch=catch)

    # 2系
    example_abc2_check <- calc_abc2(example_data)
    png("graph.png")
    graph_abc2 <- plot_abc2(example_abc2_check,fishseason = 0,detABC = 1)
    dev.off()

    # 上記結果の読み込み
    load(system.file("extdata","res_example_abc2.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","AAV","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("example_abc2_check$",testcontents[i]))), eval(parse(text=paste("example_abc2$",testcontents[i]))),label=c(testcontents[i]))
    }

    # 3系
    example_abc3_check <- calc_abc3(example_data)
    png("graph.png")
    graph_abc3 <- plot_abc3(example_abc3_check)
    dev.off()

    # 上記結果の読み込み
    load(system.file("extdata","res_example_abc3.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("example_abc3_check$",testcontents[i]))), eval(parse(text=paste("example_abc3$",testcontents[i]))),label=c(testcontents[i]))
    }


    # data_akaを使ったテスト
    data("data_aka")
    # 2系
    aka_abc2_check <- calc_abc2(data_aka)
    #graph_aka_abc2 <- plot_abc2(aka_abc2_check)
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc2.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","AAV","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc2_check$",testcontents[i]))), eval(parse(text=paste("aka_abc2$",testcontents[i]))),label=c(testcontents[i]))
    }

    # 3系
    aka_abc3_check <- calc_abc3(data_aka)
    #graph_aka_abc3 <- plot_abc3(aka_abc3_check)
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc3.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc3_check$",testcontents[i]))), eval(parse(text=paste("aka_abc3$",testcontents[i]))),label=c(testcontents[i]))
    }

})


context("type2_func")

test_that("alpha check",{
    expect_equal(0, type2_func(0, 1:10))
    expect_equal(1, type2_func(0.8, 1:10))
    expect_equal(1.105171/type2_func(1, 1:10),c(1),tol=1.0e-6)
})

context("type23 option check")

test_that("type23 option check",{

    # 経験分布(ecdf関数)利用のテスト
    # data_akaを使ったテスト
    aka_abc2_empir_check <- calc_abc2(data_aka,empir.dist = T,summary_abc = F)
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc2_empir.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","AAV","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc2_empir_check$",testcontents[i]))), eval(parse(text=paste("aka_abc2_empir$",testcontents[i]))),label=c(testcontents[i]))
    }
    # 経験分布(シンプルな経験分布[min(cpue)=0,max(cpue=1)]関数)利用のテスト
    aka_abc2_empir_simple_check <- calc_abc2(data_aka,empir.dist = T,simple.empir = T, summary_abc = F)
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc2_empir_simple.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","AAV","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc2_empir_simple_check$",testcontents[i]))), eval(parse(text=paste("aka_abc2_empir_simple$",testcontents[i]))),label=c(testcontents[i]))
    }

    # CPUE時系列データの最終年を３年移動平均で平滑化
    aka_abc2_sm3cpue_check <- calc_abc2(data_aka,smooth.cpue = T,n.cpue=3,summary_abc = F)
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc2_sm3cpue.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Recent_Status","AAV","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc2_sm3cpue_check$",testcontents[i]))), eval(parse(text=paste("aka_abc2_sm3cpue$",testcontents[i]))),label=c(testcontents[i]))
    }

    # CPUE時系列データ全体を３年移動平均で平滑化
    aka_abc2_sm3cpuedist_check <- calc_abc2(data_aka,smooth.dist = T,n.cpue=3,summary_abc = F)
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc2_sm3cpuedist.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Recent_Status","AAV","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc2_sm3cpuedist_check$",testcontents[i]))), eval(parse(text=paste("aka_abc2_sm3cpuedist$",testcontents[i]))),label=c(testcontents[i]))
    }


    # 管理基準計算年固定
    aka_abc2_bt2010_check <- calc_abc2(data_aka,BTyear = 2010,summary_abc = F)
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc2_bt2010.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","AAV","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc2_bt2010_check$",testcontents[i]))), eval(parse(text=paste("aka_abc2_bt2010$",testcontents[i]))),label=c(testcontents[i]))
    }
})

