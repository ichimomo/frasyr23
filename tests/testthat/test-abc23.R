library(frasyr23)

context("type23 check")

test_that("type23 check",{
    # example dataを使ったテスト
    catch <- c(10,5,10,3,2,1,3)
    cpue <- c(11,10,2,3,2,5,2)
    example_data <- data.frame(year=1:7,cpue=cpue,catch=catch)

    # 2系
    example_abc2_check <- calc_abc2(example_data)
    #graph_abc2 <- plot_abc2(example_abc2_check)

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
    #graph_abc3 <- plot_abc3(example_abc3_check)

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
    # 上記結果の読み込み
    load(system.file("extdata","res_aka_abc3.rda",package = "frasyr23"))
    # テスト内容
    testcontents<-c("BRP","Obs_BRP","Current_Status","ABC")
    # 結果の照合
    for(i in 1:length(testcontents)){
        expect_equal(eval(parse(text=paste("aka_abc3_check$",testcontents[i]))), eval(parse(text=paste("aka_abc3$",testcontents[i]))),label=c(testcontents[i]))
    }

})


