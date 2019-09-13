library(frasyr23)

context("type23 check")

test_that("type23 check",{
    catch <- c(10,5,10,3,2,1,3)
    cpue <- c(11,10,2,3,2,5,2)
    example_data <- data.frame(year=1:7,cpue=cpue,catch=catch)

    # 2系
    example_abc2 <- calc_abc2(example_data)
    graph_abc2 <- plot_abc2(example_abc2)
    expect_equal(round(example_abc2$ABC,6), 1.894191)

    # 3系
    example_abc3 <- calc_abc3(example_data,BT=0.1,PL=2,PB=10,tune.par=c(1.5,2))
    graph_abc3 <- plot_abc3(example_abc3)
    expect_equal(round(example_abc3$ABC,6), 2.658756)

})


