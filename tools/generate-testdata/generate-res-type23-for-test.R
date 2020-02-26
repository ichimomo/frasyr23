source("./tools/generate-testdata/abc_t23.r")

#data_akaをつかう。
data("data_aka")
#head(data_aka)

# 2系
aka_abc2 <- abc_t23(catch=data_aka$catch,cpue = data_aka$cpue)
save(aka_abc2,file="./inst/extdata/res_aka_abc2.rda")

# 3系
aka_abc3 <- abc_t23(catch=data_aka$catch)
save(aka_abc3,file="./inst/extdata/res_aka_abc3.rda")

# exampleを使う
catch <- c(10,5,10,3,2,1,3)
cpue <- c(11,10,2,3,2,5,2)
example_data <- data.frame(year=1:7,cpue=cpue,catch=catch)

# 2系
example_abc2 <- abc_t23(catch=example_data$catch,cpue=example_data$cpue)
save(example_abc2,file="./inst/extdata/res_example_abc2.rda")
# 3系
example_abc3 <- abc_t23(catch=example_data$catch)
save(example_abc3,file="./inst/extdata/res_example_abc3.rda")
