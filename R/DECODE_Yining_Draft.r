H <- 1
T <- 1

p1 <- c()
p2 <- c()
for (i in 1:10000) {
    pis_prior_unnorm <- runif(H+T, 0, 1)
    pis_prior <- pis_prior_unnorm / sum(pis_prior_unnorm)
    p1 <- c(p1, pis_prior[1])   
    p2 <- c(p2, pis_prior[2])
}
par(mfrow = c(1, 2))
hist(p1, main = "Histogram of P1", xlab = "P1")
hist(p2, main = "Histogram of P2", xlab = "P2")

###########################
# 设置参数
library(MCMCpack)
H <- 1
T <- 1
alpha <- rep(1, H + T)  # Dirichlet 分布的参数，所有元素为 1 表示均匀分布

# 初始化存储结果的向量
p1 <- c()
p2 <- c()

# 生成 10000 个和为 1 的随机向量
for (i in 1:10000) {
    pis_prior <- rdirichlet(1, alpha)
    p1 <- c(p1, pis_prior[1])
    p2 <- c(p2, pis_prior[2])
}

# 绘制直方图
par(mfrow = c(1, 2))
hist(p1, main = "Histogram of P1", xlab = "P1")
hist(p2, main = "Histogram of P2", xlab = "P2")