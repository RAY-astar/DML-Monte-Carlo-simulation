# 双重机器学习(DML)处理效应估计的蒙特卡洛模拟BOOST

# 清空环境
rm(list=ls())

# 加载必要的包
if (!requireNamespace("kableExtra", quietly = TRUE)) {
  install.packages("kableExtra")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}

library(caret)
library(gbm)
library(sandwich)
library(lmtest)
library(reshape2)
library(ggplot2)
library(knitr)
library(kableExtra)
library(dplyr)
library(scales)
library(openxlsx)
##############################################
# 第一步：定义五种不同的数据生成情形
##############################################

# 数据生成函数
generate_data <- function(scenario, n, theta0=1) {
  # 根据不同情景生成数据
  if (scenario == 1) {
    # 情形1: 简单线性模拟
    X1 <- runif(n, min=0, max=1)
    X2 <- runif(n, min=0, max=1)
    
    m0 <- (X1 + X2) / 2
    f0 <- -6 + 12 * X1 + 6 * X2
    
    X <- cbind(X1, X2)
    colnames(X) <- c("X1", "X2")
  } else if (scenario == 2) {
    # 情形2: 线性模拟添加交互项
    X1 <- rnorm(n, mean=0, sd=1)
    X2 <- rnorm(n, mean=0, sd=1)
    X3 <- X1 * X2  # 交互项
    
    D_temp <- X1 + X2 + X3
    m0 <- 1 / (1 + exp(-D_temp))
    f0 <- -6 + 12 * X1 + 6 * X2 + 2 * X3
    
    X <- cbind(X1, X2, X3)
    colnames(X) <- c("X1", "X2", "X3")
  } else if (scenario == 3) {
    # 情形3: 线性模拟添加平方项
    X1 <- rnorm(n, mean=0, sd=1)
    X2 <- rnorm(n, mean=0, sd=1)
    X3 <- X1 * X1  # 平方项
    X4 <- X2 * X2  # 平方项
    
    D_temp <- X1 + X2 + X3 + X4
    m0 <- 1 / (1 + exp(-D_temp))
    f0 <- -6 + 12 * X1 + 6 * X2 + 2 * X3 + 3 * X4
    
    X <- cbind(X1, X2, X3, X4)
    colnames(X) <- c("X1", "X2", "X3", "X4")
  } else if (scenario == 4) {
    # 情形4: 线性模拟添加交互项和平方项
    X1 <- rnorm(n, mean=0, sd=1)
    X2 <- rnorm(n, mean=0, sd=1)
    X3 <- X1 * X1  # 平方项
    X4 <- X2 * X2  # 平方项
    X5 <- X1 * X2  # 交互项
    
    D_temp <- X1 + X2 + X3 + X4 + X5
    m0 <- 1 / (1 + exp(-D_temp))
    f0 <- -6 + 12 * X1 + 6 * X2 + 2 * X3 + 3 * X4 + 3 * X5
    
    X <- cbind(X1, X2, X3, X4, X5)
    colnames(X) <- c("X1", "X2", "X3", "X4", "X5")
  } else if (scenario == 5) {
    # 情形5: 指数型非线性 - 优化处理以提高数值稳定性
    X1 <- rnorm(n, mean=0, sd=1)
    X2 <- rnorm(n, mean=0, sd=1)
    X3 <- rnorm(n, mean=0, sd=1)
    X4 <- rnorm(n, mean=0, sd=1)
    X5 <- rnorm(n, mean=0, sd=1)
    X_matrix <- cbind(X1, X2, X3, X4, X5)
    
    # 计算f0(Xi) - 添加数值稳定性保护措施
    f0 <- rep(1, n)
    denominator <- sqrt(max(0.001, -0.5 * exp(2) + 2 * exp(1) - 1.5))
    
    for (j in 1:5) {
      # 限制指数增长以防止数值溢出
      exp_val <- pmin(exp(X_matrix[,j]), 1e10)  
      f0 <- f0 + (4 * (exp_val - exp(1) + 1)) / denominator
    }
    
    # 计算m0(Xi) - 使用数值稳定的方法
    m0 <- rep(0.5, n)
    log2 <- log(2)  # 预计算log(2)提高效率
    
    for (j in 1:5) {
      # 防止log(负数)或log(0)
      x_safe <- pmax(1e-10, X_matrix[,j] + 1)
      log_term <- log(x_safe) / log2
      m0 <- m0 + (-0.1 + 0.2 * (log2 - log_term))
    }
    
    # 确保m0严格在(0,1)范围内，避免极端值
    m0 <- pmax(0.01, pmin(0.99, m0))
    
    X <- X_matrix
    colnames(X) <- c("X1", "X2", "X3", "X4", "X5")
  }
  
  # 生成误差项和处理变量
  U <- rnorm(n, mean=0, sd=1)
  D <- rbinom(n, size=1, prob=m0)
  
  # 生成结果变量Y
  Y <- D * theta0 + f0 + U
  
  # 返回数据
  data <- data.frame(Y=Y, D=D, X)
  return(list(data=data, y=Y, d=D, X=X))
}

##############################################
# 第二步：实现优化的双重机器学习 - 提升法(Boosting)
##############################################

dml_boost <- function(data, y, d, X, seed=123) {
  # 设置随机种子
  set.seed(seed)
  
  # 获取样本量和特征数
  n <- nrow(data)
  p <- ncol(X)
  
  # 计算适当的模型参数
  n_trees <- max(500, min(5000, 200 * p))  # 根据特征数量动态设置树的数量
  depth <- min(5, max(3, floor(sqrt(p))))  # 根据特征数量动态设置树的深度
  shrinkage_val <- 0.01                    # 较小的学习率提高稳定性
  
  # 使用分层抽样将样本分成两部分，确保处理组和对照组的比例保持一致
  trainIndex <- createDataPartition(d, p=0.5, list=FALSE)
  data1 <- data[trainIndex, ]
  data2 <- data[-trainIndex, ]
  y1 <- y[trainIndex]
  y2 <- y[-trainIndex]
  d1 <- d[trainIndex]
  d2 <- d[-trainIndex]
  X1 <- X[trainIndex, ]
  X2 <- X[-trainIndex, ]
  
  # 定义函数用于拟合提升模型并返回最佳迭代次数
  fit_boosting <- function(y_var, x_data, distribution="gaussian") {
    model <- gbm(y_var ~ ., 
                 data=data.frame(x_data), 
                 distribution=distribution, 
                 n.trees=n_trees, 
                 interaction.depth=depth, 
                 shrinkage=shrinkage_val,
                 bag.fraction=0.8,
                 train.fraction=0.8,
                 verbose=FALSE)
    
    # 获取最佳迭代次数
    best_iter <- gbm.perf(model, method="OOB", plot.it=FALSE)
    return(list(model=model, best_iter=best_iter))
  }
  
  # 第一部分数据上拟合Y~X和D~X
  y1_model <- fit_boosting(y1, X1)
  d1_model <- fit_boosting(d1, X1, distribution="bernoulli")
  
  # 预测第二部分数据
  yhat1 <- predict(y1_model$model, newdata=data.frame(X2), n.trees=y1_model$best_iter)
  dhat1 <- predict(d1_model$model, newdata=data.frame(X2), n.trees=d1_model$best_iter, type="response")
  dhat1 <- pmax(0.01, pmin(0.99, dhat1))  # 确保预测值在合理范围
  
  # 计算残差
  res_y1 <- y2 - yhat1
  res_d1 <- d2 - dhat1
  
  # 第二部分数据上拟合Y~X和D~X
  y2_model <- fit_boosting(y2, X2)
  d2_model <- fit_boosting(d2, X2, distribution="bernoulli")
  
  # 预测第一部分数据
  yhat2 <- predict(y2_model$model, newdata=data.frame(X1), n.trees=y2_model$best_iter)
  dhat2 <- predict(d2_model$model, newdata=data.frame(X1), n.trees=d2_model$best_iter, type="response")
  dhat2 <- pmax(0.01, pmin(0.99, dhat2))  # 确保预测值在合理范围
  
  # 计算残差
  res_y2 <- y1 - yhat2
  res_d2 <- d1 - dhat2
  
  # 检测和处理异常值
  is_outlier <- function(x) {
    if(length(unique(x)) <= 3) return(rep(FALSE, length(x)))  # 如果取值太少，不做检测
    
    quantiles <- quantile(x, probs=c(0.25, 0.75), na.rm=TRUE)
    IQR <- quantiles[2] - quantiles[1]
    if(IQR < 1e-10) return(rep(FALSE, length(x)))  # 如果IQR过小，不做检测
    
    lower_bound <- quantiles[1] - 3 * IQR
    upper_bound <- quantiles[2] + 3 * IQR
    return(x < lower_bound | x > upper_bound)
  }
  
  # 对两部分数据分别进行异常值检测
  outliers1 <- is_outlier(res_y1) | is_outlier(res_d1)
  outliers2 <- is_outlier(res_y2) | is_outlier(res_d2)
  
  # 如果过滤后样本太少，则不过滤
  if(sum(!outliers1) < length(res_y1) * 0.9) outliers1 <- rep(FALSE, length(res_y1))
  if(sum(!outliers2) < length(res_y2) * 0.9) outliers2 <- rep(FALSE, length(res_y2))
  
  # 计算有效样本量
  n1_eff <- sum(!outliers1)
  n2_eff <- sum(!outliers2)
  
  # 回归残差获取处理效应估计
  tryCatch({
    lm_fit1 <- lm(res_y1[!outliers1] ~ res_d1[!outliers1])
    lm_fit2 <- lm(res_y2[!outliers2] ~ res_d2[!outliers2])
    
    # 计算稳健标准误
    est1 <- coeftest(lm_fit1, vcov=vcovHC, type="HC0")
    est2 <- coeftest(lm_fit2, vcov=vcovHC, type="HC0")
    
    # 方法1: 分别回归并求加权平均
    w1 <- n1_eff / (n1_eff + n2_eff)
    w2 <- n2_eff / (n1_eff + n2_eff)
    
    b1 <- est1[2, 1]
    b2 <- est2[2, 1]
    be1 <- w1 * b1 + w2 * b2
    
    se1 <- est1[2, 2]
    se2 <- est2[2, 2]
    sig2 <- w1 * (se1^2 + (b1 - be1)^2) + w2 * (se2^2 + (b2 - be1)^2)
    se_avg <- sqrt(sig2)
    
    # 方法2: 合并残差进行回归
    res_y_combined <- c(res_y1[!outliers1], res_y2[!outliers2])
    res_d_combined <- c(res_d1[!outliers1], res_d2[!outliers2])
    
    lm_fit_combined <- lm(res_y_combined ~ res_d_combined)
    est_combined <- coeftest(lm_fit_combined, vcov=vcovHC, type="HC0")
    
    be2 <- est_combined[2, 1]
    se_combined <- est_combined[2, 2]
    
    # 检查估计结果是否合理
    if(abs(be1) > 10 || abs(be2) > 10) {
      be1 <- NA
      be2 <- NA
      se_avg <- NA
      se_combined <- NA
    }
  }, error=function(e) {
    be1 <- NA
    be2 <- NA
    se_avg <- NA
    se_combined <- NA
  })
  
  # 返回估计结果
  result <- list(
    estimate1 = be1,
    se1 = se_avg,
    estimate2 = be2,
    se2 = se_combined
  )
  
  return(result)
}

##############################################
# 第三步：蒙特卡洛模拟主函数
##############################################

monte_carlo_dml_boost <- function(scenario, n, n_rep=1000, seed=123) {
  # 初始化结果存储矩阵
  results <- matrix(NA, nrow=n_rep, ncol=4)
  colnames(results) <- c("Boost_est1", "Boost_se1", "Boost_est2", "Boost_se2")
  
  # 设置全局随机种子
  set.seed(seed)
  
  # 重复实验n_rep次
  successful_runs <- 0
  
  for (i in 1:n_rep) {
    if (i %% 50 == 0) {
      cat(sprintf("情形 %d, 样本量 %d, 迭代 %d/%d\n", scenario, n, i, n_rep))
    }
    
    # 生成数据
    current_seed <- seed + i
    tryCatch({
      # 生成模拟数据
      sim_data <- generate_data(scenario=scenario, n=n, theta0=1)
      
      # 使用提升法进行DML估计
      result_boost <- dml_boost(
        data=sim_data$data, 
        y=sim_data$y, 
        d=sim_data$d, 
        X=sim_data$X, 
        seed=current_seed
      )
      
      # 存储结果
      results[i, ] <- c(result_boost$estimate1, result_boost$se1, 
                        result_boost$estimate2, result_boost$se2)
      
      # 检查是否为有效结果
      if(!is.na(result_boost$estimate1) && !is.na(result_boost$estimate2)) {
        successful_runs <- successful_runs + 1
      }
    }, error=function(e) {
      results[i, ] <- NA
    })
  }
  
  # 输出成功率
  success_rate <- successful_runs / n_rep * 100
  cat(sprintf("情形 %d, 样本量 %d 的成功率: %.1f%%\n", scenario, n, success_rate))
  
  return(list(results=results, success_rate=success_rate))
}

# 计算汇总统计量
summarize_results <- function(results, true_value=1) {
  # 删除NA行
  valid_rows <- complete.cases(results)
  valid_results <- results[valid_rows, ]
  
  if(nrow(valid_results) == 0) {
    return(data.frame(
      Method = c("Boosting-Est1", "Boosting-Est2"),
      Mean = c(NA, NA),
      SD = c(NA, NA),
      Bias = c(NA, NA),
      RMSE = c(NA, NA),
      Success_Rate = 0
    ))
  }
  
  # 计算均值
  boost_mean1 <- mean(valid_results[, "Boost_est1"])
  boost_mean2 <- mean(valid_results[, "Boost_est2"])
  
  # 计算标准差
  boost_sd1 <- sd(valid_results[, "Boost_est1"])
  boost_sd2 <- sd(valid_results[, "Boost_est2"])
  
  # 计算偏差(Bias)
  boost_bias1 <- boost_mean1 - true_value
  boost_bias2 <- boost_mean2 - true_value
  
  # 计算均方误差(MSE)和均方根误差(RMSE)
  boost_mse1 <- mean((valid_results[, "Boost_est1"] - true_value)^2)
  boost_mse2 <- mean((valid_results[, "Boost_est2"] - true_value)^2)
  boost_rmse1 <- sqrt(boost_mse1)
  boost_rmse2 <- sqrt(boost_mse2)
  
  # 成功率
  success_rate <- sum(valid_rows) / nrow(results)
  
  # 组织结果
  summary_table <- data.frame(
    Method = c("Boosting-Est1", "Boosting-Est2"),
    Mean = c(boost_mean1, boost_mean2),
    SD = c(boost_sd1, boost_sd2),
    Bias = c(boost_bias1, boost_bias2),
    RMSE = c(boost_rmse1, boost_rmse2),
    Success_Rate = success_rate * 100  # 将成功率转换为百分比
  )
  
  return(summary_table)
}

##############################################
# 执行蒙特卡洛模拟
##############################################

# 定义样本容量和模拟情景
sample_sizes <- c(100, 200,500,1000,2000,5000,10000)
scenarios <- 1:5  # 五种不同的数据生成过程

# 初始化最终结果存储
all_results <- list()

# 设置模拟参数 (实际应用中可增大n_rep)
n_rep_demo <- 50

# 执行所有情形和样本量的模拟
for (scenario in scenarios) {
  scenario_results <- list()
  
  cat("\n========== 开始情形", scenario, "的模拟 ==========\n")
  
  for (n in sample_sizes) {
    cat("\n开始样本量", n, "的模拟...\n")
    
    # 执行蒙特卡洛模拟
    start_time <- Sys.time()
    mc_results <- monte_carlo_dml_boost(scenario=scenario, n=n, n_rep=n_rep_demo, seed=123)
    end_time <- Sys.time()
    
    # 计算汇总统计量
    summary <- summarize_results(mc_results$results)
    scenario_results[[as.character(n)]] <- summary
    
    # 输出当前情形和样本量的结果
    cat("用时:", round(difftime(end_time, start_time, units="mins"), 2), "分钟\n")
    print(summary)
  }
  
  all_results[[paste0("scenario_", scenario)]] <- scenario_results
}

##############################################
# 结果展示
##############################################

# 生成表格展示结果
create_summary_table <- function(all_results, scenario_id) {
  # 获取特定情形的结果
  scenario_results <- all_results[[paste0("scenario_", scenario_id)]]
  
  # 合并所有样本量的结果
  result_table <- data.frame()
  for (n in names(scenario_results)) {
    temp_df <- scenario_results[[n]]
    temp_df$SampleSize <- as.numeric(n)
    result_table <- rbind(result_table, temp_df)
  }
  
  # 整理表格格式
  result_table <- result_table[, c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate")]
  
  # 按样本量和方法排序
  result_table <- result_table[order(result_table$SampleSize, result_table$Method), ]
  
  return(result_table)
}

# 为每个情形创建结果表格并显示
for(i in 1:5) {
  table_scenario <- create_summary_table(all_results, i)
  print(kable(table_scenario, digits=4, caption=paste0("情形", i, "的结果")))
}

##############################################
# 可视化结果
##############################################

# 将所有情形的结果合并为一个数据框
plot_data <- data.frame()

for (scenario in 1:5) {
  scenario_name <- paste0("情形", scenario)
  
  for (n in names(all_results[[paste0("scenario_", scenario)]])) {
    temp_df <- all_results[[paste0("scenario_", scenario)]][[n]]
    temp_df$Scenario <- scenario_name
    temp_df$SampleSize <- as.numeric(n)
    plot_data <- rbind(plot_data, temp_df)
  }
}

# 绘制均值随样本量变化的图表
mean_plot <- ggplot(plot_data, aes(x=SampleSize, y=Mean, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  geom_hline(yintercept=1, linetype="dashed", color="black") +  # 添加真实值参考线
  labs(title="提升法(Boosting)在不同情形和样本量下的处理效应估计",
       x="样本容量", y="估计值均值") +
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 绘制RMSE随样本量变化的图表
rmse_plot <- ggplot(plot_data, aes(x=SampleSize, y=RMSE, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="提升法(Boosting)在不同情形和样本量下的均方根误差(RMSE)",
       x="样本容量", y="RMSE") +
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 绘制偏差随样本量变化的图表
bias_plot <- ggplot(plot_data, aes(x=SampleSize, y=Bias, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="提升法(Boosting)在不同情形和样本量下的偏差(Bias)",
       x="样本容量", y="Bias") +
  geom_hline(yintercept=0, linetype="dashed", color="black") +  # 添加零偏差参考线
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 绘制成功率随样本量变化的图表
success_rate_plot <- ggplot(plot_data, aes(x=SampleSize, y=Success_Rate, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="提升法(Boosting)在不同情形和样本量下的成功率",
       x="样本容量", y="成功率 (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # 将成功率显示为百分比
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 创建成功率热图
# 首先确保数据格式正确
heatmap_data <- plot_data %>%
  select(Method, Scenario, SampleSize, Success_Rate) %>%
  distinct()  # 确保没有重复

success_heatmap <- ggplot(heatmap_data, aes(x=SampleSize, y=Method, fill=Success_Rate)) +
  geom_tile() +
  scale_fill_gradient(low="red", high="green", 
                      labels=scales::percent_format(scale=1),
                      limits=c(0, 100)) +
  labs(title="提升法(Boosting)成功率热图",
       x="样本容量", y="方法", fill="成功率 (%)") +
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 显示图表
print(mean_plot)
print(rmse_plot)
print(bias_plot)
print(success_rate_plot)
print(success_heatmap)

# 将结果保存到文件
save(all_results, plot_data, file="dml_boosting_results.RData")

# 可选：保存图表为图片文件
ggsave("dml_boost_mean.png", mean_plot, width=10, height=8)
ggsave("dml_boost_rmse.png", rmse_plot, width=10, height=8)
ggsave("dml_boost_bias.png", bias_plot, width=10, height=8)
ggsave("dml_boost_success_rate.png", success_rate_plot, width=10, height=8)
ggsave("dml_boost_success_heatmap.png", success_heatmap, width=12, height=6)
# 合并所有结果并导出表格

# 创建包含所有信息的总表
create_complete_summary_table <- function(all_results) {
  # 初始化空数据框
  complete_table <- data.frame()
  
  # 遍历所有情形
  for (scenario in 1:5) {
    scenario_name <- paste0("Scenario_", scenario)
    scenario_results <- all_results[[paste0("scenario_", scenario)]]
    
    # 遍历每个样本量
    for (n in names(scenario_results)) {
      temp_df <- scenario_results[[n]]
      temp_df$Scenario <- scenario_name
      temp_df$SampleSize <- as.numeric(n)
      complete_table <- rbind(complete_table, temp_df)
    }
  }
  
  # 重新排列列的顺序
  complete_table <- complete_table[, c("Scenario", "SampleSize", "Method", "Mean", "SD", "Bias", "RMSE", "Success_Rate")]
  
  # 按情形、样本量和方法排序
  complete_table <- complete_table[order(complete_table$Scenario, complete_table$SampleSize, complete_table$Method), ]
  
  return(complete_table)
}

# 创建更宽的汇总表格（以样本量为列）
create_wide_summary_table <- function(all_results) {
  wide_table <- data.frame()
  
  for (scenario in 1:5) {
    scenario_name <- paste0("Scenario_", scenario)
    scenario_results <- all_results[[paste0("scenario_", scenario)]]
    
    # 为每个方法创建一行
    for (method in c("Boosting-Est1", "Boosting-Est2")) {
      row_data <- data.frame(Scenario = scenario_name, Method = method)
      
      # 为每个样本量添加列
      for (n in names(scenario_results)) {
        temp_df <- scenario_results[[n]]
        method_data <- temp_df[temp_df$Method == method, ]
        
        # 创建带有样本量标签的列名
        n_label <- paste0("n", n)
        row_data[[paste0(n_label, "_Mean")]] <- method_data$Mean
        row_data[[paste0(n_label, "_SD")]] <- method_data$SD
        row_data[[paste0(n_label, "_Bias")]] <- method_data$Bias
        row_data[[paste0(n_label, "_RMSE")]] <- method_data$RMSE
        row_data[[paste0(n_label, "_Success")]] <- method_data$Success_Rate
      }
      
      wide_table <- rbind(wide_table, row_data)
    }
  }
  
  return(wide_table)
}

# 创建总表
complete_summary <- create_complete_summary_table(all_results)
wide_summary <- create_wide_summary_table(all_results)

# 查看总表
print("完整汇总表格（长格式）:")
print(kable(complete_summary, digits = 4, caption = "所有情形和样本量的完整结果"))

# 导出为CSV文件
write.csv(complete_summary, "dml_boosting_complete_results.csv", row.names = FALSE)
write.csv(wide_summary, "dml_boosting_wide_results.csv", row.names = FALSE)

# 导出为Excel文件（需要openxlsx包）
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  
  # 创建工作簿
  wb <- createWorkbook()
  
  # 添加长格式表格
  addWorksheet(wb, "Complete_Results")
  writeData(wb, "Complete_Results", complete_summary)
  
  # 添加宽格式表格
  addWorksheet(wb, "Wide_Results")
  writeData(wb, "Wide_Results", wide_summary)
  
  # 为每个情形创建单独的工作表
  for (scenario in 1:5) {
    scenario_data <- complete_summary[complete_summary$Scenario == paste0("Scenario_", scenario), ]
    addWorksheet(wb, paste0("Scenario_", scenario))
    writeData(wb, paste0("Scenario_", scenario), scenario_data)
  }
  
  # 创建汇总统计工作表
  summary_stats <- data.frame(
    Scenario = unique(complete_summary$Scenario),
    Avg_Mean = sapply(unique(complete_summary$Scenario), function(s) {
      mean(complete_summary[complete_summary$Scenario == s, "Mean"], na.rm = TRUE)
    }),
    Avg_Bias = sapply(unique(complete_summary$Scenario), function(s) {
      mean(complete_summary[complete_summary$Scenario == s, "Bias"], na.rm = TRUE)
    }),
    Avg_RMSE = sapply(unique(complete_summary$Scenario), function(s) {
      mean(complete_summary[complete_summary$Scenario == s, "RMSE"], na.rm = TRUE)
    }),
    Avg_Success_Rate = sapply(unique(complete_summary$Scenario), function(s) {
      mean(complete_summary[complete_summary$Scenario == s, "Success_Rate"], na.rm = TRUE)
    })
  )
  
  addWorksheet(wb, "Summary_Statistics")
  writeData(wb, "Summary_Statistics", summary_stats)
  
  # 保存Excel文件
  saveWorkbook(wb, "dml_boosting_results_all.xlsx", overwrite = TRUE)
  
  print("Excel文件已保存为: dml_boosting_results_all.xlsx")
} else {
  print("请安装openxlsx包以导出Excel文件：install.packages('openxlsx')")
}

# 创建更漂亮的HTML表格（可选）
if (requireNamespace("kableExtra", quietly = TRUE)) {
  library(kableExtra)
  
  # 创建样式化的表格
  html_table <- complete_summary %>%
    kable(format = "html", digits = 4, 
          caption = "双重机器学习提升法(Boosting)完整结果") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                  full_width = FALSE) %>%
    scroll_box(width = "100%", height = "500px")
  
  # 保存HTML表格
  save_kable(html_table, file = "dml_boosting_results_table.html")
  
  print("HTML表格已保存为: dml_boosting_results_table.html")
}

# 创建汇总图表
library(ggplot2)
library(tidyr)

# 创建汇总性能图表
performance_plot <- ggplot(complete_summary, 
                           aes(x = SampleSize, y = RMSE, color = Method, shape = Scenario)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "提升法性能比较 - 所有情形",
       x = "样本容量",
       y = "均方根误差 (RMSE)") +
  scale_x_continuous(breaks = unique(complete_summary$SampleSize)) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 保存图表
ggsave("dml_boost_performance_summary.png", performance_plot, width = 12, height = 8)

# 创建成功率汇总图表
success_plot <- ggplot(complete_summary, 
                       aes(x = SampleSize, y = Success_Rate, fill = Scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Method) +
  labs(title = "提升法成功率比较 - 所有情形",
       x = "样本容量",
       y = "成功率 (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 保存图表
ggsave("dml_boost_success_summary.png", success_plot, width = 12, height = 8)

# 创建估计均值汇总图
mean_summary_plot <- ggplot(complete_summary, 
                            aes(x = SampleSize, y = Mean, color = Method, shape = Scenario)) +
  geom_line() +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # 添加真实值参考线
  labs(title = "提升法估计均值比较 - 所有情形",
       x = "样本容量",
       y = "估计均值") +
  scale_x_continuous(breaks = unique(complete_summary$SampleSize)) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 保存图表
ggsave("dml_boost_mean_summary.png", mean_summary_plot, width = 12, height = 8)

# 创建估计偏差汇总图
bias_summary_plot <- ggplot(complete_summary, 
                            aes(x = SampleSize, y = Bias, color = Method, shape = Scenario)) +
  geom_line() +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 添加零偏差参考线
  labs(title = "提升法估计偏差比较 - 所有情形",
       x = "样本容量",
       y = "估计偏差") +
  scale_x_continuous(breaks = unique(complete_summary$SampleSize)) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 保存图表
ggsave("dml_boost_bias_summary.png", bias_summary_plot, width = 12, height = 8)

# 更新打印文件保存列表
cat("\n已创建以下文件:\n")
cat("1. dml_boosting_complete_results.csv - 长格式完整结果\n")
cat("2. dml_boosting_wide_results.csv - 宽格式完整结果\n")
cat("3. dml_boosting_results_all.xlsx - Excel工作簿，包含多个工作表\n")
cat("4. dml_boosting_results_table.html - 交互式HTML表格\n")
cat("5. dml_boost_mean.png - 均值图表\n")
cat("6. dml_boost_rmse.png - 均方根误差图表\n")
cat("7. dml_boost_bias.png - 偏差图表\n")
cat("8. dml_boost_success_rate.png - 成功率图表\n")
cat("9. dml_boost_success_heatmap.png - 成功率热图\n")
cat("10. dml_boost_performance_summary.png - RMSE性能汇总图\n")
cat("11. dml_boost_mean_summary.png - 估计均值汇总图\n")
cat("12. dml_boost_bias_summary.png - 估计偏差汇总图\n")
cat("13. dml_boost_success_summary.png - 成功率汇总图\n")


# 返回表格数据
list(
  complete_table = complete_summary,
  wide_table = wide_summary,
  summary_stats = summary_stats
)