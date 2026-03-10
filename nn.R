##############################################
# 双重机器学习(DML)处理效应估计的蒙特卡洛模拟-神经网络(NN)
##############################################

# 清空环境
rm(list=ls())

# 加载必要的包
required_packages <- c("caret", "nnet", "neuralnet", "sandwich", "lmtest", 
                       "reshape2", "ggplot2", "knitr", "dplyr", 
                       "scales", "openxlsx", "tidyr")

# 检查并安装缺失的包
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

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
# 第二步：实现神经网络(NN)的双重机器学习
##############################################

# 双重机器学习 神经网络 函数
dml_nn <- function(data, y, d, X, seed=123, size=1, maxit=300) {
  # 设置随机种子
  set.seed(seed)
  
  # 获取样本量和特征数
  n <- nrow(data)
  p <- ncol(X)
  
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
  
  # 将变量名构建为公式
  x_names <- paste(colnames(X), collapse="+")
  formula_y1 <- as.formula(paste("y1 ~", x_names))
  formula_d1 <- as.formula(paste("d1 ~", x_names))
  formula_y2 <- as.formula(paste("y2 ~", x_names))
  formula_d2 <- as.formula(paste("d2 ~", x_names))
  
  # 创建包含结果变量的数据框
  data1_y <- data.frame(y1=y1, X1)
  data1_d <- data.frame(d1=d1, X1)
  data2_y <- data.frame(y2=y2, X2)
  data2_d <- data.frame(d2=d2, X2)
  
  # 使用tryCatch捕获可能的错误
  tryCatch({
    # 第一部分数据上拟合Y~X
    nn_y1 <- nnet(formula_y1, data=data1_y, size=size, maxit=maxit, linout=TRUE, trace=FALSE)
    
    # 第一部分数据上拟合D~X
    nn_d1 <- nnet(formula_d1, data=data1_d, size=size, maxit=maxit, linout=TRUE, trace=FALSE)
    
    # 预测第二部分数据
    yhat1 <- predict(nn_y1, newdata=as.data.frame(X2))
    dhat1 <- predict(nn_d1, newdata=as.data.frame(X2))
    
    # 计算残差
    res_y1 <- y2 - yhat1
    res_d1 <- d2 - dhat1
    
    # 第二部分数据上拟合Y~X
    nn_y2 <- nnet(formula_y2, data=data2_y, size=size, maxit=maxit, linout=TRUE, trace=FALSE)
    
    # 第二部分数据上拟合D~X
    nn_d2 <- nnet(formula_d2, data=data2_d, size=size, maxit=maxit, linout=TRUE, trace=FALSE)
    
    # 预测第一部分数据
    yhat2 <- predict(nn_y2, newdata=as.data.frame(X1))
    dhat2 <- predict(nn_d2, newdata=as.data.frame(X1))
    
    # 计算残差
    res_y2 <- y1 - yhat2
    res_d2 <- d1 - dhat2
    
    # 检测和处理异常值
    is_outlier <- function(x) {
      quantiles <- quantile(x, probs=c(0.25, 0.75), na.rm=TRUE)
      IQR <- quantiles[2] - quantiles[1]
      if(IQR < 1e-10) return(rep(FALSE, length(x)))
      
      lower_bound <- quantiles[1] - 3 * IQR
      upper_bound <- quantiles[2] + 3 * IQR
      return(x < lower_bound | x > upper_bound)
    }
    
    # 对两部分数据分别进行异常值检测
    outliers1 <- is_outlier(res_y1) | is_outlier(res_d1) | is.na(res_y1) | is.na(res_d1)
    outliers2 <- is_outlier(res_y2) | is_outlier(res_d2) | is.na(res_y2) | is.na(res_d2)
    
    # 如果过滤后样本太少，则不过滤
    if(sum(!outliers1) < length(res_y1) * 0.9) outliers1 <- rep(FALSE, length(res_y1))
    if(sum(!outliers2) < length(res_y2) * 0.9) outliers2 <- rep(FALSE, length(res_y2))
    
    # 计算有效样本量
    n1_eff <- sum(!outliers1)
    n2_eff <- sum(!outliers2)
    
    if(n1_eff == 0 || n2_eff == 0) {
      return(list(
        estimate1 = NA,
        se1 = NA,
        estimate2 = NA,
        se2 = NA
      ))
    }
    
    # 回归残差获取处理效应估计
    lm_fit1 <- lm(res_y1[!outliers1] ~ res_d1[!outliers1])
    lm_fit2 <- lm(res_y2[!outliers2] ~ res_d2[!outliers2])
    
    # 计算稳健标准误
    est1 <- coeftest(lm_fit1, vcov=vcovHC, type="HC0")
    est2 <- coeftest(lm_fit2, vcov=vcovHC, type="HC0")
    
    # 方法1: 分别回归并求加权平均 (DML1)
    w1 <- n1_eff / (n1_eff + n2_eff)
    w2 <- n2_eff / (n1_eff + n2_eff)
    
    b1 <- est1[2, 1]
    b2 <- est2[2, 1]
    be1 <- w1 * b1 + w2 * b2
    
    se1 <- est1[2, 2]
    se2 <- est2[2, 2]
    sig2 <- w1 * (se1^2 + (b1 - be1)^2) + w2 * (se2^2 + (b2 - be1)^2)
    se_avg <- sqrt(sig2)
    
    # 方法2: 合并残差进行回归 (DML2)
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
    cat("错误:", conditionMessage(e), "\n")
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

monte_carlo_dml_nn <- function(scenario, n, n_rep=1000, seed=123, size=1, maxit=300) {
  # 初始化结果存储矩阵
  results <- matrix(NA, nrow=n_rep, ncol=4)
  colnames(results) <- c("NN_est1", "NN_se1", "NN_est2", "NN_se2")
  
  # 设置全局随机种子
  set.seed(seed)
  
  # 重复实验n_rep次
  successful_runs <- 0
  
  for (i in 1:n_rep) {
    if (i %% 10 == 0) {
      cat(sprintf("情形 %d, 样本量 %d, 迭代 %d/%d\n", scenario, n, i, n_rep))
    }
    
    # 生成数据
    current_seed <- seed + i
    tryCatch({
      # 生成模拟数据
      sim_data <- generate_data(scenario=scenario, n=n, theta0=1)
      
      # 使用神经网络进行DML估计
      result_nn <- dml_nn(
        data=sim_data$data, 
        y=sim_data$y, 
        d=sim_data$d, 
        X=sim_data$X, 
        seed=current_seed,
        size=size,
        maxit=maxit
      )
      
      # 存储结果
      results[i, ] <- c(result_nn$estimate1, result_nn$se1, 
                        result_nn$estimate2, result_nn$se2)
      
      # 检查是否为有效结果
      if(!is.na(result_nn$estimate1) && !is.na(result_nn$estimate2)) {
        successful_runs <- successful_runs + 1
      }
    }, error=function(e) {
      cat("错误:", conditionMessage(e), "\n")
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
      Method = c("NN-Est1", "NN-Est2"),
      Mean = c(NA, NA),
      SD = c(NA, NA),
      Bias = c(NA, NA),
      RMSE = c(NA, NA),
      Success_Rate = 0
    ))
  }
  
  # 计算均值
  nn_mean1 <- mean(valid_results[, "NN_est1"])
  nn_mean2 <- mean(valid_results[, "NN_est2"])
  
  # 计算标准差
  nn_sd1 <- sd(valid_results[, "NN_est1"])
  nn_sd2 <- sd(valid_results[, "NN_est2"])
  
  # 计算偏差(Bias)
  nn_bias1 <- nn_mean1 - true_value
  nn_bias2 <- nn_mean2 - true_value
  
  # 计算均方误差(MSE)和均方根误差(RMSE)
  nn_mse1 <- mean((valid_results[, "NN_est1"] - true_value)^2)
  nn_mse2 <- mean((valid_results[, "NN_est2"] - true_value)^2)
  nn_rmse1 <- sqrt(nn_mse1)
  nn_rmse2 <- sqrt(nn_mse2)
  
  # 成功率
  success_rate <- sum(valid_rows) / nrow(results)
  
  # 组织结果
  summary_table <- data.frame(
    Method = c("NN-Est1", "NN-Est2"),
    Mean = c(nn_mean1, nn_mean2),
    SD = c(nn_sd1, nn_sd2),
    Bias = c(nn_bias1, nn_bias2),
    RMSE = c(nn_rmse1, nn_rmse2),
    Success_Rate = success_rate * 100  # 将成功率转换为百分比
  )
  
  return(summary_table)
}

##############################################
# 执行蒙特卡洛模拟
##############################################

# 定义样本容量和模拟情景
sample_sizes <- c(100, 200, 500, 1000, 2000, 5000, 10000)
scenarios <- 1:5  # 五种不同的数据生成过程

# 初始化最终结果存储
all_results <- list()

# 设置模拟参数 (为了速度，使用较小的样本和重复次数)
n_rep_demo <- 20  # 重复次数
nn_size <- 1      # 神经网络隐藏层节点数量
nn_maxit <- 300   # 神经网络最大迭代次数

# 执行所有情形和样本量的模拟
for (scenario in scenarios) {
  scenario_results <- list()
  
  cat("\n========== 开始情形", scenario, "的模拟 ==========\n")
  
  for (n in sample_sizes) {
    cat("\n开始样本量", n, "的模拟...\n")
    
    # 执行蒙特卡洛模拟
    start_time <- Sys.time()
    mc_results <- monte_carlo_dml_nn(
      scenario=scenario, 
      n=n, 
      n_rep=n_rep_demo, 
      seed=123, 
      size=nn_size, 
      maxit=nn_maxit
    )
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
  
  ##############################################
  # 可视化结果 - 九张图
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
  
  # 1. 绘制均值随样本量变化的图表
  mean_plot <- ggplot(plot_data, aes(x=SampleSize, y=Mean, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
    geom_line() +
    geom_point(size=3) +
    geom_hline(yintercept=1, linetype="dashed", color="black") +  # 添加真实值参考线
    labs(title="神经网络(NN)在不同情形和样本量下的处理效应估计",
         x="样本容量", y="估计值均值") +
    scale_x_continuous(breaks=sample_sizes) +
    theme_minimal() +
    facet_wrap(~ Scenario)
  
  # 2. 绘制RMSE随样本量变化的图表
  rmse_plot <- ggplot(plot_data, aes(x=SampleSize, y=RMSE, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
    geom_line() +
    geom_point(size=3) +
    labs(title="神经网络(NN)在不同情形和样本量下的均方根误差(RMSE)",
         x="样本容量", y="RMSE") +
    scale_x_continuous(breaks=sample_sizes) +
    theme_minimal() +
    facet_wrap(~ Scenario)
  
  # 3. 绘制偏差随样本量变化的图表
  bias_plot <- ggplot(plot_data, aes(x=SampleSize, y=Bias, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
    geom_line() +
    geom_point(size=3) +
    labs(title="神经网络(NN)在不同情形和样本量下的偏差(Bias)",
         x="样本容量", y="Bias") +
    geom_hline(yintercept=0, linetype="dashed", color="black") +  # 添加零偏差参考线
    scale_x_continuous(breaks=sample_sizes) +
    theme_minimal() +
    facet_wrap(~ Scenario)
  
  # 4. 绘制成功率随样本量变化的图表
  success_rate_plot <- ggplot(plot_data, aes(x=SampleSize, y=Success_Rate, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
    geom_line() +
    geom_point(size=3) +
    labs(title="神经网络(NN)在不同情形和样本量下的成功率",
         x="样本容量", y="成功率 (%)") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # 将成功率显示为百分比
    scale_x_continuous(breaks=sample_sizes) +
    theme_minimal() +
    facet_wrap(~ Scenario)
  
  # 5. 创建成功率热图
  # 首先确保数据格式正确
  heatmap_data <- plot_data %>%
    select(Method, Scenario, SampleSize, Success_Rate) %>%
    distinct()  # 确保没有重复
  
  success_heatmap <- ggplot(heatmap_data, aes(x=SampleSize, y=Method, fill=Success_Rate)) +
    geom_tile() +
    scale_fill_gradient(low="red", high="green", 
                        labels=scales::percent_format(scale=1),
                        limits=c(0, 100)) +
    labs(title="神经网络(NN)成功率热图",
         x="样本容量", y="方法", fill="成功率 (%)") +
    scale_x_continuous(breaks=sample_sizes) +
    theme_minimal() +
    facet_wrap(~ Scenario)
  
  # 6. 创建性能汇总图表
  performance_plot <- ggplot(plot_data, 
                             aes(x=SampleSize, y=RMSE, color=Method, shape=Scenario)) +
    geom_line(aes(group=interaction(Method, Scenario))) +
    geom_point(size=3) +
    labs(title="神经网络(NN)性能比较 - 所有情形",
         x="样本容量",
         y="均方根误差 (RMSE)") +
    scale_x_continuous(breaks=unique(plot_data$SampleSize)) +
    theme_minimal() +
    theme(legend.position="bottom")
  
  # 7. 创建成功率汇总图表（使用柱状图）
  success_plot <- ggplot(plot_data, 
                         aes(x=SampleSize, y=Success_Rate, fill=Scenario)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~ Method) +
    labs(title="神经网络(NN)成功率比较 - 所有情形",
         x="样本容量",
         y="成功率 (%)") +
    scale_y_continuous(labels=scales::percent_format(scale=1)) +
    theme_minimal() +
    theme(legend.position="bottom")
  
  # 8. 创建估计均值汇总图
  mean_summary_plot <- ggplot(plot_data, 
                              aes(x=SampleSize, y=Mean, color=Method, shape=Scenario)) +
    geom_line(aes(group=interaction(Method, Scenario))) +
    geom_point(size=3) +
    geom_hline(yintercept=1, linetype="dashed", color="black") +  # 添加真实值参考线
    labs(title="神经网络(NN)估计均值比较 - 所有情形",
         x="样本容量",
         y="估计均值") +
    scale_x_continuous(breaks=unique(plot_data$SampleSize)) +
    theme_minimal() +
    theme(legend.position="bottom")
  
  # 9. 创建估计偏差汇总图
  bias_summary_plot <- ggplot(plot_data, 
                              aes(x=SampleSize, y=Bias, color=Method, shape=Scenario)) +
    geom_line(aes(group=interaction(Method, Scenario))) +
    geom_point(size=3) +
    geom_hline(yintercept=0, linetype="dashed", color="black") +  # 添加零偏差参考线
    labs(title="神经网络(NN)估计偏差比较 - 所有情形",
         x="样本容量",
         y="估计偏差") +
    scale_x_continuous(breaks=unique(plot_data$SampleSize)) +
    theme_minimal() +
    theme(legend.position="bottom")
  
  # 保存所有九张图
  ggsave("dml_nn_mean.png", mean_plot, width=10, height=8)
  ggsave("dml_nn_rmse.png", rmse_plot, width=10, height=8)
  ggsave("dml_nn_bias.png", bias_plot, width=10, height=8)
  ggsave("dml_nn_success_rate.png", success_rate_plot, width=10, height=8)
  ggsave("dml_nn_success_heatmap.png", success_heatmap, width=12, height=6)
  ggsave("dml_nn_performance_summary.png", performance_plot, width=12, height=8)
  ggsave("dml_nn_success_summary.png", success_plot, width=12, height=8)
  ggsave("dml_nn_mean_summary.png", mean_summary_plot, width=12, height=8)
  ggsave("dml_nn_bias_summary.png", bias_summary_plot, width=12, height=8)
  
  # 导出汇总表格为CSV文件
  complete_summary <- data.frame()
  
  for (scenario in 1:5) {
    scenario_name <- paste0("Scenario_", scenario)
    scenario_results <- all_results[[paste0("scenario_", scenario)]]
    
    for (n in names(scenario_results)) {
      temp_df <- scenario_results[[n]]
      temp_df$Scenario <- scenario_name
      temp_df$SampleSize <- as.numeric(n)
      complete_summary <- rbind(complete_summary, temp_df)
    }
  }
  
  # 重新排列列的顺序
  complete_summary <- complete_summary[, c("Scenario", "SampleSize", "Method", "Mean", "SD", "Bias", "RMSE", "Success_Rate")]
  
  # 按情形、样本量和方法排序
  complete_summary <- complete_summary[order(complete_summary$Scenario, complete_summary$SampleSize, complete_summary$Method), ]
  
  # 导出为CSV文件
  write.csv(complete_summary, "dml_nn_complete_results.csv", row.names = FALSE)
  
  # 保存R数据
  save(all_results, plot_data, file="dml_nn_results.RData")
  
  # 导出为Excel文件
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- createWorkbook()
    addWorksheet(wb, "NN_Results")
    writeData(wb, "NN_Results", complete_summary)
    
    # 为每个情形创建单独的工作表
    for (scenario in 1:5) {
      scenario_data <- complete_summary[complete_summary$Scenario == paste0("Scenario_", scenario), ]
      addWorksheet(wb, paste0("Scenario_", scenario))
      writeData(wb, paste0("Scenario_", scenario), scenario_data)
    }
    
    saveWorkbook(wb, "dml_nn_results.xlsx", overwrite = TRUE)
  }
  
  cat("\n\n===============================================\n")
  cat("双重机器学习神经网络(NN)模拟完成！\n")
  cat("===============================================\n\n")
  
  cat("已创建以下文件:\n")
  cat("1. dml_nn_complete_results.csv - 长格式完整结果\n")
  cat("2. dml_nn_results.xlsx - Excel工作簿，包含多个工作表\n")
  cat("3. dml_nn_mean.png - 均值图表\n")
  cat("4. dml_nn_rmse.png - RMSE图表\n")
  cat("5. dml_nn_bias.png - 偏差图表\n")
  cat("6. dml_nn_success_rate.png - 成功率图表\n")
  cat("7. dml_nn_success_heatmap.png - 成功率热图\n")
  cat("8. dml_nn_performance_summary.png - 性能汇总图\n")
  cat("9. dml_nn_success_summary.png - 成功率汇总图\n")
  cat("10. dml_nn_mean_summary.png - 估计均值汇总图\n")
  cat("11. dml_nn_bias_summary.png - 估计偏差汇总图\n")
  cat("12. dml_nn_results.RData - R数据文件，包含所有结果\n\n")
  
  # 显示最终性能汇总
  cat("神经网络(NN)性能汇总:\n")
  summary_by_scenario <- aggregate(
    cbind(Mean, RMSE, Success_Rate) ~ Scenario + Method, 
    data=complete_summary, 
    FUN=mean
  )
  print(summary_by_scenario)
  
  # 返回完成信息
  cat("\n模拟过程已成功完成。\n")
}