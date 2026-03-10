##############################################
# 双重机器学习(DML)处理效应估计的蒙特卡洛模拟BAGGING
##############################################

# 清空环境
rm(list=ls())

# 加载必要的包
required_packages <- c("caret", "rpart", "sandwich", "lmtest", 
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
# 第二步：实现纯袋装法(Bagging)的双重机器学习
##############################################

# 自助抽样函数
bootstrap_sample <- function(data, size=nrow(data)) {
  idx <- sample(1:nrow(data), size=size, replace=TRUE)
  return(data[idx, ])
}

# 基于树的预测函数 - 用于回归问题(连续Y)
bagging_regression <- function(X, y, n_models=50, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  # 将X和y合并成一个数据框
  data <- data.frame(y=y, X)
  
  # 创建多个基础模型
  models <- list()
  
  for(i in 1:n_models) {
    # 自助抽样获取训练数据
    boot_data <- bootstrap_sample(data)
    
    # 拟合决策树模型
    tree_model <- rpart(y ~ ., data=boot_data, method="anova", 
                        control=rpart.control(cp=0.01, minsplit=10))
    
    # 存储模型
    models[[i]] <- tree_model
  }
  
  # 返回袋装模型集合
  return(models)
}

# 基于树的预测函数 - 用于分类问题(二分类D)
bagging_classification <- function(X, d, n_models=50, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  # 将X和d合并成一个数据框
  data <- data.frame(d=factor(d), X)
  
  # 创建多个基础模型
  models <- list()
  
  for(i in 1:n_models) {
    # 自助抽样获取训练数据
    boot_data <- bootstrap_sample(data)
    
    # 拟合决策树模型
    tree_model <- rpart(d ~ ., data=boot_data, method="class", 
                        control=rpart.control(cp=0.01, minsplit=10))
    
    # 存储模型
    models[[i]] <- tree_model
  }
  
  # 返回袋装模型集合
  return(models)
}

# 袋装集合模型预测函数 - 回归
predict_bagging_regression <- function(models, new_data) {
  # 计算每个模型的预测值
  predictions <- matrix(NA, nrow=nrow(new_data), ncol=length(models))
  
  for(i in 1:length(models)) {
    predictions[, i] <- predict(models[[i]], newdata=new_data)
  }
  
  # 计算平均预测值
  avg_predictions <- rowMeans(predictions, na.rm=TRUE)
  
  return(avg_predictions)
}

# 袋装集合模型预测函数 - 分类
predict_bagging_classification <- function(models, new_data) {
  # 计算每个模型的预测概率值
  predictions <- matrix(NA, nrow=nrow(new_data), ncol=length(models))
  
  for(i in 1:length(models)) {
    # 预测类别概率
    pred <- predict(models[[i]], newdata=new_data, type="prob")
    
    # 取第二列(类别为1的概率)，如果有两列
    if(ncol(pred) >= 2) {
      predictions[, i] <- pred[, 2]
    } else {
      # 如果只有一列，直接使用
      predictions[, i] <- pred
    }
  }
  
  # 计算平均预测概率
  avg_predictions <- rowMeans(predictions, na.rm=TRUE)
  
  return(avg_predictions)
}

# 双重机器学习 Bagging 函数
dml_bagging <- function(data, y, d, X, seed=123, n_models=50) {
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
  
  # 数据框转换
  X1_df <- as.data.frame(X1)
  X2_df <- as.data.frame(X2)
  
  tryCatch({
    # 第一部分数据上拟合Y~X和D~X
    y1_models <- bagging_regression(X1, y1, n_models=n_models, seed=seed)
    d1_models <- bagging_classification(X1, d1, n_models=n_models, seed=seed+1)
    
    # 预测第二部分数据
    yhat1 <- predict_bagging_regression(y1_models, X2_df)
    dhat1 <- predict_bagging_classification(d1_models, X2_df)
    dhat1 <- pmax(0.01, pmin(0.99, dhat1))  # 确保预测值在合理范围
    
    # 计算残差
    res_y1 <- y2 - yhat1
    res_d1 <- d2 - dhat1
    
    # 第二部分数据上拟合Y~X和D~X
    y2_models <- bagging_regression(X2, y2, n_models=n_models, seed=seed+2)
    d2_models <- bagging_classification(X2, d2, n_models=n_models, seed=seed+3)
    
    # 预测第一部分数据
    yhat2 <- predict_bagging_regression(y2_models, X1_df)
    dhat2 <- predict_bagging_classification(d2_models, X1_df)
    dhat2 <- pmax(0.01, pmin(0.99, dhat2))  # 确保预测值在合理范围
    
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

monte_carlo_dml_bagging <- function(scenario, n, n_rep=1000, seed=123, n_models=50) {
  # 初始化结果存储矩阵
  results <- matrix(NA, nrow=n_rep, ncol=4)
  colnames(results) <- c("Bagging_est1", "Bagging_se1", "Bagging_est2", "Bagging_se2")
  
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
      
      # 使用袋装法进行DML估计
      result_bagging <- dml_bagging(
        data=sim_data$data, 
        y=sim_data$y, 
        d=sim_data$d, 
        X=sim_data$X, 
        seed=current_seed,
        n_models=n_models
      )
      
      # 存储结果
      results[i, ] <- c(result_bagging$estimate1, result_bagging$se1, 
                        result_bagging$estimate2, result_bagging$se2)
      
      # 检查是否为有效结果
      if(!is.na(result_bagging$estimate1) && !is.na(result_bagging$estimate2)) {
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
      Method = c("Bagging-Est1", "Bagging-Est2"),
      Mean = c(NA, NA),
      SD = c(NA, NA),
      Bias = c(NA, NA),
      RMSE = c(NA, NA),
      Success_Rate = 0
    ))
  }
  
  # 计算均值
  bagging_mean1 <- mean(valid_results[, "Bagging_est1"])
  bagging_mean2 <- mean(valid_results[, "Bagging_est2"])
  
  # 计算标准差
  bagging_sd1 <- sd(valid_results[, "Bagging_est1"])
  bagging_sd2 <- sd(valid_results[, "Bagging_est2"])
  
  # 计算偏差(Bias)
  bagging_bias1 <- bagging_mean1 - true_value
  bagging_bias2 <- bagging_mean2 - true_value
  
  # 计算均方误差(MSE)和均方根误差(RMSE)
  bagging_mse1 <- mean((valid_results[, "Bagging_est1"] - true_value)^2)
  bagging_mse2 <- mean((valid_results[, "Bagging_est2"] - true_value)^2)
  bagging_rmse1 <- sqrt(bagging_mse1)
  bagging_rmse2 <- sqrt(bagging_mse2)
  
  # 成功率
  success_rate <- sum(valid_rows) / nrow(results)
  
  # 组织结果
  summary_table <- data.frame(
    Method = c("Bagging-Est1", "Bagging-Est2"),
    Mean = c(bagging_mean1, bagging_mean2),
    SD = c(bagging_sd1, bagging_sd2),
    Bias = c(bagging_bias1, bagging_bias2),
    RMSE = c(bagging_rmse1, bagging_rmse2),
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

# 设置模拟参数 (为了速度，使用较小的样本和重复次数)
n_rep_demo <- 20  # 重复次数
n_models <- 30    # 每个袋装模型的树数量

# 执行所有情形和样本量的模拟
for (scenario in scenarios) {
  scenario_results <- list()
  
  cat("\n========== 开始情形", scenario, "的模拟 ==========\n")
  
  for (n in sample_sizes) {
    cat("\n开始样本量", n, "的模拟...\n")
    
    # 执行蒙特卡洛模拟
    start_time <- Sys.time()
    mc_results <- monte_carlo_dml_bagging(
      scenario=scenario, 
      n=n, 
      n_rep=n_rep_demo, 
      seed=123, 
      n_models=n_models
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
}

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
  labs(title="袋装法(Bagging)在不同情形和样本量下的处理效应估计",
       x="样本容量", y="估计值均值") +
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 2. 绘制RMSE随样本量变化的图表
rmse_plot <- ggplot(plot_data, aes(x=SampleSize, y=RMSE, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="袋装法(Bagging)在不同情形和样本量下的均方根误差(RMSE)",
       x="样本容量", y="RMSE") +
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 3. 绘制偏差随样本量变化的图表
bias_plot <- ggplot(plot_data, aes(x=SampleSize, y=Bias, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="袋装法(Bagging)在不同情形和样本量下的偏差(Bias)",
       x="样本容量", y="Bias") +
  geom_hline(yintercept=0, linetype="dashed", color="black") +  # 添加零偏差参考线
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 4. 绘制成功率随样本量变化的图表
success_rate_plot <- ggplot(plot_data, aes(x=SampleSize, y=Success_Rate, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="袋装法(Bagging)在不同情形和样本量下的成功率",
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
  labs(title="袋装法(Bagging)成功率热图",
       x="样本容量", y="方法", fill="成功率 (%)") +
  scale_x_continuous(breaks=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 保存结果和图表
ggsave("dml_bagging_mean.png", mean_plot, width=10, height=8)
ggsave("dml_bagging_rmse.png", rmse_plot, width=10, height=8)
ggsave("dml_bagging_bias.png", bias_plot, width=10, height=8)
ggsave("dml_bagging_success_rate.png", success_rate_plot, width=10, height=8)
ggsave("dml_bagging_success_heatmap.png", success_heatmap, width=12, height=6)
save(all_results, plot_data, file="dml_bagging_results.RData")

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
write.csv(complete_summary, "dml_bagging_complete_results.csv", row.names = FALSE)

# 导出为Excel文件
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- createWorkbook()
  addWorksheet(wb, "Bagging_Results")
  writeData(wb, "Bagging_Results", complete_summary)
  
  # 为每个情形创建单独的工作表
  for (scenario in 1:5) {
    scenario_data <- complete_summary[complete_summary$Scenario == paste0("Scenario_", scenario), ]
    addWorksheet(wb, paste0("Scenario_", scenario))
    writeData(wb, paste0("Scenario_", scenario), scenario_data)
  }
  
  saveWorkbook(wb, "dml_bagging_results.xlsx", overwrite = TRUE)
}

# 6. 创建性能汇总图表
performance_plot <- ggplot(complete_summary, 
                           aes(x=SampleSize, y=RMSE, color=Method, shape=Scenario)) +
  geom_line(aes(group=interaction(Method, Scenario))) +
  geom_point(size=3) +
  labs(title="袋装法(Bagging)性能比较 - 所有情形",
       x="样本容量",
       y="均方根误差 (RMSE)") +
  scale_x_continuous(breaks=unique(complete_summary$SampleSize)) +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("dml_bagging_performance_summary.png", performance_plot, width=12, height=8)

# 7. 创建成功率汇总图表（使用柱状图）
success_plot <- ggplot(complete_summary, 
                       aes(x=SampleSize, y=Success_Rate, fill=Scenario)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~ Method) +
  labs(title="袋装法(Bagging)成功率比较 - 所有情形",
       x="样本容量",
       y="成功率 (%)") +
  scale_y_continuous(labels=scales::percent_format(scale=1)) +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("dml_bagging_success_summary.png", success_plot, width=12, height=8)

# 8. 创建估计均值汇总图
mean_summary_plot <- ggplot(complete_summary, 
                            aes(x=SampleSize, y=Mean, color=Method, shape=Scenario)) +
  geom_line(aes(group=interaction(Method, Scenario))) +
  geom_point(size=3) +
  geom_hline(yintercept=1, linetype="dashed", color="black") +  # 添加真实值参考线
  labs(title="袋装法(Bagging)估计均值比较 - 所有情形",
       x="样本容量",
       y="估计均值") +
  scale_x_continuous(breaks=unique(complete_summary$SampleSize)) +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("dml_bagging_mean_summary.png", mean_summary_plot, width=12, height=8)

# 9. 创建估计偏差汇总图
bias_summary_plot <- ggplot(complete_summary, 
                            aes(x=SampleSize, y=Bias, color=Method, shape=Scenario)) +
  geom_line(aes(group=interaction(Method, Scenario))) +
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +  # 添加零偏差参考线
  labs(title="袋装法(Bagging)估计偏差比较 - 所有情形",
       x="样本容量",
       y="估计偏差") +
  scale_x_continuous(breaks=unique(complete_summary$SampleSize)) +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("dml_bagging_bias_summary.png", bias_summary_plot, width=12, height=8)

cat("\n\n===============================================\n")
cat("双重机器学习袋装法(Bagging)模拟完成！\n")
cat("===============================================\n\n")

cat("已创建以下文件:\n")
cat("1. dml_bagging_complete_results.csv - 长格式完整结果\n")
cat("2. dml_bagging_results.xlsx - Excel工作簿，包含多个工作表\n")
cat("3. dml_bagging_mean.png - 均值图表\n")
cat("4. dml_bagging_rmse.png - RMSE图表\n")
cat("5. dml_bagging_bias.png - 偏差图表\n")
cat("6. dml_bagging_success_rate.png - 成功率图表\n")
cat("7. dml_bagging_success_heatmap.png - 成功率热图\n")
cat("8. dml_bagging_performance_summary.png - 性能汇总图\n")
cat("9. dml_bagging_success_summary.png - 成功率汇总图\n")
cat("10. dml_bagging_mean_summary.png - 估计均值汇总图\n")
cat("11. dml_bagging_bias_summary.png - 估计偏差汇总图\n")
cat("12. dml_bagging_results.RData - R数据文件，包含所有结果\n\n")

# 显示最终性能汇总
cat("袋装法(Bagging)性能汇总:\n")
summary_by_scenario <- aggregate(
  cbind(Mean, RMSE, Success_Rate) ~ Scenario + Method, 
  data=complete_summary, 
  FUN=mean
)
print(summary_by_scenario)

# 返回完成信息
cat("\n模拟过程已成功完成。\n")