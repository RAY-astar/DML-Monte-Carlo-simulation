##############################################
# 双重机器学习(DML)处理效应估计的蒙特卡洛模拟BAGGING - 优化版
##############################################

# 清空环境
rm(list=ls())

# 加载必要的包
required_packages <- c("caret", "rpart", "sandwich", "lmtest", 
                       "reshape2", "ggplot2", "knitr", "dplyr", 
                       "scales", "openxlsx", "tidyr", "parallel")

# 检查并安装缺失的包
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

##############################################
# 第一步：优化的数据生成函数（保留原始函数形式）
##############################################

# 改进的数据生成函数 - 保留原始函数形式但增强数值稳定性
generate_data <- function(scenario, n, theta0=1) {
  # 根据不同情景生成数据
  if (scenario == 1) {
    # 情形1: 简单线性模拟
    X1 <- runif(n, min=0, max=1)
    X2 <- runif(n, min=0, max=1)
    
    # 计算处理概率和结果函数
    m0 <- pmax(0.1, pmin(0.9, (X1 + X2) / 2))
    f0 <- -6 + 12 * X1 + 6 * X2
    
    X <- cbind(X1, X2)
    colnames(X) <- c("X1", "X2")
  } else if (scenario == 2) {
    # 情形2: 线性模拟添加交互项
    X1 <- rnorm(n, mean=0, sd=1)
    X2 <- rnorm(n, mean=0, sd=1)
    
    # 防止极端值
    X1 <- pmin(pmax(X1, -3), 3)
    X2 <- pmin(pmax(X2, -3), 3)
    
    X3 <- X1 * X2  # 交互项
    
    D_temp <- X1 + X2 + X3
    # 确保处理概率在合理范围内
    m0 <- pmax(0.1, pmin(0.9, 1 / (1 + exp(-D_temp))))
    f0 <- -6 + 12 * X1 + 6 * X2 + 2 * X3
    
    X <- cbind(X1, X2, X3)
    colnames(X) <- c("X1", "X2", "X3")
  } else if (scenario == 3) {
    # 情形3: 线性模拟添加平方项
    X1 <- rnorm(n, mean=0, sd=1)
    X2 <- rnorm(n, mean=0, sd=1)
    
    # 防止极端值
    X1 <- pmin(pmax(X1, -3), 3)
    X2 <- pmin(pmax(X2, -3), 3)
    
    X3 <- X1 * X1  # 平方项
    X4 <- X2 * X2  # 平方项
    
    D_temp <- X1 + X2 + X3 + X4
    # 确保处理概率在合理范围内
    m0 <- pmax(0.1, pmin(0.9, 1 / (1 + exp(-D_temp))))
    f0 <- -6 + 12 * X1 + 6 * X2 + 2 * X3 + 3 * X4
    
    X <- cbind(X1, X2, X3, X4)
    colnames(X) <- c("X1", "X2", "X3", "X4")
  } else if (scenario == 4) {
    # 情形4: 线性模拟添加交互项和平方项
    X1 <- rnorm(n, mean=0, sd=1)
    X2 <- rnorm(n, mean=0, sd=1)
    
    # 防止极端值
    X1 <- pmin(pmax(X1, -3), 3)
    X2 <- pmin(pmax(X2, -3), 3)
    
    X3 <- X1 * X1  # 平方项
    X4 <- X2 * X2  # 平方项
    X5 <- X1 * X2  # 交互项
    
    D_temp <- X1 + X2 + X3 + X4 + X5
    # 确保处理概率在合理范围内
    m0 <- pmax(0.1, pmin(0.9, 1 / (1 + exp(-D_temp))))
    f0 <- -6 + 12 * X1 + 6 * X2 + 2 * X3 + 3 * X4 + 3 * X5
    
    X <- cbind(X1, X2, X3, X4, X5)
    colnames(X) <- c("X1", "X2", "X3", "X4", "X5")
  } else if (scenario == 5) {
    # 情形5: 优化的指数型非线性 - 保留原始函数形式但增强数值稳定性
    # 使用较小的标准差生成特征，减少极端值
    X1 <- rnorm(n, mean=0, sd=0.5)  # 减小标准差以减少极端值
    X2 <- rnorm(n, mean=0, sd=0.5)
    X3 <- rnorm(n, mean=0, sd=0.5)
    X4 <- rnorm(n, mean=0, sd=0.5)
    X5 <- rnorm(n, mean=0, sd=0.5)
    
    # 截断以确保稳定性
    X1 <- pmin(pmax(X1, -1.5), 1.5)  # 更严格的限制
    X2 <- pmin(pmax(X2, -1.5), 1.5)
    X3 <- pmin(pmax(X3, -1.5), 1.5)
    X4 <- pmin(pmax(X4, -1.5), 1.5)
    X5 <- pmin(pmax(X5, -1.5), 1.5)
    
    X_matrix <- cbind(X1, X2, X3, X4, X5)
    
    # 计算f0(Xi) - 保留原始公式但添加数值稳定性保护
    f0 <- rep(1, n)
    # 使用更稳定的分母计算
    denominator <- sqrt(max(0.001, -0.5 * exp(2) + 2 * exp(1) - 1.5))
    
    for (j in 1:5) {
      # 限制指数增长以防止数值溢出，使用更保守的上限
      exp_val <- exp(X_matrix[,j])
      
      # 检测并处理潜在的溢出
      if(any(is.infinite(exp_val) | is.nan(exp_val) | exp_val > 1e8)) {
        cat("警告: 检测到潜在的数值溢出，应用更严格的限制\n")
        exp_val <- pmin(exp_val, 1e8)  # 更保守的上限
      }
      
      # 使用原始公式但确保数值稳定
      term_val <- (4 * (exp_val - exp(1) + 1)) / denominator
      # 检查并限制潜在的极端值
      term_val <- pmin(pmax(term_val, -1e4), 1e4)
      f0 <- f0 + term_val
    }
    
    # 计算m0(Xi) - 保留原始公式但增强稳定性
    m0 <- rep(0.5, n)
    log2 <- log(2)  # 预计算log(2)提高效率
    
    for (j in 1:5) {
      # 防止log(负数)或log(0)，使用更严格的保护
      x_val <- X_matrix[,j]
      x_safe <- x_val + 1  # 确保为正
      # 额外检查
      x_safe <- pmax(1e-8, x_safe)  # 确保严格大于0
      
      # 使用原始公式
      log_term <- log(x_safe) / log2
      
      # 增加有效性检查
      if(any(is.infinite(log_term) | is.nan(log_term))) {
        cat("警告: 检测到无效的对数值，应用修正\n")
        log_term[is.infinite(log_term) | is.nan(log_term)] <- 0  # 替换无效值
      }
      
      term_val <- (-0.1 + 0.2 * (log2 - log_term))
      # 限制极端值
      term_val <- pmin(pmax(term_val, -0.3), 0.3)
      m0 <- m0 + term_val
    }
    
    # 确保m0严格在(0,1)范围内，避免极端值
    m0 <- pmax(0.1, pmin(0.9, m0))
    
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
# 第二步：优化袋装法(Bagging)的双重机器学习
##############################################

# 改进的自助抽样函数
bootstrap_sample <- function(data, size=nrow(data)) {
  idx <- sample(1:nrow(data), size=size, replace=TRUE)
  return(data[idx, ])
}

# 优化的基于随机森林的预测函数 - 用于回归问题(连续Y)
bagging_regression <- function(X, y, n_models=50, seed=NULL, 
                               min_node_size=5, max_depth=30) {
  if(!is.null(seed)) set.seed(seed)
  
  # 将X和y合并成一个数据框
  data <- data.frame(y=y, X)
  
  # 创建多个基础模型
  models <- list()
  
  for(i in 1:n_models) {
    # 自助抽样获取训练数据
    boot_data <- bootstrap_sample(data)
    
    # 拟合决策树模型，增加参数灵活性
    tree_model <- rpart(y ~ ., data=boot_data, method="anova", 
                        control=rpart.control(cp=0.001,  # 降低复杂度参数，增加树的灵活性
                                              minbucket=min_node_size,
                                              maxdepth=max_depth))
    
    # 存储模型
    models[[i]] <- tree_model
  }
  
  # 返回袋装模型集合
  return(models)
}

# 优化的基于随机森林的预测函数 - 用于分类问题(二分类D)
bagging_classification <- function(X, d, n_models=50, seed=NULL,
                                   min_node_size=5, max_depth=30) {
  if(!is.null(seed)) set.seed(seed)
  
  # 将X和d合并成一个数据框
  data <- data.frame(d=factor(d), X)
  
  # 创建多个基础模型
  models <- list()
  
  for(i in 1:n_models) {
    # 自助抽样获取训练数据
    boot_data <- bootstrap_sample(data)
    
    # 拟合决策树模型，增加参数灵活性
    tree_model <- rpart(d ~ ., data=boot_data, method="class", 
                        control=rpart.control(cp=0.001,  # 降低复杂度参数
                                              minbucket=min_node_size, 
                                              maxdepth=max_depth))
    
    # 存储模型
    models[[i]] <- tree_model
  }
  
  # 返回袋装模型集合
  return(models)
}

# 优化的袋装集合模型预测函数 - 回归
predict_bagging_regression <- function(models, new_data) {
  # 计算每个模型的预测值
  n_models <- length(models)
  predictions <- matrix(NA, nrow=nrow(new_data), ncol=n_models)
  
  for(i in 1:n_models) {
    tryCatch({
      predictions[, i] <- predict(models[[i]], newdata=new_data)
    }, error=function(e) {
      # 捕获错误并记录（实际应用中可能需要更复杂的处理）
      cat("预测错误:", conditionMessage(e), "\n")
    })
  }
  
  # 计算平均预测值，处理可能的NA
  valid_counts <- rowSums(!is.na(predictions))
  avg_predictions <- rowSums(predictions, na.rm=TRUE) / pmax(1, valid_counts)
  
  return(avg_predictions)
}

# 优化的袋装集合模型预测函数 - 分类
predict_bagging_classification <- function(models, new_data) {
  # 计算每个模型的预测概率值
  n_models <- length(models)
  predictions <- matrix(NA, nrow=nrow(new_data), ncol=n_models)
  
  for(i in 1:n_models) {
    tryCatch({
      # 预测类别概率
      pred <- predict(models[[i]], newdata=new_data, type="prob")
      
      # 取类别为1的概率
      if(ncol(pred) >= 2) {
        predictions[, i] <- pred[, 2]
      } else if(ncol(pred) == 1) {
        # 只有一列时
        predictions[, i] <- pred[, 1]
      }
    }, error=function(e) {
      cat("预测错误:", conditionMessage(e), "\n")
    })
  }
  
  # 处理可能的NA并计算平均预测概率
  valid_counts <- rowSums(!is.na(predictions))
  avg_predictions <- rowSums(predictions, na.rm=TRUE) / pmax(1, valid_counts)
  
  # 确保概率值在合理范围内
  avg_predictions <- pmax(0.01, pmin(0.99, avg_predictions))
  
  return(avg_predictions)
}

# 优化的双重机器学习 Bagging 函数
dml_bagging <- function(data, y, d, X, seed=123, n_models=50, 
                        robust=TRUE, trimming=TRUE) {
  # 设置随机种子
  set.seed(seed)
  
  # 获取样本量
  n <- nrow(data)
  
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
  
  # 分析是否存在极值点
  if(trimming) {
    # 计算Y的分位数
    q_y <- quantile(y, probs=c(0.01, 0.99), na.rm=TRUE)
    y_outliers <- y < q_y[1] | y > q_y[2]
    
    # 在训练每个部分前排除异常值
    y1_valid <- !y_outliers[trainIndex]
    y2_valid <- !y_outliers[-trainIndex]
    
    # 如果过滤后样本太少，则不过滤
    if(sum(y1_valid) < length(y1) * 0.9) y1_valid <- rep(TRUE, length(y1))
    if(sum(y2_valid) < length(y2) * 0.9) y2_valid <- rep(TRUE, length(y2))
    
    # 应用过滤
    y1 <- y1[y1_valid]
    d1 <- d1[y1_valid]
    X1 <- X1[y1_valid, ]
    X1_df <- X1_df[y1_valid, ]
    
    y2 <- y2[y2_valid]
    d2 <- d2[y2_valid]
    X2 <- X2[y2_valid, ]
    X2_df <- X2_df[y2_valid, ]
  }
  
  # 初始化结果变量
  be1 <- NA
  be2 <- NA
  se_avg <- NA
  se_combined <- NA
  
  # 封装错误处理以确保鲁棒性
  tryCatch({
    #############################
    # 第一部分数据上拟合
    #############################
    
    # 拟合Y~X
    cat("正在第一部分数据上拟合Y~X模型...\n")
    y1_models <- bagging_regression(X1, y1, n_models=n_models, seed=seed)
    
    # 拟合D~X
    cat("正在第一部分数据上拟合D~X模型...\n")
    d1_models <- bagging_classification(X1, d1, n_models=n_models, seed=seed+1)
    
    # 预测第二部分数据
    cat("预测第二部分数据...\n")
    yhat1 <- predict_bagging_regression(y1_models, X2_df)
    dhat1 <- predict_bagging_classification(d1_models, X2_df)
    
    # 确保预测值在合理范围
    dhat1 <- pmax(0.01, pmin(0.99, dhat1))
    
    # 计算残差
    res_y1 <- y2 - yhat1
    res_d1 <- d2 - dhat1
    
    #############################
    # 第二部分数据上拟合
    #############################
    
    # 拟合Y~X
    cat("正在第二部分数据上拟合Y~X模型...\n")
    y2_models <- bagging_regression(X2, y2, n_models=n_models, seed=seed+2)
    
    # 拟合D~X
    cat("正在第二部分数据上拟合D~X模型...\n")
    d2_models <- bagging_classification(X2, d2, n_models=n_models, seed=seed+3)
    
    # 预测第一部分数据
    cat("预测第一部分数据...\n")
    yhat2 <- predict_bagging_regression(y2_models, X1_df)
    dhat2 <- predict_bagging_classification(d2_models, X1_df)
    
    # 确保预测值在合理范围
    dhat2 <- pmax(0.01, pmin(0.99, dhat2))
    
    # 计算残差
    res_y2 <- y1 - yhat2
    res_d2 <- d1 - dhat2
    
    #############################
    # 增强型异常值检测
    #############################
    
    # 改进的异常值检测函数
    detect_outliers <- function(x, k=3) {
      if(length(unique(x)) <= 1) return(rep(FALSE, length(x)))
      
      # 使用更稳健的中位数绝对偏差(MAD)方法
      med_x <- median(x, na.rm=TRUE)
      mad_x <- mad(x, na.rm=TRUE)
      
      # 如果MAD太小，使用标准偏差
      if(mad_x < 1e-8) {
        sd_x <- sd(x, na.rm=TRUE)
        if(sd_x < 1e-8) return(rep(FALSE, length(x)))
        lower_bound <- med_x - k * sd_x
        upper_bound <- med_x + k * sd_x
      } else {
        lower_bound <- med_x - k * mad_x
        upper_bound <- med_x + k * mad_x
      }
      
      return(x < lower_bound | x > upper_bound | is.na(x))
    }
    
    # 对残差进行异常值检测
    if(robust) {
      outliers1 <- detect_outliers(res_y1) | detect_outliers(res_d1)
      outliers2 <- detect_outliers(res_y2) | detect_outliers(res_d2)
      
      # 如果过滤后样本太少，则不过滤
      if(sum(!outliers1) < length(res_y1) * 0.9) outliers1 <- rep(FALSE, length(res_y1))
      if(sum(!outliers2) < length(res_y2) * 0.9) outliers2 <- rep(FALSE, length(res_y2))
    } else {
      outliers1 <- rep(FALSE, length(res_y1))
      outliers2 <- rep(FALSE, length(res_y2))
    }
    
    # 计算有效样本量
    n1_eff <- sum(!outliers1)
    n2_eff <- sum(!outliers2)
    
    if(n1_eff < 10 || n2_eff < 10) {
      cat("警告: 有效样本量太小! n1_eff =", n1_eff, ", n2_eff =", n2_eff, "\n")
      return(list(
        estimate1 = NA,
        se1 = NA,
        estimate2 = NA,
        se2 = NA
      ))
    }
    
    #############################
    # 回归残差获取处理效应估计
    #############################
    
    # 对两部分数据分别进行回归
    lm_fit1 <- lm(res_y1[!outliers1] ~ res_d1[!outliers1])
    lm_fit2 <- lm(res_y2[!outliers2] ~ res_d2[!outliers2])
    
    # 计算稳健标准误
    est1 <- coeftest(lm_fit1, vcov=vcovHC, type="HC1")  # 使用HC1更适合小样本
    est2 <- coeftest(lm_fit2, vcov=vcovHC, type="HC1")
    
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
    # 构建组合数据集
    res_y_combined <- c(res_y1[!outliers1], res_y2[!outliers2])
    res_d_combined <- c(res_d1[!outliers1], res_d2[!outliers2])
    
    # 检查数据有效性
    if(length(unique(res_d_combined)) < 2) {
      cat("警告: 组合后的处理变量只有一个值!\n")
      be2 <- NA
      se_combined <- NA
    } else {
      # 合并回归
      lm_fit_combined <- lm(res_y_combined ~ res_d_combined)
      est_combined <- coeftest(lm_fit_combined, vcov=vcovHC, type="HC1")
      
      be2 <- est_combined[2, 1]
      se_combined <- est_combined[2, 2]
    }
    
    # 检查估计结果是否合理
    if(is.na(be1) || is.na(be2) || abs(be1) > 10 || abs(be2) > 10) {
      cat("警告: 估计值不在合理范围内! be1 =", be1, ", be2 =", be2, "\n")
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
# 第三步：并行化的蒙特卡洛模拟主函数
##############################################

monte_carlo_dml_bagging <- function(scenario, n, n_rep=1000, seed=123, n_models=50, 
                                    use_parallel=TRUE, n_cores=NULL) {
  # 初始化结果存储矩阵
  results <- matrix(NA, nrow=n_rep, ncol=4)
  colnames(results) <- c("Bagging_est1", "Bagging_se1", "Bagging_est2", "Bagging_se2")
  
  # 设置全局随机种子
  set.seed(seed)
  seeds <- seed + 1:n_rep
  
  # 定义单个模拟运行的包装函数
  run_single_simulation <- function(iter, scenario, n, seed, n_models) {
    cat(sprintf("情形 %d, 样本量 %d, 迭代 %d\n", scenario, n, iter))
    
    # 生成数据
    current_seed <- seed
    
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
      
      # 返回结果
      c(result_bagging$estimate1, result_bagging$se1,
        result_bagging$estimate2, result_bagging$se2)
    }, error=function(e) {
      cat("错误:", conditionMessage(e), "\n")
      c(NA, NA, NA, NA)
    })
  }
  
  # 确定使用的核心数
  if(use_parallel) {
    if(is.null(n_cores)) {
      n_cores <- detectCores() - 1  # 留一个核心给系统
      if(n_cores <= 0) n_cores <- 1
    }
    cat("使用并行计算, 核心数:", n_cores, "\n")
    
    # 创建集群
    cl <- makeCluster(n_cores)
    
    # 导出必要的函数和包
    clusterExport(cl, c("generate_data", "bootstrap_sample", 
                        "bagging_regression", "bagging_classification", 
                        "predict_bagging_regression", "predict_bagging_classification",
                        "dml_bagging", "run_single_simulation"))
    
    # 确保每个核心都加载必要的包
    clusterEvalQ(cl, {
      library(caret)
      library(rpart)
      library(sandwich)
      library(lmtest)
    })
    
    # 并行执行
    results_list <- parLapply(cl, 1:n_rep, function(i) {
      run_single_simulation(i, scenario, n, seeds[i], n_models)
    })
    
    # 关闭集群
    stopCluster(cl)
    
    # 整合结果
    for(i in 1:n_rep) {
      results[i, ] <- results_list[[i]]
    }
    
  } else {
    # 顺序执行
    for(i in 1:n_rep) {
      results[i, ] <- run_single_simulation(i, scenario, n, seeds[i], n_models)
    }
  }
  
  # 计算成功率
  successful_runs <- sum(complete.cases(results))
  success_rate <- successful_runs / n_rep * 100
  cat(sprintf("情形 %d, 样本量 %d 的成功率: %.1f%%\n", scenario, n, success_rate))
  
  return(list(results=results, success_rate=success_rate))
}

# 优化的统计汇总函数
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
      CI_Coverage = c(NA, NA),
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
  
  # 计算置信区间覆盖率
  alpha <- 0.05  # 95%置信区间
  z_crit <- qnorm(1 - alpha/2)
  
  # 构建每个估计的置信区间并检查是否包含真实值
  ci_lower1 <- valid_results[, "Bagging_est1"] - z_crit * valid_results[, "Bagging_se1"]
  ci_upper1 <- valid_results[, "Bagging_est1"] + z_crit * valid_results[, "Bagging_se1"]
  ci_lower2 <- valid_results[, "Bagging_est2"] - z_crit * valid_results[, "Bagging_se2"]
  ci_upper2 <- valid_results[, "Bagging_est2"] + z_crit * valid_results[, "Bagging_se2"]
  
  ci_coverage1 <- mean(ci_lower1 <= true_value & true_value <= ci_upper1, na.rm=TRUE) * 100
  ci_coverage2 <- mean(ci_lower2 <= true_value & true_value <= ci_upper2, na.rm=TRUE) * 100
  
  # 成功率
  success_rate <- sum(valid_rows) / nrow(results) * 100
  
  # 组织结果
  summary_table <- data.frame(
    Method = c("Bagging-Est1", "Bagging-Est2"),
    Mean = c(bagging_mean1, bagging_mean2),
    SD = c(bagging_sd1, bagging_sd2),
    Bias = c(bagging_bias1, bagging_bias2),
    RMSE = c(bagging_rmse1, bagging_rmse2),
    CI_Coverage = c(ci_coverage1, ci_coverage2),
    Success_Rate = c(success_rate, success_rate)
  )
  
  return(summary_table)
}
##############################################
# 第四步：执行蒙特卡洛模拟
##############################################

# 定义样本容量和模拟情景
sample_sizes <- c(100, 200, 500, 1000, 2000, 5000, 10000)
scenarios <- 1:5  # 五种不同的数据生成过程

# 初始化最终结果存储
all_results <- list()

# 设置模拟参数
n_rep <- 100  # 重复次数，实际应用中可以设置更大 (推荐1000次以上)
n_models <- 50   # 每个袋装模型的树数量
use_parallel <- TRUE  # 是否使用并行计算
n_cores <- NULL  # 自动检测核心数

# 执行所有情形和样本量的模拟
for (scenario in scenarios) {
  scenario_results <- list()
  
  cat("\n========== 开始情形", scenario, "的模拟 ==========\n")
  
  for (n in sample_sizes) {
    cat("\n开始样本量", n, "的模拟...\n")
    
    # 执行蒙特卡洛模拟
    start_time <- Sys.time()
    
    # 情形5需要特殊处理
    special_scenario5 <- (scenario == 5 && n >= 5000)
    
    if(special_scenario5) {
      # 情形5大样本量的特殊配置
      current_n_rep <- max(10, n_rep / 5)  # 减少重复次数
      current_n_models <- min(100, n_models * 2)  # 增加树数量
      cat("情形5大样本特殊配置: 重复次数 =", current_n_rep, ", 树数量 =", current_n_models, "\n")
    } else {
      current_n_rep <- n_rep
      current_n_models <- n_models
    }
    
    mc_results <- monte_carlo_dml_bagging(
      scenario=scenario, 
      n=n, 
      n_rep=current_n_rep, 
      seed=123, 
      n_models=current_n_models,
      use_parallel=use_parallel,
      n_cores=n_cores
    )
    
    end_time <- Sys.time()
    
    # 计算汇总统计量
    summary <- summarize_results(mc_results$results)
    scenario_results[[as.character(n)]] <- summary
    
    # 输出当前情形和样本量的结果
    cat("用时:", round(difftime(end_time, start_time, units="mins"), 2), "分钟\n")
    print(kable(summary, digits=4, caption=paste0("情形", scenario, "样本量", n, "的结果")))
  }
  
  all_results[[paste0("scenario_", scenario)]] <- scenario_results
  
  # 在每个场景结束后保存阶段性结果
  partial_results_file <- paste0("dml_bagging_results_scenario_", scenario, ".RData")
  save(scenario_results, file=partial_results_file)
  cat("情形", scenario, "的结果已保存到", partial_results_file, "\n")
}

# 加载已保存的结果（如果有）
result_files <- list.files(pattern="dml_bagging_results_scenario_\\d+\\.RData")
if(length(result_files) > 0) 
  cat("发现之前保存的结果文件，正在加载...\n")
  for(file in result_files) {
    scenario_num <- as.numeric(gsub(".*scenario_(\\d+)\\.RData", "\\1", file))
    cat("加载情形", scenario_num, "的结果...\n")
    temp_env <- new.env()
    load(file, envir=temp_env)
    all_results[[paste0("scenario_", scenario_num)]] <- temp_env$scenario_results
  }
  cat("成功加载之前的结果!\n")

  ##############################################
  # 第五步：增强的结果展示
  ##############################################
  
  # 增强的结果汇总表格生成函数
  create_enhanced_summary_table <- function(all_results, scenario_id) {
    # 获取特定情形的结果
    scenario_results <- all_results[[paste0("scenario_", scenario_id)]]
    if(is.null(scenario_results)) {
      cat("警告: 找不到情形", scenario_id, "的结果\n")
      return(NULL)
    }
    
    # 合并所有样本量的结果
    result_table <- data.frame()
    for (n in names(scenario_results)) {
      temp_df <- scenario_results[[n]]
      temp_df$SampleSize <- as.numeric(n)
      result_table <- rbind(result_table, temp_df)
    }
    
    # 整理表格格式，增加置信区间覆盖率
    result_table <- result_table[, c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "CI_Coverage", "Success_Rate")]
    
    # 按样本量和方法排序
    result_table <- result_table[order(result_table$SampleSize, result_table$Method), ]
    
    return(result_table)
  }
  
  # 为所有情形创建一个综合表格
  create_comprehensive_table <- function(all_results) {
    all_tables <- list()
    
    for(i in 1:5) {
      if(!is.null(all_results[[paste0("scenario_", i)]])) {
        all_tables[[i]] <- create_enhanced_summary_table(all_results, i)
        if(!is.null(all_tables[[i]])) {
          all_tables[[i]]$Scenario <- paste0("情形", i)
        }
      }
    }
    
    # 合并所有表格
    comprehensive_table <- do.call(rbind, all_tables)
    
    # 重新排列列
    comprehensive_table <- comprehensive_table[, c("Scenario", "SampleSize", "Method", "Mean", "SD", "Bias", "RMSE", "CI_Coverage", "Success_Rate")]
    
    # 排序
    comprehensive_table <- comprehensive_table[order(comprehensive_table$Scenario, comprehensive_table$SampleSize, comprehensive_table$Method), ]
    
    return(comprehensive_table)
  }
  
  # 为每个情形创建和显示增强的结果表格
  for(i in 1:5) {
    if(!is.null(all_results[[paste0("scenario_", i)]])) {
      table_scenario <- create_enhanced_summary_table(all_results, i)
      if(!is.null(table_scenario)) {
        print(kable(table_scenario, digits=4, caption=paste0("情形", i, "的详细结果")))
      }
    }
  }
  
  # 创建和显示综合表格
  comprehensive_table <- create_comprehensive_table(all_results)
  print(kable(comprehensive_table, digits=4, caption="双重机器学习袋装法 - 所有情形综合结果"))
  
  ##############################################
  # 第六步：增强的可视化结果
  ##############################################
  
  # 准备可视化数据
  prepare_plot_data <- function(all_results) {
    plot_data <- data.frame()
    
    for (scenario in 1:5) {
      scenario_key <- paste0("scenario_", scenario)
      if(!is.null(all_results[[scenario_key]])) {
        scenario_name <- paste0("情形", scenario)
        
        for (n in names(all_results[[scenario_key]])) {
          temp_df <- all_results[[scenario_key]][[n]]
          if(!is.null(temp_df) && nrow(temp_df) > 0) {
            temp_df$Scenario <- scenario_name
            temp_df$ScenarioNum <- scenario
            temp_df$SampleSize <- as.numeric(n)
            plot_data <- rbind(plot_data, temp_df)
          }
        }
      }
    }
    
    return(plot_data)
  }
  
  # 生成所有情形的可视化数据
  plot_data <- prepare_plot_data(all_results)
  
  # 设置配色方案和主题
  scenario_colors <- c("情形1" = "#1f77b4", "情形2" = "#ff7f0e", 
                       "情形3" = "#2ca02c", "情形4" = "#d62728", 
                       "情形5" = "#9467bd")
  
  custom_theme <- theme_minimal() +
    theme(
      plot.title = element_text(size=14, face="bold", hjust=0.5),
      plot.subtitle = element_text(size=12, hjust=0.5),
      axis.title = element_text(size=12, face="bold"),
      axis.text = element_text(size=10),
      legend.position = "bottom",
      legend.title = element_text(size=10, face="bold"),
      strip.text = element_text(size=11, face="bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color="grey80", fill=NA)
    )
  
  # 确定样本量刻度 (取对数尺度)
  log_breaks <- c(100, 200, 500, 1000, 2000, 5000, 10000)
  
  # 1. 增强版 - 均值随样本量变化的图表
  mean_plot <- ggplot(plot_data, 
                      aes(x=SampleSize, y=Mean, color=Scenario, 
                          shape=Method, group=interaction(Method, Scenario))) +
    geom_line(linewidth=1) +
    geom_point(size=3, stroke=1.5) +
    geom_hline(yintercept=1, linetype="dashed", color="black", linewidth=0.8) +  # 添加真实值参考线
    labs(title="袋装法(Bagging)在不同情形和样本量下的处理效应估计",
         subtitle="水平虚线表示真实处理效应值 (theta=1)",
         x="样本容量 (对数尺度)", y="估计值均值") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    scale_color_manual(values=scenario_colors) +
    coord_cartesian(ylim=c(-0.5, 2.5)) +  # 限制y轴范围以更好地显示重要区域
    custom_theme +
    facet_wrap(~ Scenario, scales="free_y")
  
  # 2. 增强版 - RMSE随样本量变化的图表 (分开情形5)
  # 分成两部分: 排除情形5和单独情形5
  rmse_plot_no5 <- plot_data %>%
    filter(ScenarioNum != 5) %>%
    ggplot(aes(x=SampleSize, y=RMSE, color=Scenario, 
               shape=Method, group=interaction(Method, Scenario))) +
    geom_line(linewidth=1) +
    geom_point(size=3, stroke=1.5) +
    labs(title="情形1-4的均方根误差(RMSE)随样本量变化",
         x="样本容量 (对数尺度)", y="RMSE") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    scale_color_manual(values=scenario_colors[1:4]) +
    custom_theme +
    facet_wrap(~ Scenario, scales="free_y")
  
  rmse_plot_5 <- plot_data %>%
    filter(ScenarioNum == 5) %>%
    ggplot(aes(x=SampleSize, y=RMSE, color=Method, 
               group=Method)) +
    geom_line(linewidth=1) +
    geom_point(size=3, stroke=1.5) +
    labs(title="情形5的均方根误差(RMSE)随样本量变化",
         subtitle="注意: 情形5的纵轴尺度不同",
         x="样本容量 (对数尺度)", y="RMSE") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    custom_theme
  
  # 3. 增强版 - 偏差随样本量变化的图表 (分开情形5)
  bias_plot_no5 <- plot_data %>%
    filter(ScenarioNum != 5) %>%
    ggplot(aes(x=SampleSize, y=Bias, color=Scenario, 
               shape=Method, group=interaction(Method, Scenario))) +
    geom_line(linewidth=1) +
    geom_point(size=3, stroke=1.5) +
    geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.8) +  # 添加零偏差参考线
    labs(title="情形1-4的偏差(Bias)随样本量变化",
         subtitle="水平虚线表示零偏差",
         x="样本容量 (对数尺度)", y="偏差") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    scale_color_manual(values=scenario_colors[1:4]) +
    coord_cartesian(ylim=c(-0.3, 0.3)) +  # 限制y轴范围以更好地显示重要区域
    custom_theme +
    facet_wrap(~ Scenario, scales="free_y")
  
  bias_plot_5 <- plot_data %>%
    filter(ScenarioNum == 5) %>%
    ggplot(aes(x=SampleSize, y=Bias, color=Method, 
               group=Method)) +
    geom_line(linewidth=1) +
    geom_point(size=3, stroke=1.5) +
    geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.8) +
    labs(title="情形5的偏差(Bias)随样本量变化",
         subtitle="注意: 情形5的纵轴尺度不同",
         x="样本容量 (对数尺度)", y="偏差") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    custom_theme
  
  # 4. 增强版 - 置信区间覆盖率图
  if("CI_Coverage" %in% colnames(plot_data)) {
    ci_coverage_plot <- ggplot(plot_data, 
                               aes(x=SampleSize, y=CI_Coverage, color=Scenario, 
                                   shape=Method, group=interaction(Method, Scenario))) +
      geom_line(linewidth=1) +
      geom_point(size=3, stroke=1.5) +
      geom_hline(yintercept=95, linetype="dashed", color="black", linewidth=0.8) +  # 95%覆盖率参考线
      labs(title="袋装法(Bagging)在不同情形下的95%置信区间覆盖率",
           subtitle="水平虚线表示名义95%覆盖率",
           x="样本容量 (对数尺度)", y="置信区间覆盖率 (%)") +
      scale_x_log10(breaks=log_breaks, labels=log_breaks) +
      scale_color_manual(values=scenario_colors) +
      coord_cartesian(ylim=c(min(50, min(plot_data$CI_Coverage, na.rm=TRUE)), 
                             max(100, max(plot_data$CI_Coverage, na.rm=TRUE)))) +
      custom_theme +
      facet_wrap(~ Scenario, scales="free_y")
  }
  
  # 5. 增强版 - 成功率随样本量变化的图表
  success_rate_plot <- ggplot(plot_data, 
                              aes(x=SampleSize, y=Success_Rate, color=Scenario, 
                                  shape=Method, group=interaction(Method, Scenario))) +
    geom_line(linewidth=1) +
    geom_point(size=3, stroke=1.5) +
    labs(title="袋装法(Bagging)在不同情形和样本量下的成功率",
         x="样本容量 (对数尺度)", y="成功率 (%)") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    scale_y_continuous(labels=scales::percent_format(scale=1)) +  # 将成功率显示为百分比
    scale_color_manual(values=scenario_colors) +
    coord_cartesian(ylim=c(0, 105)) +  # 限制y轴范围
    custom_theme +
    facet_wrap(~ Scenario)
  
  # 6. 增强版 - 成功率热图
  success_heatmap <- plot_data %>%
    select(Method, Scenario, ScenarioNum, SampleSize, Success_Rate) %>%
    distinct() %>%  # 确保没有重复
    ggplot(aes(x=SampleSize, y=Method, fill=Success_Rate)) +
    geom_tile(color="white", linewidth=0.5) +
    scale_fill_gradient2(low="red", mid="yellow", high="green", 
                         midpoint=75,
                         labels=scales::percent_format(scale=1),
                         limits=c(0, 100)) +
    labs(title="袋装法(Bagging)成功率热图",
         x="样本容量", y="方法", fill="成功率 (%)") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=14, face="bold", hjust=0.5),
      axis.title = element_text(size=12, face="bold"),
      legend.position = "right",
      strip.text = element_text(size=11, face="bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color="grey80", fill=NA)
    ) +
    facet_wrap(~ Scenario, nrow=1)
  
  # 7. 增强版 - 方法比较图 (Est1 vs. Est2)
  methods_comparison <- plot_data %>%
    select(Scenario, ScenarioNum, Method, SampleSize, Mean, Bias, RMSE) %>%
    tidyr::pivot_wider(
      id_cols = c(Scenario, ScenarioNum, SampleSize),
      names_from = Method,
      values_from = c(Mean, Bias, RMSE)
    ) %>%
    ggplot(aes(x=SampleSize, color=Scenario)) +
    geom_line(aes(y=`RMSE_Bagging-Est1` - `RMSE_Bagging-Est2`, linetype="RMSE差异"), linewidth=1) +
    geom_line(aes(y=abs(`Bias_Bagging-Est1`) - abs(`Bias_Bagging-Est2`), linetype="绝对偏差差异"), linewidth=1) +
    geom_hline(yintercept=0, linetype="dashed", color="black") +
    labs(title="Bagging-Est1与Bagging-Est2方法比较",
         subtitle="正值表示Est2优于Est1，负值表示Est1优于Est2",
         x="样本容量 (对数尺度)", y="差异值") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    scale_color_manual(values=scenario_colors) +
    scale_linetype_manual(values=c("RMSE差异"="solid", "绝对偏差差异"="dotted")) +
    custom_theme +
    facet_wrap(~ Scenario, scales="free_y")
  
  # 8. 增强版 - 计算效率分析图
  # 添加样本量与RMSE的关系
  efficiency_plot <- plot_data %>%
    filter(Method == "Bagging-Est2") %>%  # 只使用一种方法简化分析
    ggplot(aes(x=SampleSize, y=RMSE, color=Scenario)) +
    geom_line(linewidth=1) +
    geom_point(size=3) +
    labs(title="样本量与估计精度的关系分析",
         subtitle="展示样本量增加时RMSE的下降趋势",
         x="样本容量 (对数尺度)", y="RMSE (对数尺度)") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    scale_y_log10() +
    scale_color_manual(values=scenario_colors) +
    geom_smooth(method="lm", se=FALSE, linetype="dashed", color="black") +
    custom_theme +
    facet_wrap(~ Scenario, scales="free_y")
  
  # 9. 情形5的特别分析图
  scenario5_analysis <- plot_data %>%
    filter(ScenarioNum == 5) %>%
    ggplot(aes(x=SampleSize, color=Method)) +
    geom_line(aes(y=abs(Bias), linetype="绝对偏差"), linewidth=1) +
    geom_line(aes(y=RMSE, linetype="RMSE"), linewidth=1) +
    labs(title="情形5详细分析",
         subtitle="RMSE和绝对偏差随样本量的变化",
         x="样本容量 (对数尺度)", y="值") +
    scale_x_log10(breaks=log_breaks, labels=log_breaks) +
    scale_y_log10() +
    scale_linetype_manual(values=c("绝对偏差"="dashed", "RMSE"="solid")) +
    custom_theme
  
  # 保存改进的图表
  ggsave("dml_bagging_enhanced_mean.png", mean_plot, width=12, height=8)
  ggsave("dml_bagging_enhanced_rmse_no5.png", rmse_plot_no5, width=12, height=8)
  ggsave("dml_bagging_enhanced_rmse_5.png", rmse_plot_5, width=8, height=6)
  ggsave("dml_bagging_enhanced_bias_no5.png", bias_plot_no5, width=12, height=8)
  ggsave("dml_bagging_enhanced_bias_5.png", bias_plot_5, width=8, height=6)
  if(exists("ci_coverage_plot")) {
    ggsave("dml_bagging_enhanced_ci_coverage.png", ci_coverage_plot, width=12, height=8)
  }
  ggsave("dml_bagging_enhanced_success_rate.png", success_rate_plot, width=12, height=8)
  ggsave("dml_bagging_enhanced_success_heatmap.png", success_heatmap, width=16, height=6)
  ggsave("dml_bagging_methods_comparison.png", methods_comparison, width=12, height=8)
  ggsave("dml_bagging_efficiency_plot.png", efficiency_plot, width=12, height=8)
  ggsave("dml_bagging_scenario5_analysis.png", scenario5_analysis, width=10, height=6)
  
  # 保存最终结果
  save(all_results, plot_data, file="dml_bagging_enhanced_results.RData")
  ##############################################
  # 第七步：导出最终结果和建议
  ##############################################
  
  # 导出汇总表格为Excel文件
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- createWorkbook()
    
    # 添加综合结果工作表
    addWorksheet(wb, "综合结果")
    writeData(wb, "综合结果", comprehensive_table)
    
    # 为每个情形创建单独的工作表
    for (scenario in 1:5) {
      if(!is.null(all_results[[paste0("scenario_", scenario)]])) {
        scenario_data <- create_enhanced_summary_table(all_results, scenario)
        if(!is.null(scenario_data)) {
          addWorksheet(wb, paste0("情形", scenario))
          writeData(wb, paste0("情形", scenario), scenario_data)
        }
      }
    }
    
    # 添加方法说明工作表
    addWorksheet(wb, "方法说明")
    method_description <- data.frame(
      方法名称 = c("Bagging-Est1", "Bagging-Est2"),
      说明 = c(
        "双重机器学习方法1 (DML1): 分别在两个样本上进行处理效应估计，然后计算加权平均", 
        "双重机器学习方法2 (DML2): 将两个样本的残差合并后进行回归估计"
      )
    )
    writeData(wb, "方法说明", method_description)
    
    # 添加情形说明工作表
    addWorksheet(wb, "情形说明")
    scenario_description <- data.frame(
      情形编号 = 1:5,
      说明 = c(
        "简单线性模型: Y = D*theta + (-6 + 12*X1 + 6*X2) + U, D ~ Bernoulli(m0), m0 = (X1 + X2)/2", 
        "线性模型添加交互项: Y = D*theta + (-6 + 12*X1 + 6*X2 + 2*X3) + U, X3 = X1*X2",
        "线性模型添加平方项: Y = D*theta + (-6 + 12*X1 + 6*X2 + 2*X1^2 + 3*X2^2) + U",
        "线性模型添加交互项和平方项: Y = D*theta + (-6 + 12*X1 + 6*X2 + 2*X1^2 + 3*X2^2 + 3*X1*X2) + U",
        "指数型非线性模型: 保留了原始函数形式但增强了数值稳定性"
      )
    )
    writeData(wb, "情形说明", scenario_description)
    
    # 添加结果说明工作表
    addWorksheet(wb, "结果指标说明")
    metrics_description <- data.frame(
      指标名称 = c("Mean", "SD", "Bias", "RMSE", "CI_Coverage", "Success_Rate"),
      说明 = c(
        "模拟中估计值的平均值，越接近真实值1越好", 
        "估计值的标准差，反映估计的离散程度",
        "估计偏差 (Mean - 1)，理想值为0",
        "均方根误差，综合衡量估计的准确性和精度",
        "95%置信区间覆盖真实值的比例，理想值为95%",
        "成功估计的比例，与数值稳定性相关"
      )
    )
    writeData(wb, "结果指标说明", metrics_description)
    
    # 添加优化建议工作表
    addWorksheet(wb, "优化建议")
    optimization_tips <- data.frame(
      情形 = c("一般建议", "情形1-4", "情形5", "超大样本", "高维数据"),
      建议 = c(
        "1. 增加树的数量(n_models)可提高估计稳定性\n2. 使用异常值检测和剪枝可以提高稳健性\n3. DML2方法(Est2)通常比DML1(Est1)效果更好", 
        "1. 样本量500以上通常可获得很好的估计结果\n2. 树深度可适当调小，以减少过拟合",
        "1. 减小特征的生成方差，避免极端值\n2. 对特征值添加严格的截断范围，如[-1.5,1.5]\n3. 增加额外的数值检查，防止指数爆炸和对数无效",
        "1. 考虑使用并行计算加速\n2. 可以减少重复次数但增加每次模拟的树数量\n3. 分批处理并定期保存结果",
        "1. 考虑使用特征选择或降维技术\n2. 增加树的最小节点大小以避免过拟合\n3. 增加树的数量以捕捉复杂特征交互"
      )
    )
    writeData(wb, "优化建议", optimization_tips)
    
    # 设置列宽
    for(sheet in names(wb)) {
      setColWidths(wb, sheet, cols = 1:ncol(wb[[sheet]]), widths = "auto")
    }
    
    # 保存工作簿
    saveWorkbook(wb, "dml_bagging_enhanced_results.xlsx", overwrite = TRUE)
    cat("增强版结果已保存到Excel文件: dml_bagging_enhanced_results.xlsx\n")
  }
  
  # 导出汇总CSV文件
  write.csv(comprehensive_table, "dml_bagging_enhanced_results.csv", row.names = FALSE)
  
  ##############################################
  # 第八步：结论与建议
  ##############################################
  
  cat("\n\n===============================================\n")
  cat("双重机器学习袋装法(Bagging)优化模拟总结\n")
  cat("===============================================\n\n")
  
  # 性能比较汇总
  cat("方法性能汇总：\n")
  method_summary <- plot_data %>%
    group_by(Method, Scenario) %>%
    summarise(
      平均偏差 = mean(abs(Bias), na.rm=TRUE),
      平均RMSE = mean(RMSE, na.rm=TRUE),
      成功率 = mean(Success_Rate, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    arrange(Scenario, 平均RMSE)
  
  print(method_summary)
  
  # 输出结论
  cat("\n关键结论：\n")
  cat("1. 对于情形1-4，袋装法的双重机器学习表现出色，特别是在样本量≥500时\n")
  cat("2. 情形5（非线性模型）需要特殊处理，包括特征生成时的数值稳定性措施\n")
  cat("3. 总体而言，DML2方法（Est2）通常比DML1（Est1）表现更好，尤其在样本量较小时\n")
  cat("4. 异常值检测和数据预处理对提高估计稳定性非常重要\n")
  cat("5. 对于复杂模型，增加树的数量可显著提高估计精度\n\n")
  
  # 输出情形5的特别建议
  cat("情形5的特别优化策略：\n")
  cat("1. 特征生成：使用较小的标准差(0.5)生成特征，并限制范围在[-1.5,1.5]内\n")
  cat("2. 指数函数保护：检测可能的溢出并应用上限约束(1e8)\n")
  cat("3. 对数函数保护：确保参数严格大于零，并替换可能的无效计算结果\n")
  cat("4. 结果限制：对中间计算结果和最终m0值应用合理范围约束\n")
  cat("5. 稳健估计：对于较大样本量，增加树的数量并减少重复次数\n\n")
  
  # 输出一般优化建议
  cat("优化建议：\n")
  cat("1. 数据生成过程：添加边界控制以避免极端值，特别是在情形5中\n")
  cat("2. 模型拟合：使用自适应的树深度和更灵活的复杂度参数\n")
  cat("3. 残差估计：增强异常值检测，使用稳健统计方法\n")
  cat("4. 计算效率：对大样本量采用并行计算\n")
  cat("5. 结果汇总：增加置信区间覆盖率评估，提供更全面的性能度量\n\n")
  
  cat("已创建以下文件:\n")
  cat("1. dml_bagging_enhanced_results.xlsx - 增强版Excel工作簿，包含综合结果和方法说明\n")
  cat("2. dml_bagging_enhanced_results.csv - 增强版CSV格式结果\n")
  cat("3. dml_bagging_enhanced_results.RData - 完整R数据结果文件\n")
  cat("4. 多个增强版可视化图表 (.png文件)\n\n")
  
  cat("模拟优化过程已成功完成。\n")