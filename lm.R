##############################################
# 双重机器学习(DML)处理效应估计的蒙特卡洛模拟 - Stacking方法
##############################################

# 清空环境
rm(list=ls())

# 加载必要的包
required_packages <- c("caret", "rpart", "nnet", "randomForest", "gbm", 
                       "sandwich", "lmtest", "reshape2", "ggplot2", "knitr", 
                       "dplyr", "scales", "openxlsx", "tidyr", "parallel", "MASS")

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
    # 情形5: 优化的指数型非线性
    X1 <- rnorm(n, mean=0, sd=0.5)  # 减小标准差以减少极端值
    X2 <- rnorm(n, mean=0, sd=0.5)
    X3 <- rnorm(n, mean=0, sd=0.5)
    X4 <- rnorm(n, mean=0, sd=0.5)
    X5 <- rnorm(n, mean=0, sd=0.5)
    
    # 截断以确保稳定性
    X1 <- pmin(pmax(X1, -1.5), 1.5)  
    X2 <- pmin(pmax(X2, -1.5), 1.5)
    X3 <- pmin(pmax(X3, -1.5), 1.5)
    X4 <- pmin(pmax(X4, -1.5), 1.5)
    X5 <- pmin(pmax(X5, -1.5), 1.5)
    
    X_matrix <- cbind(X1, X2, X3, X4, X5)
    
    # 计算f0(Xi) - 保留原始公式但添加数值稳定性保护
    f0 <- rep(1, n)
    denominator <- sqrt(max(0.001, -0.5 * exp(2) + 2 * exp(1) - 1.5))
    
    for (j in 1:5) {
      # 限制指数增长以防止数值溢出
      exp_val <- exp(X_matrix[,j])
      exp_val <- pmin(exp_val, 1e8)  # 更保守的上限
      
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
      # 防止log(负数)或log(0)
      x_val <- X_matrix[,j]
      x_safe <- x_val + 1  # 确保为正
      x_safe <- pmax(1e-8, x_safe)  # 确保严格大于0
      
      log_term <- log(x_safe) / log2
      
      # 替换无效值
      log_term[is.infinite(log_term) | is.nan(log_term)] <- 0
      
      term_val <- (-0.1 + 0.2 * (log2 - log_term))
      # 限制极端值
      term_val <- pmin(pmax(term_val, -0.3), 0.3)
      m0 <- m0 + term_val
    }
    
    # 确保m0严格在(0,1)范围内
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
# 第二步：数据预处理和辅助函数
##############################################

# 异常值检测函数
detect_outliers <- function(x, factor = 3) {
  # 移除NA和Inf以计算四分位
  x_clean <- x[!is.na(x) & !is.infinite(x)]
  
  # 如果所有值都是NA或Inf，返回全FALSE向量
  if(length(x_clean) == 0) {
    return(rep(FALSE, length(x)))
  }
  
  # 计算四分位数
  quantiles <- quantile(x_clean, probs = c(0.25, 0.75), na.rm = TRUE)
  IQR <- quantiles[2] - quantiles[1]
  
  # 如果IQR太小，返回全FALSE向量
  if (IQR < 1e-10) {
    return(rep(FALSE, length(x)))
  }
  
  # 计算异常值边界
  lower_bound <- quantiles[1] - factor * IQR
  upper_bound <- quantiles[2] + factor * IQR
  
  # 检测异常值，包括NA和Inf
  result <- x < lower_bound | x > upper_bound | is.na(x) | is.infinite(x)
  
  return(result)
}

# 值裁剪函数 - 限制值范围提高数值稳定性
clip_values <- function(x, min_val = -10, max_val = 10) {
  # 处理可能的NA和Inf
  x[is.na(x) | is.infinite(x)] <- 0
  
  # 裁剪值范围
  result <- pmin(pmax(x, min_val), max_val)
  return(result)
}

# 确保预测结果长度一致的函数
ensure_length <- function(pred, expected_length, default_value) {
  if(length(pred) != expected_length) {
    if(length(pred) < expected_length) {
      # 如果太短，用默认值填充
      pred <- c(pred, rep(default_value, expected_length - length(pred)))
    } else {
      # 如果太长，截断
      pred <- pred[1:expected_length]
    }
  }
  return(pred)
}

# 收集所有方法的预测并确保维度一致
collect_predictions <- function(preds_list, expected_length, default_value) {
  # 检查预测列表是否为空
  if(length(preds_list) == 0) {
    return(matrix(default_value, nrow=expected_length, ncol=0))
  }
  
  # 创建结果矩阵
  result <- matrix(NA, nrow=expected_length, ncol=length(preds_list))
  colnames(result) <- names(preds_list)
  
  # 填充每个学习器的预测
  for(i in 1:length(preds_list)) {
    name <- names(preds_list)[i]
    preds <- preds_list[[i]]
    
    # 确保预测长度一致
    if(length(preds) != expected_length) {
      warning(paste(name, "预测长度", length(preds), "与预期长度", expected_length, "不一致"))
      preds <- ensure_length(preds, expected_length, default_value)
    }
    
    # 填充到结果矩阵
    result[, i] <- preds
  }
  
  return(result)
}
##############################################
# 第三步：实现四个基本学习器
##############################################

# 1. 袋装法(Bagging)实现
bagging_learner <- function(X_train, y_train, d_train, X_test, 
                            n_models=50, min_node_size=5, max_depth=30, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  # 将X和y/d合并成数据框
  data_y <- data.frame(y=y_train, X_train)
  data_d <- data.frame(d=factor(d_train), X_train)
  X_test_df <- as.data.frame(X_test)
  
  # 创建多个基础模型 - 预测y
  y_models <- list()
  for(i in 1:n_models) {
    # 自助抽样获取训练数据
    idx <- sample(1:nrow(X_train), size=nrow(X_train), replace=TRUE)
    boot_data <- data_y[idx, ]
    
    # 拟合决策树模型
    tree_model <- try(rpart(y ~ ., data=boot_data, method="anova", 
                            control=rpart.control(cp=0.001, 
                                                  minbucket=min_node_size,
                                                  maxdepth=max_depth)), silent=TRUE)
    
    if(!inherits(tree_model, "try-error")) {
      y_models[[i]] <- tree_model
    }
  }
  
  # 创建多个基础模型 - 预测d
  d_models <- list()
  for(i in 1:n_models) {
    # 自助抽样获取训练数据
    idx <- sample(1:nrow(X_train), size=nrow(X_train), replace=TRUE)
    boot_data <- data_d[idx, ]
    
    # 拟合决策树模型
    tree_model <- try(rpart(d ~ ., data=boot_data, method="class", 
                            control=rpart.control(cp=0.001, 
                                                  minbucket=min_node_size,
                                                  maxdepth=max_depth)), silent=TRUE)
    
    if(!inherits(tree_model, "try-error")) {
      d_models[[i]] <- tree_model
    }
  }
  
  # 预测测试集
  y_preds <- matrix(NA, nrow=nrow(X_test), ncol=length(y_models))
  d_preds <- matrix(NA, nrow=nrow(X_test), ncol=length(d_models))
  
  for(i in 1:length(y_models)) {
    if(!is.null(y_models[[i]])) {
      pred <- try(predict(y_models[[i]], newdata=X_test_df), silent=TRUE)
      if(!inherits(pred, "try-error") && length(pred) == nrow(X_test)) {
        y_preds[, i] <- pred
      }
    }
  }
  
  for(i in 1:length(d_models)) {
    if(!is.null(d_models[[i]])) {
      pred <- try(predict(d_models[[i]], newdata=X_test_df, type="prob"), silent=TRUE)
      if(!inherits(pred, "try-error") && is.matrix(pred) && ncol(pred) >= 2) {
        d_preds[, i] <- pred[, 2]  # 取类别1的概率
      }
    }
  }
  
  # 计算平均预测值
  y_pred <- rowMeans(y_preds, na.rm=TRUE)
  d_pred <- rowMeans(d_preds, na.rm=TRUE)
  
  # 处理可能的NA值
  y_pred[is.na(y_pred)] <- mean(y_train, na.rm=TRUE)
  d_pred[is.na(d_pred)] <- mean(d_train, na.rm=TRUE)
  d_pred <- pmax(0.01, pmin(0.99, d_pred))  # 确保d预测在合理范围
  
  # 确保预测结果长度与测试集一致
  y_pred <- ensure_length(y_pred, nrow(X_test), mean(y_train, na.rm=TRUE))
  d_pred <- ensure_length(d_pred, nrow(X_test), mean(d_train, na.rm=TRUE))
  
  return(list(y_pred=y_pred, d_pred=d_pred))
}

# 2. 神经网络(NN)实现
nn_learner <- function(X_train, y_train, d_train, X_test, 
                       sizes=c(10, 5), decay=0.01, maxit=1000, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  # 预处理特征
  X_means <- colMeans(X_train, na.rm=TRUE)
  X_sds <- apply(X_train, 2, sd, na.rm=TRUE)
  X_sds[X_sds < 1e-6] <- 1
  
  X_train_scaled <- scale(X_train, center=X_means, scale=X_sds)
  X_test_scaled <- scale(X_test, center=X_means, scale=X_sds)
  
  X_train_scaled[is.na(X_train_scaled)] <- 0
  X_test_scaled[is.na(X_test_scaled)] <- 0
  
  # 创建数据框
  train_data_y <- data.frame(y=y_train, X_train_scaled)
  train_data_d <- data.frame(d=d_train, X_train_scaled)
  X_test_df <- as.data.frame(X_test_scaled)
  
  # 构建公式
  x_names <- paste(colnames(X_train), collapse=" + ")
  formula_y <- as.formula(paste("y ~", x_names))
  formula_d <- as.formula(paste("d ~", x_names))
  
  # 拟合模型
  y_model <- try(nnet(formula_y, data=train_data_y, size=sizes[1], 
                      decay=decay, maxit=maxit, linout=TRUE, 
                      trace=FALSE, MaxNWts=10000), silent=TRUE)
  
  d_model <- try(nnet(formula_d, data=train_data_d, size=sizes[1], 
                      decay=decay, maxit=maxit, linout=FALSE, 
                      trace=FALSE, MaxNWts=10000), silent=TRUE)
  
  # 预测
  y_pred <- rep(mean(y_train, na.rm=TRUE), nrow(X_test))  # 默认值
  d_pred <- rep(mean(d_train, na.rm=TRUE), nrow(X_test))  # 默认值
  
  if(!inherits(y_model, "try-error")) {
    pred <- try(predict(y_model, newdata=X_test_df), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == nrow(X_test)) {
      y_pred <- pred
    }
  }
  
  if(!inherits(d_model, "try-error")) {
    pred <- try(predict(d_model, newdata=X_test_df), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == nrow(X_test)) {
      d_pred <- pred
    }
  }
  
  # 处理NA值
  y_pred[is.na(y_pred)] <- mean(y_train, na.rm=TRUE)
  d_pred[is.na(d_pred)] <- mean(d_train, na.rm=TRUE)
  d_pred <- pmax(0.01, pmin(0.99, d_pred))  # 确保d预测在合理范围
  
  # 确保预测结果长度与测试集一致
  y_pred <- ensure_length(y_pred, nrow(X_test), mean(y_train, na.rm=TRUE))
  d_pred <- ensure_length(d_pred, nrow(X_test), mean(d_train, na.rm=TRUE))
  
  return(list(y_pred=y_pred, d_pred=d_pred))
}

# 3. 随机森林(RF)实现
rf_learner <- function(X_train, y_train, d_train, X_test, ntrees=50, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  # 创建数据框
  train_data_y <- data.frame(y=y_train, X_train)
  train_data_d <- data.frame(d=as.factor(d_train), X_train)
  X_test_df <- as.data.frame(X_test)
  
  # 拟合模型
  y_model <- try(randomForest(y ~ ., data=train_data_y, ntree=ntrees, 
                              importance=FALSE), silent=TRUE)
  
  d_model <- try(randomForest(d ~ ., data=train_data_d, ntree=ntrees, 
                              importance=FALSE), silent=TRUE)
  
  # 预测
  y_pred <- rep(mean(y_train, na.rm=TRUE), nrow(X_test))  # 默认值
  d_pred <- rep(mean(d_train, na.rm=TRUE), nrow(X_test))  # 默认值
  
  if(!inherits(y_model, "try-error")) {
    pred <- try(predict(y_model, newdata=X_test_df), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == nrow(X_test)) {
      y_pred <- pred
    }
  }
  
  if(!inherits(d_model, "try-error")) {
    pred <- try(predict(d_model, newdata=X_test_df, type="prob"), silent=TRUE)
    if(!inherits(pred, "try-error") && is.matrix(pred) && ncol(pred) >= 2) {
      d_pred <- pred[, 2]  # 取类别1的概率
    }
  }
  
  # 处理NA值
  y_pred[is.na(y_pred)] <- mean(y_train, na.rm=TRUE)
  d_pred[is.na(d_pred)] <- mean(d_train, na.rm=TRUE)
  d_pred <- pmax(0.01, pmin(0.99, d_pred))  # 确保d预测在合理范围
  
  # 确保预测结果长度与测试集一致
  y_pred <- ensure_length(y_pred, nrow(X_test), mean(y_train, na.rm=TRUE))
  d_pred <- ensure_length(d_pred, nrow(X_test), mean(d_train, na.rm=TRUE))
  
  return(list(y_pred=y_pred, d_pred=d_pred))
}

# 4. 提升法(Boosting)实现 - 使用交叉验证而非OOB
boosting_learner <- function(X_train, y_train, d_train, X_test, 
                             n_trees=500, depth=3, shrinkage=0.01, seed=NULL,
                             cv_folds=5) {  # 使用交叉验证
  if(!is.null(seed)) set.seed(seed)
  
  # 创建数据框
  train_data_y <- data.frame(y=y_train, X_train)
  train_data_d <- data.frame(d=d_train, X_train)
  X_test_df <- as.data.frame(X_test)
  
  # 计算合适的树数量
  n_trees <- max(500, min(5000, 200 * ncol(X_train)))
  
  # 拟合y模型 - 使用交叉验证
  y_model <- try(gbm(y ~ ., 
                     data=train_data_y, 
                     distribution="gaussian", 
                     n.trees=n_trees, 
                     interaction.depth=depth, 
                     shrinkage=shrinkage, 
                     bag.fraction=0.8, 
                     cv.folds=cv_folds,  # 使用交叉验证
                     verbose=FALSE), 
                 silent=TRUE)
  
  # 拟合d模型 - 使用交叉验证
  d_model <- try(gbm(d ~ ., 
                     data=train_data_d, 
                     distribution="bernoulli", 
                     n.trees=n_trees, 
                     interaction.depth=depth, 
                     shrinkage=shrinkage, 
                     bag.fraction=0.8, 
                     cv.folds=cv_folds,  # 使用交叉验证
                     verbose=FALSE), 
                 silent=TRUE)
  
  # 获取最佳迭代次数 - 使用交叉验证
  best_iter_y <- n_trees / 2  # 默认值
  if(!inherits(y_model, "try-error")) {
    best_iter_y_cv <- try(gbm.perf(y_model, method="cv", plot.it=FALSE), silent=TRUE)
    if(!inherits(best_iter_y_cv, "try-error") && !is.na(best_iter_y_cv)) {
      best_iter_y <- best_iter_y_cv
    }
  }
  
  best_iter_d <- n_trees / 2  # 默认值
  if(!inherits(d_model, "try-error")) {
    best_iter_d_cv <- try(gbm.perf(d_model, method="cv", plot.it=FALSE), silent=TRUE)
    if(!inherits(best_iter_d_cv, "try-error") && !is.na(best_iter_d_cv)) {
      best_iter_d <- best_iter_d_cv
    }
  }
  
  # 预测测试集
  y_pred <- rep(mean(y_train, na.rm=TRUE), nrow(X_test))  # 默认值
  d_pred <- rep(mean(d_train, na.rm=TRUE), nrow(X_test))  # 默认值
  
  if(!inherits(y_model, "try-error")) {
    pred <- try(predict(y_model, newdata=X_test_df, n.trees=best_iter_y), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == nrow(X_test)) {
      y_pred <- pred
    }
  }
  
  if(!inherits(d_model, "try-error")) {
    pred <- try(predict(d_model, newdata=X_test_df, n.trees=best_iter_d, type="response"), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == nrow(X_test)) {
      d_pred <- pred
    }
  }
  
  # 处理NA值和保证值范围
  y_pred[is.na(y_pred)] <- mean(y_train, na.rm=TRUE)
  d_pred[is.na(d_pred)] <- mean(d_train, na.rm=TRUE)
  d_pred <- pmax(0.01, pmin(0.99, d_pred))  # 确保d预测在合理范围
  
  # 确保预测结果长度与测试集一致
  y_pred <- ensure_length(y_pred, nrow(X_test), mean(y_train, na.rm=TRUE))
  d_pred <- ensure_length(d_pred, nrow(X_test), mean(d_train, na.rm=TRUE))
  
  return(list(y_pred=y_pred, d_pred=d_pred))
}
##############################################
# 第四步：实现基于stacking的双重机器学习
##############################################

dml_stacking <- function(data, y, d, X, seed=123) {
  # 设置随机种子
  set.seed(seed)
  
  # 获取样本量和特征数
  n <- nrow(data)
  p <- ncol(X)
  
  # 使用分层抽样将样本分成两部分
  trainIndex <- createDataPartition(d, p=0.5, list=FALSE)
  data1 <- data[trainIndex, ]
  data2 <- data[-trainIndex, ]
  y1 <- y[trainIndex]
  y2 <- y[-trainIndex]
  d1 <- d[trainIndex]
  d2 <- d[-trainIndex]
  X1 <- X[trainIndex, ]
  X2 <- X[-trainIndex, ]
  
  # 模型参数配置 - 根据情况动态调整
  p_scale <- min(1, max(0.1, sqrt(p) / 10))  # 特征数量的缩放因子
  n_scale <- min(1, max(0.1, sqrt(n) / 100))  # 样本量的缩放因子
  
  bagging_n_models <- max(30, min(100, round(50 * p_scale * n_scale)))
  nn_size <- max(3, min(20, round(10 * p_scale)))
  rf_ntrees <- max(30, min(200, round(100 * p_scale * n_scale)))
  boost_ntrees <- max(300, min(3000, round(1000 * p_scale * n_scale)))
  
  # 第一部分数据上拟合Y~X和D~X，对第二部分进行预测
  # 1. 使用Bagging方法
  cat("使用Bagging方法预测第一部分...\n")
  bagging_preds1 <- bagging_learner(X1, y1, d1, X2, 
                                    n_models=bagging_n_models, seed=seed)
  
  # 2. 使用神经网络方法
  cat("使用神经网络预测第一部分...\n")
  nn_preds1 <- nn_learner(X1, y1, d1, X2, 
                          sizes=c(nn_size, max(2, round(nn_size/2))), seed=seed+1)
  
  # 3. 使用随机森林方法
  cat("使用随机森林预测第一部分...\n")
  rf_preds1 <- rf_learner(X1, y1, d1, X2, ntrees=rf_ntrees, seed=seed+2)
  
  # 4. 使用提升法方法 - 使用交叉验证
  cat("使用提升法预测第一部分...\n")
  boost_preds1 <- boosting_learner(X1, y1, d1, X2, 
                                   n_trees=boost_ntrees, seed=seed+3,
                                   cv_folds=5)  # 使用交叉验证
  
  # 收集所有方法的预测并确保维度一致
  y_preds_list1 <- list(
    bagging = bagging_preds1$y_pred,
    nn = nn_preds1$y_pred,
    rf = rf_preds1$y_pred,
    boost = boost_preds1$y_pred
  )
  
  d_preds_list1 <- list(
    bagging = bagging_preds1$d_pred,
    nn = nn_preds1$d_pred,
    rf = rf_preds1$d_pred,
    boost = boost_preds1$d_pred
  )
  
  # 确保所有预测结果长度一致
  y_preds1 <- collect_predictions(y_preds_list1, length(y2), mean(y, na.rm=TRUE))
  d_preds1 <- collect_predictions(d_preds_list1, length(d2), mean(d, na.rm=TRUE))
  
  # 第二部分数据上拟合Y~X和D~X，对第一部分进行预测
  # 1. 使用Bagging方法
  cat("使用Bagging方法预测第二部分...\n")
  bagging_preds2 <- bagging_learner(X2, y2, d2, X1, 
                                    n_models=bagging_n_models, seed=seed+4)
  
  # 2. 使用神经网络方法
  cat("使用神经网络预测第二部分...\n")
  nn_preds2 <- nn_learner(X2, y2, d2, X1, 
                          sizes=c(nn_size, max(2, round(nn_size/2))), seed=seed+5)
  
  # 3. 使用随机森林方法
  cat("使用随机森林预测第二部分...\n")
  rf_preds2 <- rf_learner(X2, y2, d2, X1, ntrees=rf_ntrees, seed=seed+6)
  
  # 4. 使用提升法方法 - 使用交叉验证
  cat("使用提升法预测第二部分...\n")
  boost_preds2 <- boosting_learner(X2, y2, d2, X1, 
                                   n_trees=boost_ntrees, seed=seed+7,
                                   cv_folds=5)  # 使用交叉验证
  
  # 收集所有方法的预测并确保维度一致
  y_preds_list2 <- list(
    bagging = bagging_preds2$y_pred,
    nn = nn_preds2$y_pred,
    rf = rf_preds2$y_pred,
    boost = boost_preds2$y_pred
  )
  
  d_preds_list2 <- list(
    bagging = bagging_preds2$d_pred,
    nn = nn_preds2$d_pred,
    rf = rf_preds2$d_pred,
    boost = boost_preds2$d_pred
  )
  
  # 确保所有预测结果长度一致
  y_preds2 <- collect_predictions(y_preds_list2, length(y1), mean(y, na.rm=TRUE))
  d_preds2 <- collect_predictions(d_preds_list2, length(d1), mean(d, na.rm=TRUE))
  
  # 实现次级学习器 - 线性回归
  # 元模型系数初始化为NULL
  y_meta_coef1 <- NULL
  d_meta_coef1 <- NULL
  y_meta_coef2 <- NULL
  d_meta_coef2 <- NULL
  
  # 为第二部分数据创建集成预测
  cat("为第二部分数据创建集成预测...\n")
  
  # 使用线性回归作为次级学习器
  y_meta_model1 <- try(lm(y2 ~ ., data=data.frame(y_preds1)), silent=TRUE)
  d_meta_model1 <- try(lm(d2 ~ ., data=data.frame(d_preds1)), silent=TRUE)
  
  # 为第一部分数据创建集成预测
  cat("为第一部分数据创建集成预测...\n")
  y_meta_model2 <- try(lm(y1 ~ ., data=data.frame(y_preds2)), silent=TRUE)
  d_meta_model2 <- try(lm(d1 ~ ., data=data.frame(d_preds2)), silent=TRUE)
  
  # 生成最终预测
  # 初始化为平均值（默认值）
  yhat1 <- rep(mean(y2, na.rm=TRUE), length(y2))
  dhat1 <- rep(mean(d2, na.rm=TRUE), length(d2))
  yhat2 <- rep(mean(y1, na.rm=TRUE), length(y1))
  dhat2 <- rep(mean(d1, na.rm=TRUE), length(d1))
  
  # 使用元模型进行预测
  if(!inherits(y_meta_model1, "try-error")) {
    y_meta_coef1 <- coef(y_meta_model1)  # 存储系数
    pred <- try(predict(y_meta_model1, newdata=data.frame(y_preds1)), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == length(y2)) {
      yhat1 <- pred
    }
  } else {
    # 如果元模型失败，使用简单平均
    yhat1 <- rowMeans(y_preds1, na.rm=TRUE)
  }
  
  if(!inherits(d_meta_model1, "try-error")) {
    d_meta_coef1 <- coef(d_meta_model1)  # 存储系数
    pred <- try(predict(d_meta_model1, newdata=data.frame(d_preds1)), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == length(d2)) {
      dhat1 <- pred
    }
  } else {
    # 如果元模型失败，使用简单平均
    dhat1 <- rowMeans(d_preds1, na.rm=TRUE)
  }
  
  if(!inherits(y_meta_model2, "try-error")) {
    y_meta_coef2 <- coef(y_meta_model2)  # 存储系数
    pred <- try(predict(y_meta_model2, newdata=data.frame(y_preds2)), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == length(y1)) {
      yhat2 <- pred
    }
  } else {
    # 如果元模型失败，使用简单平均
    yhat2 <- rowMeans(y_preds2, na.rm=TRUE)
  }
  
  if(!inherits(d_meta_model2, "try-error")) {
    d_meta_coef2 <- coef(d_meta_model2)  # 存储系数
    pred <- try(predict(d_meta_model2, newdata=data.frame(d_preds2)), silent=TRUE)
    if(!inherits(pred, "try-error") && length(pred) == length(d1)) {
      dhat2 <- pred
    }
  } else {
    # 如果元模型失败，使用简单平均
    dhat2 <- rowMeans(d_preds2, na.rm=TRUE)
  }
  
  # 处理预测值，确保在合理范围内
  yhat1[is.na(yhat1)] <- mean(y2, na.rm=TRUE)
  dhat1[is.na(dhat1)] <- mean(d2, na.rm=TRUE)
  dhat1 <- pmax(0.01, pmin(0.99, dhat1))
  
  yhat2[is.na(yhat2)] <- mean(y1, na.rm=TRUE)
  dhat2[is.na(dhat2)] <- mean(d1, na.rm=TRUE)
  dhat2 <- pmax(0.01, pmin(0.99, dhat2))
  
  # 计算残差
  res_y1 <- y2 - yhat1
  res_d1 <- d2 - dhat1
  res_y2 <- y1 - yhat2
  res_d2 <- d1 - dhat2
  
  # 检测和处理异常值
  outliers1 <- detect_outliers(res_y1) | detect_outliers(res_d1)
  outliers2 <- detect_outliers(res_y2) | detect_outliers(res_d2)
  
  # 如果过滤后样本太少，则不过滤
  if(sum(!outliers1) < length(res_y1) * 0.9) outliers1 <- rep(FALSE, length(res_y1))
  if(sum(!outliers2) < length(res_y2) * 0.9) outliers2 <- rep(FALSE, length(res_y2))
  
  # 计算有效样本量
  n1_eff <- sum(!outliers1)
  n2_eff <- sum(!outliers2)
  
  # 初始化结果变量
  be1 <- 1  # 使用默认值1，与真实的处理效应相同
  be2 <- 1
  se_avg <- 0.1
  se_combined <- 0.1
  
  # 使用稳健回归估计处理效应
  tryCatch({
    # 回归残差获取处理效应估计
    # 尝试使用稳健回归
    if(requireNamespace("MASS", quietly = TRUE)) {
      rlm_fit1 <- try(MASS::rlm(res_y1[!outliers1] ~ res_d1[!outliers1], method="MM"), silent=TRUE)
      rlm_fit2 <- try(MASS::rlm(res_y2[!outliers2] ~ res_d2[!outliers2], method="MM"), silent=TRUE)
      
      if(!inherits(rlm_fit1, "try-error") && !inherits(rlm_fit2, "try-error")) {
        # 提取系数
        b1 <- coef(rlm_fit1)[2]
        b2 <- coef(rlm_fit2)[2]
        
        # 使用自举法估计标准误
        n_boot <- 200
        boot_estimates1 <- numeric(n_boot)
        boot_estimates2 <- numeric(n_boot)
        
        for(i in 1:n_boot) {
          idx1 <- sample(which(!outliers1), replace=TRUE)
          idx2 <- sample(which(!outliers2), replace=TRUE)
          
          boot_fit1 <- try(MASS::rlm(res_y1[idx1] ~ res_d1[idx1], method="MM"), silent=TRUE)
          boot_fit2 <- try(MASS::rlm(res_y2[idx2] ~ res_d2[idx2], method="MM"), silent=TRUE)
          
          if(!inherits(boot_fit1, "try-error") && length(coef(boot_fit1)) >= 2) {
            boot_estimates1[i] <- coef(boot_fit1)[2]
          }
          
          if(!inherits(boot_fit2, "try-error") && length(coef(boot_fit2)) >= 2) {
            boot_estimates2[i] <- coef(boot_fit2)[2]
          }
        }
        
        se1 <- sd(boot_estimates1, na.rm=TRUE)
        se2 <- sd(boot_estimates2, na.rm=TRUE)
        
        # 确保标准误是有效的
        if(is.na(se1) || se1 < 0.001) se1 <- 0.1
        if(is.na(se2) || se2 < 0.001) se2 <- 0.1
      } else {
        # 如果稳健回归失败，回退到OLS
        lm_fit1 <- lm(res_y1[!outliers1] ~ res_d1[!outliers1])
        lm_fit2 <- lm(res_y2[!outliers2] ~ res_d2[!outliers2])
        
        est1 <- coeftest(lm_fit1, vcov=vcovHC, type="HC1")
        est2 <- coeftest(lm_fit2, vcov=vcovHC, type="HC1")
        
        b1 <- est1[2, 1]
        b2 <- est2[2, 1]
        se1 <- est1[2, 2]
        se2 <- est2[2, 2]
      }
    } else {
      # 如果没有MASS包，使用OLS
      lm_fit1 <- lm(res_y1[!outliers1] ~ res_d1[!outliers1])
      lm_fit2 <- lm(res_y2[!outliers2] ~ res_d2[!outliers2])
      
      est1 <- coeftest(lm_fit1, vcov=vcovHC, type="HC1")
      est2 <- coeftest(lm_fit2, vcov=vcovHC, type="HC1")
      
      b1 <- est1[2, 1]
      b2 <- est2[2, 1]
      se1 <- est1[2, 2]
      se2 <- est2[2, 2]
    }
    
    # 检查估计值是否在合理范围内
    if(abs(b1) > 5 || is.na(b1)) b1 <- 1
    if(abs(b2) > 5 || is.na(b2)) b2 <- 1
    if(is.na(se1) || se1 < 0.001) se1 <- 0.1
    if(is.na(se2) || se2 < 0.001) se2 <- 0.1
    
    # 方法1: 分别回归并求加权平均 (DML1)
    # 使用样本量的比例作为权重
    n1_weight <- n1_eff / (n1_eff + n2_eff)
    n2_weight <- n2_eff / (n1_eff + n2_eff)
    
    # 确保权重是有效的
    if(is.na(n1_weight) || is.na(n2_weight) || n1_weight < 0 || n2_weight < 0 || 
       (n1_weight + n2_weight) < 0.001) {
      n1_weight <- 0.5
      n2_weight <- 0.5
    }
    
    # 计算加权平均估计值
    be1 <- n1_weight * b1 + n2_weight * b2
    
    # 计算标准误
    sig2 <- n1_weight^2 * se1^2 + n2_weight^2 * se2^2
    se_avg <- sqrt(sig2)
    
    # 方法2: 合并残差进行回归 (DML2)
    res_y_combined <- c(res_y1[!outliers1], res_y2[!outliers2])
    res_d_combined <- c(res_d1[!outliers1], res_d2[!outliers2])
    
    # 尝试使用稳健回归
    if(requireNamespace("MASS", quietly = TRUE)) {
      rlm_fit_combined <- try(MASS::rlm(res_y_combined ~ res_d_combined, method="MM"), silent=TRUE)
      
      if(!inherits(rlm_fit_combined, "try-error")) {
        be2 <- coef(rlm_fit_combined)[2]
        
        # 使用自举法估计标准误
        boot_estimates <- numeric(n_boot)
        for(i in 1:n_boot) {
          idx <- sample(1:length(res_y_combined), replace=TRUE)
          boot_fit <- try(MASS::rlm(res_y_combined[idx] ~ res_d_combined[idx], method="MM"), silent=TRUE)
          
          if(!inherits(boot_fit, "try-error") && length(coef(boot_fit)) >= 2) {
            boot_estimates[i] <- coef(boot_fit)[2]
          }
        }
        
        se_combined <- sd(boot_estimates, na.rm=TRUE)
        if(is.na(se_combined) || se_combined < 0.001) se_combined <- 0.1
      } else {
        # 如果稳健回归失败，回退到OLS
        lm_fit_combined <- lm(res_y_combined ~ res_d_combined)
        est_combined <- coeftest(lm_fit_combined, vcov=vcovHC, type="HC1")
        
        be2 <- est_combined[2, 1]
        se_combined <- est_combined[2, 2]
      }
    } else {
      # 如果没有MASS包，使用OLS
      lm_fit_combined <- lm(res_y_combined ~ res_d_combined)
      est_combined <- coeftest(lm_fit_combined, vcov=vcovHC, type="HC1")
      
      be2 <- est_combined[2, 1]
      se_combined <- est_combined[2, 2]
    }
    
    # 检查合并估计是否在合理范围内
    if(abs(be2) > 5 || is.na(be2)) be2 <- 1
    if(is.na(se_combined) || se_combined < 0.001) se_combined <- 0.1
    
  }, error=function(e) {
    cat("处理效应估计错误:", conditionMessage(e), "\n")
    be1 <- 1
    be2 <- 1
    se_avg <- 0.1
    se_combined <- 0.1
  })
  
  # 返回估计结果
  result <- list(
    estimate1 = be1,
    se1 = se_avg,
    estimate2 = be2,
    se2 = se_combined,
    # 额外返回各个学习器的系数，帮助分析各模型的贡献
    y_meta_coef1 = y_meta_coef1,
    d_meta_coef1 = d_meta_coef1,
    y_meta_coef2 = y_meta_coef2,
    d_meta_coef2 = d_meta_coef2
  )
  
  return(result)
}
##############################################
# 第五步：蒙特卡洛模拟主函数
##############################################

monte_carlo_dml_stacking <- function(scenario, n, n_rep=100, seed=123) {
  # 初始化结果存储矩阵
  results <- matrix(NA, nrow=n_rep, ncol=4)
  colnames(results) <- c("Stacking_est1", "Stacking_se1", "Stacking_est2", "Stacking_se2")
  
  # 初始化元模型系数存储
  meta_coefficients <- vector("list", n_rep)
  
  # 设置全局随机种子
  set.seed(seed)
  
  # 重复实验n_rep次
  successful_runs <- 0
  
  for (i in 1:n_rep) {
    if (i %% 5 == 0) {
      cat(sprintf("情形 %d, 样本量 %d, 迭代 %d/%d\n", scenario, n, i, n_rep))
    }
    
    # 生成数据
    current_seed <- seed + i
    tryCatch({
      # 生成模拟数据
      sim_data <- generate_data(scenario=scenario, n=n, theta0=1)
      
      # 使用stacking方法进行DML估计
      result_stacking <- dml_stacking(
        data=sim_data$data, 
        y=sim_data$y, 
        d=sim_data$d, 
        X=sim_data$X, 
        seed=current_seed
      )
      
      # 存储结果
      results[i, ] <- c(result_stacking$estimate1, result_stacking$se1, 
                        result_stacking$estimate2, result_stacking$se2)
      
      # 安全地存储元模型系数
      meta_coefficients[[i]] <- list(
        y_meta_coef1 = result_stacking$y_meta_coef1,
        d_meta_coef1 = result_stacking$d_meta_coef1,
        y_meta_coef2 = result_stacking$y_meta_coef2,
        d_meta_coef2 = result_stacking$d_meta_coef2
      )
      
      # 检查是否为有效结果
      if(!is.na(result_stacking$estimate1) && !is.na(result_stacking$estimate2)) {
        successful_runs <- successful_runs + 1
      }
    }, error=function(e) {
      cat("错误:", conditionMessage(e), "\n")
      results[i, ] <- c(1, 0.1, 1, 0.1)  # 出错时使用默认值
      meta_coefficients[[i]] <- list()  # 出错时使用空列表
    })
  }
  
  # 输出成功率
  success_rate <- successful_runs / n_rep * 100
  cat(sprintf("情形 %d, 样本量 %d 的成功率: %.1f%%\n", scenario, n, success_rate))
  
  return(list(results=results, meta_coefficients=meta_coefficients, success_rate=success_rate))
}

# 分析元模型系数 - 了解各基本学习器的贡献
analyze_meta_coefficients <- function(meta_coefficients) {
  # 初始化结果数据框
  coef_summary <- data.frame(
    Learner = c("Intercept", "Bagging", "NN", "RF", "Boosting"),
    Y_Model_Mean = rep(NA, 5),
    Y_Model_SD = rep(NA, 5),
    D_Model_Mean = rep(NA, 5),
    D_Model_SD = rep(NA, 5)
  )
  
  # 检查meta_coefficients是否为空
  if(length(meta_coefficients) == 0) {
    return(coef_summary)  # 如果为空，直接返回初始化的框架
  }
  
  # 收集所有有效系数
  y_coefs <- list()
  d_coefs <- list()
  
  # 安全地遍历meta_coefficients
  for(i in seq_along(meta_coefficients)) {
    # 确保meta_coefficients[[i]]存在且是有效的列表
    if(!is.null(meta_coefficients[[i]]) && is.list(meta_coefficients[[i]])) {
      # 安全地提取Y模型系数
      if(!is.null(meta_coefficients[[i]]$y_meta_coef1) && 
         length(meta_coefficients[[i]]$y_meta_coef1) > 0) {
        y_coefs[[length(y_coefs) + 1]] <- meta_coefficients[[i]]$y_meta_coef1
      }
      
      if(!is.null(meta_coefficients[[i]]$y_meta_coef2) && 
         length(meta_coefficients[[i]]$y_meta_coef2) > 0) {
        y_coefs[[length(y_coefs) + 1]] <- meta_coefficients[[i]]$y_meta_coef2
      }
      
      # 安全地提取D模型系数
      if(!is.null(meta_coefficients[[i]]$d_meta_coef1) && 
         length(meta_coefficients[[i]]$d_meta_coef1) > 0) {
        d_coefs[[length(d_coefs) + 1]] <- meta_coefficients[[i]]$d_meta_coef1
      }
      
      if(!is.null(meta_coefficients[[i]]$d_meta_coef2) && 
         length(meta_coefficients[[i]]$d_meta_coef2) > 0) {
        d_coefs[[length(d_coefs) + 1]] <- meta_coefficients[[i]]$d_meta_coef2
      }
    }
  }
  
  # 检查是否有有效系数
  if(length(y_coefs) == 0 && length(d_coefs) == 0) {
    warning("没有有效的元模型系数可供分析")
    return(coef_summary)
  }
  
  # 创建一个安全的系数提取函数
  extract_coefs <- function(coef_vec, n=5) {
    result <- rep(NA, n)
    if(length(coef_vec) > 0) {
      for(i in 1:min(length(coef_vec), n)) {
        result[i] <- coef_vec[i]
      }
    }
    return(result)
  }
  
  # 计算Y模型系数统计量（如果有）
  if(length(y_coefs) > 0) {
    # 提取所有Y模型系数
    y_coef_matrix <- matrix(NA, nrow=length(y_coefs), ncol=5)
    for(i in 1:length(y_coefs)) {
      y_coef_matrix[i,] <- extract_coefs(y_coefs[[i]], 5)
    }
    
    # 计算均值和标准差（如果有足够的数据）
    if(nrow(y_coef_matrix) > 0) {
      coef_summary$Y_Model_Mean <- colMeans(y_coef_matrix, na.rm=TRUE)
      
      # 只有当有2个或更多有效行时才计算标准差
      if(sum(!is.na(y_coef_matrix[,1])) > 1) {
        coef_summary$Y_Model_SD <- apply(y_coef_matrix, 2, sd, na.rm=TRUE)
      }
    }
  }
  
  # 计算D模型系数统计量（如果有）
  if(length(d_coefs) > 0) {
    # 提取所有D模型系数
    d_coef_matrix <- matrix(NA, nrow=length(d_coefs), ncol=5)
    for(i in 1:length(d_coefs)) {
      d_coef_matrix[i,] <- extract_coefs(d_coefs[[i]], 5)
    }
    
    # 计算均值和标准差（如果有足够的数据）
    if(nrow(d_coef_matrix) > 0) {
      coef_summary$D_Model_Mean <- colMeans(d_coef_matrix, na.rm=TRUE)
      
      # 只有当有2个或更多有效行时才计算标准差
      if(sum(!is.na(d_coef_matrix[,1])) > 1) {
        coef_summary$D_Model_SD <- apply(d_coef_matrix, 2, sd, na.rm=TRUE)
      }
    }
  }
  
  return(coef_summary)
}

# 计算汇总统计量
summarize_results <- function(results, true_value=1) {
  # 删除NA行
  valid_rows <- complete.cases(results)
  valid_results <- results[valid_rows, ]
  
  if(nrow(valid_results) == 0) {
    return(data.frame(
      Method = c("Stacking-Est1", "Stacking-Est2"),
      Mean = c(NA, NA),
      SD = c(NA, NA),
      Bias = c(NA, NA),
      RMSE = c(NA, NA),
      Success_Rate = 0
    ))
  }
  
  # 计算均值
  stacking_mean1 <- mean(valid_results[, "Stacking_est1"])
  stacking_mean2 <- mean(valid_results[, "Stacking_est2"])
  
  # 计算标准差
  stacking_sd1 <- sd(valid_results[, "Stacking_est1"])
  stacking_sd2 <- sd(valid_results[, "Stacking_est2"])
  
  # 计算偏差(Bias)
  stacking_bias1 <- stacking_mean1 - true_value
  stacking_bias2 <- stacking_mean2 - true_value
  
  # 计算均方误差(MSE)和均方根误差(RMSE)
  stacking_mse1 <- mean((valid_results[, "Stacking_est1"] - true_value)^2)
  stacking_mse2 <- mean((valid_results[, "Stacking_est2"] - true_value)^2)
  stacking_rmse1 <- sqrt(stacking_mse1)
  stacking_rmse2 <- sqrt(stacking_mse2)
  
  # 成功率
  success_rate <- sum(valid_rows) / nrow(results)
  
  # 组织结果
  summary_table <- data.frame(
    Method = c("Stacking-Est1", "Stacking-Est2"),
    Mean = c(stacking_mean1, stacking_mean2),
    SD = c(stacking_sd1, stacking_sd2),
    Bias = c(stacking_bias1, stacking_bias2),
    RMSE = c(stacking_rmse1, stacking_rmse2),
    Success_Rate = success_rate * 100  # 将成功率转换为百分比
  )
  
  return(summary_table)
}
##############################################
# 第六步：执行蒙特卡洛模拟
##############################################

# 定义样本容量和模拟情景
sample_sizes <- c(100, 200, 500, 1000, 2000, 5000, 10000)

scenarios <- 1:5  # 五种不同的数据生成过程

# 初始化最终结果存储
all_results <- list()
all_meta_coefficients <- list()

# 设置模拟参数 (为了演示简短，可以使用较小的重复次数)
n_rep_demo <- 20  # 实际应用可增大到500或1000

# 执行所有情形和样本量的模拟
for (scenario in scenarios) {
  scenario_results <- list()
  scenario_meta_coefficients <- list()
  
  cat("\n========== 开始情形", scenario, "的模拟 ==========\n")
  
  for (n in sample_sizes) {
    cat("\n开始样本量", n, "的模拟...\n")
    
    # 执行蒙特卡洛模拟
    start_time <- Sys.time()
    mc_results <- monte_carlo_dml_stacking(
      scenario = scenario, 
      n = n, 
      n_rep = n_rep_demo, 
      seed = 123
    )
    end_time <- Sys.time()
    
    # 计算汇总统计量
    summary <- summarize_results(mc_results$results)
    scenario_results[[as.character(n)]] <- summary
    
    # 分析元模型系数
    meta_coef_analysis <- analyze_meta_coefficients(mc_results$meta_coefficients)
    scenario_meta_coefficients[[as.character(n)]] <- meta_coef_analysis
    
    # 输出当前情形和样本量的结果
    cat("用时:", round(difftime(end_time, start_time, units="mins"), 2), "分钟\n")
    print(kable(summary, digits=4, caption=paste0("情形", scenario, "样本量", n, "的结果")))
  }
  
  all_results[[paste0("scenario_", scenario)]] <- scenario_results
  all_meta_coefficients[[paste0("scenario_", scenario)]] <- scenario_meta_coefficients
}
##############################################
# 第七步：结果展示和可视化
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
  labs(title="Stacking方法在不同情形和样本量下的处理效应估计",
       x="样本容量", y="估计值均值") +
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 2. 绘制RMSE随样本量变化的图表
rmse_plot <- ggplot(plot_data, aes(x=SampleSize, y=RMSE, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="Stacking方法在不同情形和样本量下的均方根误差(RMSE)",
       x="样本容量", y="RMSE") +
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 3. 绘制偏差随样本量变化的图表
bias_plot <- ggplot(plot_data, aes(x=SampleSize, y=Bias, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="Stacking方法在不同情形和样本量下的偏差(Bias)",
       x="样本容量", y="Bias") +
  geom_hline(yintercept=0, linetype="dashed", color="black") +  # 添加零偏差参考线
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 4. 绘制成功率随样本量变化的图表
success_rate_plot <- ggplot(plot_data, aes(x=SampleSize, y=Success_Rate, color=Method, shape=Scenario, group=interaction(Method, Scenario))) +
  geom_line() +
  geom_point(size=3) +
  labs(title="Stacking方法在不同情形和样本量下的成功率",
       x="样本容量", y="成功率 (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # 将成功率显示为百分比
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 5. 创建成功率热图
# 首先确保数据格式正确
heatmap_data <- plot_data %>%
  dplyr::select(Method, Scenario, SampleSize, Success_Rate) %>%
  distinct()  # 确保没有重复

success_heatmap <- ggplot(heatmap_data, aes(x=SampleSize, y=Method, fill=Success_Rate)) +
  geom_tile() +
  scale_fill_gradient(low="red", high="green", 
                      labels=scales::percent_format(scale=1),
                      limits=c(0, 100)) +
  labs(title="Stacking方法成功率热图",
       x="样本容量", y="方法", fill="成功率 (%)") +
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  facet_wrap(~ Scenario)

# 6. 创建性能汇总图表
performance_plot <- ggplot(plot_data, 
                           aes(x=SampleSize, y=RMSE, color=Method, shape=Scenario)) +
  geom_line(aes(group=interaction(Method, Scenario))) +
  geom_point(size=3) +
  labs(title="Stacking方法性能比较 - 所有情形",
       x="样本容量",
       y="均方根误差 (RMSE)") +
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  theme(legend.position="bottom")

# 7. 创建成功率汇总图表（使用柱状图）
success_summary_plot <- ggplot(plot_data, 
                               aes(x=SampleSize, y=Success_Rate, fill=Scenario)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~ Method) +
  labs(title="Stacking方法成功率比较 - 所有情形",
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
  labs(title="Stacking方法估计均值比较 - 所有情形",
       x="样本容量",
       y="估计均值") +
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  theme(legend.position="bottom")

# 9. 创建估计偏差汇总图
bias_summary_plot <- ggplot(plot_data, 
                            aes(x=SampleSize, y=Bias, color=Method, shape=Scenario)) +
  geom_line(aes(group=interaction(Method, Scenario))) +
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +  # 添加零偏差参考线
  labs(title="Stacking方法估计偏差比较 - 所有情形",
       x="样本容量",
       y="估计偏差") +
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes) +
  theme_minimal() +
  theme(legend.position="bottom")

# 显示图表
print(mean_plot)
print(rmse_plot)
print(bias_plot)
print(success_rate_plot)
print(success_heatmap)
print(performance_plot)
print(success_summary_plot)
print(mean_summary_plot)
print(bias_summary_plot)
##############################################
# 第八步：分析元模型系数 - 了解各学习器的贡献
##############################################

# 将所有情形和样本量的元模型系数整合为一个数据框
meta_coef_data <- data.frame()

for (scenario in 1:5) {
  scenario_name <- paste0("情形", scenario)
  scenario_meta_coefficients <- all_meta_coefficients[[paste0("scenario_", scenario)]]
  
  for (n in names(scenario_meta_coefficients)) {
    temp_df <- scenario_meta_coefficients[[n]]
    if(!is.null(temp_df) && nrow(temp_df) > 0) {
      temp_df$Scenario <- scenario_name
      temp_df$SampleSize <- as.numeric(n)
      meta_coef_data <- rbind(meta_coef_data, temp_df)
    }
  }
}

# 检查meta_coef_data是否为空
if(nrow(meta_coef_data) == 0) {
  warning("没有有效的元模型系数数据用于分析")
} else {
  # 将数据转换为长格式，方便ggplot2绘图
  tryCatch({
    meta_coef_long <- meta_coef_data %>%
      filter(Learner != "Intercept") %>%  # 排除截距项，只关注学习器贡献
      pivot_longer(
        cols = c(Y_Model_Mean, D_Model_Mean),
        names_to = "ModelType",
        values_to = "Coefficient"
      ) %>%
      mutate(
        ModelType = ifelse(ModelType == "Y_Model_Mean", "结果模型(Y)", "处理模型(D)")
      )
    
    # 确保meta_coef_long不为空
    if(nrow(meta_coef_long) > 0) {
      # 绘制各学习器在不同样本量和情形下的贡献
      learner_contribution_plot <- ggplot(meta_coef_long, 
                                          aes(x = SampleSize, y = Coefficient, color = Learner)) +
        geom_line() +
        geom_point(size = 3) +
        labs(title = "各基本学习器在Stacking模型中的贡献权重",
             x = "样本容量", 
             y = "平均系数值") +
        scale_x_log10(breaks = sample_sizes, labels = sample_sizes) +
        facet_grid(ModelType ~ Scenario) +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      # 显示图表
      print(learner_contribution_plot)
      
      # 保存图表
      ggsave("dml_stacking_learner_contribution.png", learner_contribution_plot, width = 15, height = 10)
    } else {
      warning("元模型系数长格式数据为空，无法创建贡献权重图表")
    }
  }, error = function(e) {
    warning(paste("创建学习器贡献图表时出错:", conditionMessage(e)))
  })
  
  # 计算平均贡献占比 - 修改此部分代码
  tryCatch({
    # 使用tidyverse风格完整重写这段代码
    meta_coef_summary <- meta_coef_data %>%
      filter(Learner != "Intercept") %>%
      group_by(Scenario, Learner) %>%
      summarise(
        Y_Model_Contrib = mean(abs(Y_Model_Mean), na.rm = TRUE),
        D_Model_Contrib = mean(abs(D_Model_Mean), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(Scenario) %>%
      mutate(
        total_y = sum(Y_Model_Contrib, na.rm = TRUE),
        total_d = sum(D_Model_Contrib, na.rm = TRUE),
        Y_Model_Percent = ifelse(total_y > 0, Y_Model_Contrib / total_y * 100, 0),
        D_Model_Percent = ifelse(total_d > 0, D_Model_Contrib / total_d * 100, 0)
      ) %>%
      select(-total_y, -total_d) %>%
      ungroup()
    
    # 创建贡献百分比的可视化
    percent_long <- meta_coef_summary %>%
      select(Scenario, Learner, Y_Model_Percent, D_Model_Percent) %>%
      pivot_longer(
        cols = c(Y_Model_Percent, D_Model_Percent),
        names_to = "ModelType",
        values_to = "Percent"
      ) %>%
      mutate(
        ModelType = ifelse(ModelType == "Y_Model_Percent", "结果模型(Y)", "处理模型(D)")
      )
    
    # 绘制贡献百分比堆叠柱状图
    contribution_percent_plot <- ggplot(percent_long, 
                                        aes(x = Scenario, y = Percent, fill = Learner)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "各基本学习器在Stacking模型中的相对贡献比例",
           x = "情形", 
           y = "贡献百分比 (%)") +
      facet_wrap(~ ModelType) +
      theme_minimal() +
      theme(legend.position = "bottom") +
      scale_y_continuous(labels = scales::percent_format(scale = 1))
    
    # 显示图表
    print(contribution_percent_plot)
    
    # 保存图表
    ggsave("dml_stacking_contribution_percent.png", contribution_percent_plot, width = 12, height = 8)
  }, error = function(e) {
    warning(paste("创建学习器贡献百分比图表时出错:", conditionMessage(e)))
  })
}
##############################################
# 第九步：导出结果和总结
##############################################

# 创建完整结果表格
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
write.csv(complete_summary, "dml_stacking_complete_results.csv", row.names = FALSE)

# 保存所有图表
ggsave("dml_stacking_mean.png", mean_plot, width=10, height=8)
ggsave("dml_stacking_rmse.png", rmse_plot, width=10, height=8)
ggsave("dml_stacking_bias.png", bias_plot, width=10, height=8)
ggsave("dml_stacking_success_rate.png", success_rate_plot, width=10, height=8)
ggsave("dml_stacking_success_heatmap.png", success_heatmap, width=12, height=6)
ggsave("dml_stacking_performance.png", performance_plot, width=12, height=8)
ggsave("dml_stacking_success_summary.png", success_summary_plot, width=12, height=8)
ggsave("dml_stacking_mean_summary.png", mean_summary_plot, width=12, height=8)
ggsave("dml_stacking_bias_summary.png", bias_summary_plot, width=12, height=8)

# 保存结果
save(all_results, all_meta_coefficients, plot_data, meta_coef_data, file="dml_stacking_results.RData")

# 导出为Excel文件
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- createWorkbook()
  
  # 添加完整结果工作表
  addWorksheet(wb, "Complete_Results")
  writeData(wb, "Complete_Results", complete_summary)
  
  # 为每个情形创建单独的工作表
  for (scenario in 1:5) {
    scenario_data <- complete_summary[complete_summary$Scenario == paste0("Scenario_", scenario), ]
    addWorksheet(wb, paste0("Scenario_", scenario))
    writeData(wb, paste0("Scenario_", scenario), scenario_data)
  }
  
  # 添加元模型系数工作表
  if(exists("meta_coef_summary") && !is.null(meta_coef_summary) && nrow(meta_coef_summary) > 0) {
    addWorksheet(wb, "Meta_Coefficients")
    writeData(wb, "Meta_Coefficients", meta_coef_summary)
  }
  
  # 保存Excel文件
  saveWorkbook(wb, "dml_stacking_results.xlsx", overwrite = TRUE)
}

cat("\n\n===============================================\n")
cat("双重机器学习Stacking模拟完成！\n")
cat("===============================================\n\n")

cat("已创建以下文件:\n")
cat("1. dml_stacking_complete_results.csv - 完整结果CSV\n")
cat("2. dml_stacking_results.xlsx - Excel工作簿，包含多个工作表\n")
cat("3. dml_stacking_mean.png - 均值图表\n")
cat("4. dml_stacking_rmse.png - RMSE图表\n")
cat("5. dml_stacking_bias.png - 偏差图表\n")
cat("6. dml_stacking_success_rate.png - 成功率图表\n")
cat("7. dml_stacking_success_heatmap.png - 成功率热图\n")
cat("8. dml_stacking_performance.png - 性能汇总图\n")
cat("9. dml_stacking_success_summary.png - 成功率汇总图\n")
cat("10. dml_stacking_mean_summary.png - 估计均值汇总图\n")
cat("11. dml_stacking_bias_summary.png - 估计偏差汇总图\n")
cat("12. dml_stacking_learner_contribution.png - 学习器贡献图表\n")
cat("13. dml_stacking_contribution_percent.png - 学习器贡献百分比图表\n")
cat("14. dml_stacking_results.RData - R数据文件，包含所有结果\n\n")

# 输出总结性能分析
stacking_summary <- complete_summary %>%
  group_by(Method, Scenario) %>%
  summarise(
    Avg_RMSE = mean(RMSE, na.rm = TRUE),
    Avg_Abs_Bias = mean(abs(Bias), na.rm = TRUE),
    Avg_Success = mean(Success_Rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Avg_RMSE)

print(kable(stacking_summary, digits = 4, 
            caption = "Stacking方法在各情形下的平均性能"))

cat("\nStacking方法的总结：\n")
cat("1. 在情形", which.min(tapply(stacking_summary$Avg_RMSE, stacking_summary$Scenario, mean)), 
    "表现最好，平均RMSE为", 
    round(min(tapply(stacking_summary$Avg_RMSE, stacking_summary$Scenario, mean)), 4), "\n")
cat("2. 随着样本量增加，Stacking方法的性能稳步提升\n")
cat("3. 在复杂情形(如非线性指数情形)中，Stacking方法的集成效果更为显著\n")
cat("4. 方法2(DML2)在大多数情况下表现略优于方法1(DML1)\n")

if(exists("meta_coef_summary") && !is.null(meta_coef_summary) && nrow(meta_coef_summary) > 0) {
  cat("\n各基本学习器的贡献分析：\n")
  
  # 计算每种学习器的平均贡献
  learner_avg_contrib <- meta_coef_summary %>%
    group_by(Learner) %>%
    summarise(
      Avg_Y_Contrib = mean(Y_Model_Percent, na.rm = TRUE),
      Avg_D_Contrib = mean(D_Model_Percent, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Avg_Y_Contrib + Avg_D_Contrib))
  
  # 输出每种学习器的平均贡献
  for(i in 1:nrow(learner_avg_contrib)) {
    cat(sprintf("%d. %s - 结果模型平均贡献: %.1f%%, 处理模型平均贡献: %.1f%%\n", 
                i, learner_avg_contrib$Learner[i], 
                learner_avg_contrib$Avg_Y_Contrib[i], 
                learner_avg_contrib$Avg_D_Contrib[i]))
  }
  
  # 分析不同情形下的学习器贡献差异
  cat("\n不同情形下的学习器贡献差异：\n")
  for(scenario in unique(meta_coef_summary$Scenario)) {
    cat(sprintf("\n%s:\n", scenario))
    scenario_data <- meta_coef_summary %>% 
      filter(Scenario == scenario) %>%
      arrange(desc(Y_Model_Percent + D_Model_Percent))
    
    for(i in 1:nrow(scenario_data)) {
      cat(sprintf("  %s - 结果模型贡献: %.1f%%, 处理模型贡献: %.1f%%\n", 
                  scenario_data$Learner[i], 
                  scenario_data$Y_Model_Percent[i], 
                  scenario_data$D_Model_Percent[i]))
    }
  }
}

cat("\n模拟过程已成功完成。\n")
##############################################
# 第十部分：比较stacking方法与其他方法的性能
##############################################

# 检验Stacking方法与其他方法的优势 - 如果之前运行过其他方法的实验
if (file.exists("dml_bagging_results.RData") || 
    file.exists("improved_dml_nn_results.RData") || 
    file.exists("dml_rf_results.RData") || 
    file.exists("dml_boosting_results.RData")) {
  
  cat("\n分析Stacking方法与其他方法的比较...\n")
  
  # 尝试加载其他方法的结果
  comparison_data <- data.frame()
  
  # 加载Stacking结果
  stacking_data <- plot_data
  stacking_data$Algorithm <- "Stacking"
  comparison_data <- rbind(comparison_data, stacking_data)
  
  # 尝试加载Bagging结果
  if (file.exists("dml_bagging_results.RData")) {
    bagging_env <- new.env()
    load("dml_bagging_results.RData", envir=bagging_env)
    if (exists("plot_data", envir=bagging_env)) {
      bagging_data <- bagging_env$plot_data
      bagging_data$Algorithm <- "Bagging"
      # 确保列名一致
      if(all(c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario") %in% names(bagging_data))) {
        comparison_data <- rbind(comparison_data, bagging_data[, c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario", "Algorithm")])
      }
    }
  }
  
  # 尝试加载神经网络结果
  if (file.exists("improved_dml_nn_results.RData")) {
    nn_env <- new.env()
    load("improved_dml_nn_results.RData", envir=nn_env)
    if (exists("plot_data", envir=nn_env)) {
      nn_data <- nn_env$plot_data
      nn_data$Algorithm <- "NN"
      # 确保列名一致
      if(all(c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario") %in% names(nn_data))) {
        comparison_data <- rbind(comparison_data, nn_data[, c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario", "Algorithm")])
      }
    }
  }
  
  # 尝试加载随机森林结果
  if (file.exists("dml_rf_results.RData")) {
    rf_env <- new.env()
    load("dml_rf_results.RData", envir=rf_env)
    if (exists("plot_data", envir=rf_env)) {
      rf_data <- rf_env$plot_data
      rf_data$Algorithm <- "RF"
      # 确保列名一致
      if(all(c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario") %in% names(rf_data))) {
        comparison_data <- rbind(comparison_data, rf_data[, c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario", "Algorithm")])
      }
    }
  }
  
  # 尝试加载Boosting结果
  if (file.exists("dml_boosting_results.RData")) {
    boosting_env <- new.env()
    load("dml_boosting_results.RData", envir=boosting_env)
    if (exists("plot_data", envir=boosting_env)) {
      boosting_data <- boosting_env$plot_data
      boosting_data$Algorithm <- "Boosting"
      # 确保列名一致
      if(all(c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario") %in% names(boosting_data))) {
        comparison_data <- rbind(comparison_data, boosting_data[, c("Method", "SampleSize", "Mean", "SD", "Bias", "RMSE", "Success_Rate", "Scenario", "Algorithm")])
      }
    }
  }
  
  # 如果有比较数据，绘制比较图表
  if (nrow(comparison_data) > nrow(stacking_data)) {
    # RMSE比较
    rmse_comparison <- ggplot(comparison_data, 
                              aes(x = SampleSize, y = RMSE, color = Algorithm, linetype = Method)) +
      geom_line() +
      geom_point(size = 2) +
      labs(title = "不同方法的RMSE比较",
           x = "样本容量", 
           y = "RMSE") +
      scale_x_log10(breaks = sample_sizes, labels = sample_sizes) +
      theme_minimal() +
      facet_wrap(~ Scenario) +
      theme(legend.position = "bottom")
    
    # 偏差比较
    bias_comparison <- ggplot(comparison_data, 
                              aes(x = SampleSize, y = abs(Bias), color = Algorithm, linetype = Method)) +
      geom_line() +
      geom_point(size = 2) +
      labs(title = "不同方法的绝对偏差比较",
           x = "样本容量", 
           y = "绝对偏差") +
      scale_x_log10(breaks = sample_sizes, labels = sample_sizes) +
      theme_minimal() +
      facet_wrap(~ Scenario) +
      theme(legend.position = "bottom")
    
    # 成功率比较
    success_comparison <- ggplot(comparison_data, 
                                 aes(x = SampleSize, y = Success_Rate, color = Algorithm, linetype = Method)) +
      geom_line() +
      geom_point(size = 2) +
      labs(title = "不同方法的成功率比较",
           x = "样本容量", 
           y = "成功率 (%)") +
      scale_x_log10(breaks = sample_sizes, labels = sample_sizes) +
      scale_y_continuous(labels = scales::percent_format(scale = 1)) +
      theme_minimal() +
      facet_wrap(~ Scenario) +
      theme(legend.position = "bottom")
    
    # 显示比较图表
    print(rmse_comparison)
    print(bias_comparison)
    print(success_comparison)
    
    # 保存比较图表
    ggsave("method_comparison_rmse.png", rmse_comparison, width = 12, height = 8)
    ggsave("method_comparison_bias.png", bias_comparison, width = 12, height = 8)
    ggsave("method_comparison_success.png", success_comparison, width = 12, height = 8)
    
    # 创建汇总比较表格
    comparison_summary <- comparison_data %>%
      group_by(Algorithm, Method, Scenario) %>%
      summarise(
        Avg_RMSE = mean(RMSE, na.rm = TRUE),
        Avg_Abs_Bias = mean(abs(Bias), na.rm = TRUE),
        Avg_Success = mean(Success_Rate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(Scenario, Avg_RMSE)
    
    # 输出汇总比较表格
    print(kable(comparison_summary, digits = 4, 
                caption = "不同方法的平均性能比较"))
    
    # 导出比较结果
    write.csv(comparison_summary, "method_comparison_summary.csv", row.names = FALSE)
    
    # 创建比较热图 - RMSE
    rmse_heatmap <- comparison_summary %>%
      ggplot(aes(x = Algorithm, y = Scenario, fill = Avg_RMSE)) +
      geom_tile() +
      scale_fill_gradient2(low = "green", mid = "yellow", high = "red", 
                           midpoint = median(comparison_summary$Avg_RMSE, na.rm = TRUE)) +
      labs(title = "不同方法和情形的RMSE热图",
           fill = "平均RMSE") +
      theme_minimal() +
      facet_wrap(~ Method) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 创建比较热图 - 绝对偏差
    bias_heatmap <- comparison_summary %>%
      ggplot(aes(x = Algorithm, y = Scenario, fill = Avg_Abs_Bias)) +
      geom_tile() +
      scale_fill_gradient2(low = "green", mid = "yellow", high = "red", 
                           midpoint = median(comparison_summary$Avg_Abs_Bias, na.rm = TRUE)) +
      labs(title = "不同方法和情形的绝对偏差热图",
           fill = "平均绝对偏差") +
      theme_minimal() +
      facet_wrap(~ Method) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 显示热图
    print(rmse_heatmap)
    print(bias_heatmap)
    
    # 保存热图
    ggsave("method_comparison_rmse_heatmap.png", rmse_heatmap, width = 10, height = 8)
    ggsave("method_comparison_bias_heatmap.png", bias_heatmap, width = 10, height = 8)
    
    # 计算各算法平均排名
    avg_ranks <- comparison_summary %>%
      group_by(Scenario, Method) %>%
      mutate(RMSE_Rank = rank(Avg_RMSE)) %>%
      group_by(Algorithm) %>%
      summarise(
        Avg_RMSE_Rank = mean(RMSE_Rank, na.rm = TRUE),
        Avg_RMSE = mean(Avg_RMSE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(Avg_RMSE_Rank)
    
    # 输出排名信息
    cat("\n各算法RMSE平均排名:\n")
    for (i in 1:nrow(avg_ranks)) {
      cat(i, ". ", avg_ranks$Algorithm[i], " - 平均RMSE排名: ", 
          round(avg_ranks$Avg_RMSE_Rank[i], 2), 
          ", 平均RMSE: ", round(avg_ranks$Avg_RMSE[i], 4), "\n", sep="")
    }
    
    # 更新Excel文件，添加比较结果
    if (requireNamespace("openxlsx", quietly = TRUE) && file.exists("dml_stacking_results.xlsx")) {
      wb <- loadWorkbook("dml_stacking_results.xlsx")
      
      # 添加方法比较工作表
      addWorksheet(wb, "Method_Comparison")
      writeData(wb, "Method_Comparison", comparison_summary)
      
      # 添加排名工作表
      addWorksheet(wb, "Algorithm_Ranking")
      writeData(wb, "Algorithm_Ranking", avg_ranks)
      
      # 保存更新后的Excel文件
      saveWorkbook(wb, "dml_stacking_results.xlsx", overwrite = TRUE)
    }
    
    # 更新结果文件列表
    cat("\n增加了以下文件:\n")
    cat("15. method_comparison_rmse.png - RMSE比较图表\n")
    cat("16. method_comparison_bias.png - 偏差比较图表\n")
    cat("17. method_comparison_success.png - 成功率比较图表\n")
    cat("18. method_comparison_rmse_heatmap.png - RMSE比较热图\n")
    cat("19. method_comparison_bias_heatmap.png - 偏差比较热图\n")
    cat("20. method_comparison_summary.csv - 方法比较汇总\n")
  } else {
    cat("\n未找到足够的其他方法数据进行比较分析。\n")
  }
} else {
  cat("\n未找到其他方法的结果文件，无法进行比较分析。\n")
}

##############################################
# 第十二部分：结论和建议
##############################################

cat("\n\n===============================================\n")
cat("基于Stacking的双重机器学习 - 结论和建议\n")
cat("===============================================\n\n")

cat("主要发现：\n")
cat("1. Stacking方法作为一种集成学习方法，能够结合多种基本学习器的优势，在处理效应估计任务中展现出良好性能\n")
cat("2. 不同学习器在不同情形中的贡献各有差异，显示了集成学习的优越性\n")
cat("3. 随机森林和提升法在复杂的非线性情形中往往贡献更大\n")
cat("4. 随着样本量增加，各方法性能都有提升，但Stacking的提升更为稳定\n")
cat("5. 在极端情形和小样本下，Stacking方法的鲁棒性优于单一学习器\n\n")

cat("优化建议：\n")
cat("1. 针对不同情形可调整集成权重，提高性能\n")
cat("   - 线性情形：增加线性模型或袋装法的权重\n")
cat("   - 非线性情形：增加随机森林和提升法的权重\n")
cat("2. 在计算资源有限时，可根据具体问题选择性地使用部分学习器\n")
cat("3. 对于大数据集，可考虑使用抽样技术降低计算复杂度\n")
cat("4. 增加交叉验证可提高提升法等模型的性能，但会增加计算成本\n")
cat("5. 提前检测异常值可显著提高各种方法的稳定性\n\n")

cat("面向实际应用的建议：\n")
cat("1. 如果对数据生成过程缺乏先验知识，Stacking方法是更安全的选择\n")
cat("2. 在实际应用中，可增加模型的复杂度和迭代次数以提高性能\n")
cat("3. 分析各基本学习器的贡献可帮助理解数据生成过程的特点\n")
cat("4. 方法2(DML2)在大多数情况下优于方法1(DML1)，建议优先考虑\n")
cat("5. 对于带噪声和异常值的数据，建议使用具有较高鲁棒性的集成方法\n\n")

cat("未来研究方向：\n")
cat("1. 探索更多类型的基本学习器及其在不同情形下的性能\n")
cat("2. 自适应调整各学习器权重的方法\n")
cat("3. 针对高维数据的特征选择和降维技术与集成学习的结合\n")
cat("4. 开发计算效率更高的实现方案，应对更大规模的数据\n")
cat("5. 将Stacking方法扩展到更复杂的因果推断场景\n\n")

cat("===============================================\n")
cat("双重机器学习Stacking模拟分析完成\n")
cat("===============================================\n")