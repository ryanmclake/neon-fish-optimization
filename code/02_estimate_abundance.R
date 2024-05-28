
reach_ids <- c(unique(EF_joined$reachID))
out <- list()

for(i in 1:length(reach_ids)){
  
  for_vector <- EF_joined %>%
    filter(reachID == reach_ids[i])
  
  vector = c(for_vector$`1`, for_vector$`2`, for_vector$`3`)
  
  carlestrube_predict <- removal(vector)
  
  carlestrube_estimate <- data.frame(estimate = carlestrube_predict$est[1], 
                                     estimate_error = carlestrube_predict$est[2]) %>%
    tibble() %>%
    mutate(reachID = reach_ids[i]) %>%
    mutate(method = "Carle-Strube")
  
  moran_predict <- removal(vector,method="Moran")
  
  moran_estimate <- data.frame(estimate = moran_predict$est[1], 
                               estimate_error = moran_predict$est[2]) %>%
    tibble() %>%
    mutate(reachID = reach_ids[i]) %>%
    mutate(method = "Moran")
  
  schnute_predict <- removal(vector,method="Schnute")
  
  schnute_estimate <- data.frame(estimate = schnute_predict$est[1], 
                                 estimate_error = schnute_predict$est[2]) %>%
    tibble() %>%
    mutate(reachID = reach_ids[i]) %>%
    mutate(method = "Schnute")
  
  seber_predict <- removal(vector,method="Seber3")
  
  seber_estimate <- data.frame(estimate = seber_predict$est[1], 
                               estimate_error = seber_predict$est[2]) %>%
    tibble() %>%
    mutate(reachID = reach_ids[i]) %>%
    mutate(method = "Seber")
  
  burnham_predict <- removal(vector,method="Burnham")
  
  burnham_estimate <- data.frame(estimate = burnham_predict$est[1], 
                                 estimate_error = burnham_predict$est[2]) %>%
    tibble() %>%
    mutate(reachID = reach_ids[i]) %>%
    mutate(method = "Burnham")
  
  out[[i]] <- rbind(carlestrube_estimate, moran_estimate, 
                    schnute_estimate, seber_estimate, burnham_estimate)
  
  
}

three_pass_predict_methods_compare <- do.call("rbind", out) %>%
  mutate(estimate_low = estimate - estimate_error,
         estimate_high = estimate + estimate_error) %>%
  dplyr::select(-estimate_error) %>%
  dplyr::select(estimate, estimate_low, estimate_high, reachID, method)


three_pass_w_zero <- EF_joined %>%
  tidyr::gather(., passNumber, count, `1`:`3`, factor_key=TRUE) %>%
  group_by(reachID) %>%
  mutate(pass_sum = cumsum(count)) %>%
  mutate(pass = as.numeric(passNumber)) %>%
  group_modify(~ add_row(.x,.before=0)) %>%
  select(reachID, pass, pass_sum) %>%
  mutate(pass = ifelse(is.na(pass), 0, pass)) %>%
  mutate(pass_sum = ifelse(is.na(pass_sum), 0, pass_sum))

reach_ids <- c(unique(three_pass_w_zero$reachID))
out_hollings <- list()
three_prediction_raw <- list()

for(i in 1:length(reach_ids)){
  
  tryCatch({
    
    data_for_model <- three_pass_w_zero %>%
      filter(reachID == reach_ids[i])
    
    three_pass_predictor = nlsLM(pass_sum ~ (alpha*pass)/(beta+pass),
                                 start = list(alpha = data_for_model$pass_sum[4], 
                                              beta = 3),
                                 data = data_for_model,
                                 control = nls.lm.control(maxiter=1024))
    
    alpha_estimate <- data.frame(estimate = three_pass_predictor$m$getPars()[1], 
                                 estimate_error = summary(three_pass_predictor)$coefficients[,2][1]) %>%
      tibble() %>%
      mutate(reachID = reach_ids[i]) %>%
      mutate(method = "Type 2 Vmax") %>%
      mutate(estimate_low = estimate - estimate_error,
             estimate_high = estimate + estimate_error) %>%
      dplyr::select(-estimate_error) %>%
      dplyr::select(estimate, estimate_low, estimate_high, reachID, method)
    
    beta_estimate <- data.frame(estimate = three_pass_predictor$m$getAllPars()[2], 
                                estimate_error = summary(three_pass_predictor)$coefficients[,2][2]) %>%
      tibble() %>%
      mutate(reachID = reach_ids[i]) %>%
      mutate(method = "Beta 3-Pass") %>%
      mutate(estimate_low = estimate - estimate_error,
             estimate_high = estimate + estimate_error) %>%
      dplyr::select(-estimate_error) %>%
      dplyr::select(estimate, estimate_low, estimate_high, reachID, method)
    
    three_pass_alpha <- rnorm(500, three_pass_predictor$m$getAllPars()[1], summary(three_pass_predictor)$coefficients[,2][1])
    three_pass_beta <- rnorm(500, three_pass_predictor$m$getAllPars()[2], summary(three_pass_predictor)$coefficients[,2][2])
    
    three_pass_parms <- data.frame(alpha = three_pass_alpha, 
                                   beta = three_pass_beta) %>%
      tibble(.)
    
    
    passes <- seq(from = 3, to = 50, by = 1)
    prediction <- list()
    
    for(s in 1:length(passes)){
      
      prediction_raw <- predict_function(x = passes[s],
                                         beta = three_pass_parms$beta,
                                         alpha = three_pass_parms$alpha) %>%
        data.frame(abundance = ., pass = passes[s]) %>%
        tibble(.) %>%     
        mutate(abundance = round(abundance, digits = 0)) %>%
        mutate(abundance = ifelse(abundance < 0, 0, abundance)) %>%
        mutate(iteration = rep(row_number()))
      
      prediction[[s]] <- prediction_raw
      
    }
    
    three_prediction_raw[[i]] <- do.call("rbind", prediction) %>%
      mutate(reachID = reach_ids[i])
    
    hollings_predict <- do.call("rbind", prediction) %>%
      filter(pass == 50) %>%
      summarise(estimate = mean(abundance, na.rm = TRUE),
                estimate_low = quantile(abundance, probs = 0.05),
                estimate_high = quantile(abundance, probs = 0.95)) %>%
      mutate(reachID = reach_ids[i]) %>%
      mutate(method = "Type 2 Posterior Prediction 3-Pass")
    
    out_hollings[[i]] <-  rbind(alpha_estimate, hollings_predict, beta_estimate)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

three_pass_prediction <- do.call("rbind", out_hollings) %>%
  dplyr::select(estimate, estimate_low, estimate_high, reachID, method) %>%
  arrange(reachID)

prediction_comparisons <- bind_rows(three_pass_prediction, three_pass_predict_methods_compare) %>%
  arrange(reachID)

