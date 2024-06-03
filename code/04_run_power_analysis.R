# clear your environment and memory
rm(list=ls())
gc()

# read the exported files from 01_download_stack_save_NEON_data.R script
joined_fish <- readr::read_csv("./input/joined_fish_tables.csv")
predictions <- readr::read_csv("./input/predictions_comparisons_file.csv")

# set the time to see how long it takes
s = Sys.time()

# set the sites
sites <- c(unique(joined_fish$siteID))

three_pass_new_method <- predictions %>%
  filter(method == "Carle-Strube") %>%
  left_join(., reach_lengths, by = "reachID")


for(g in 1:length(sites)){
  
  tryCatch({
    
    # DATA FOR THTREE PASS #
    three_pass_data_wide_pwr <- three_pass_new_method %>%
      left_join(., total_prob_summary, by = c("reachID"), relationship = "many-to-many") %>%
      na.omit(.) %>%
      rename(cap_prob = value, cap_prob_low = .lower, cap_prob_high = .upper) %>%
      filter(siteID == sites[g]) %>%
      na.omit(.)
    
    mean_wetted_width_pwr <- mean_wetted_widths %>% filter(siteID == sites[g])
    
    data_for_model <- fuzzyjoin::difference_join(three_pass_data_wide_pwr, 
                                                 mean_wetted_width_pwr, 
                                                 by = c("year", "month"), 
                                                 max_dist = 1, mode = "left") %>%
      dplyr::select(-year.y, -month.y, -sd_wetted_width_m) %>%
      rename(month = month.x, year = year.x) %>%
      group_by_at(vars(-mean_wetted_width_m)) %>%
      summarise(mean_wetted_width_m = mean(mean_wetted_width_m)) %>%
      ungroup(.) %>%
      fill(measuredReachLength) |>
      arrange(year)
    
    ##Abundance model
    model.pop.threepass=glmer(estimate~(1|year)+(1|reach), 
                              data=data_for_model, family=poisson, 
                              offset=log(measuredReachLength/
                                           max(data_for_model$measuredReachLength)))
    
    ## Capture probability model
    data_for_model$p.logit=log(data_for_model$cap_prob_low/(1-data_for_model$cap_prob_low)) # capture probability on a logit scale
    data_for_model$WidthStd <- scale(data_for_model$mean_wetted_width_m) 	          # centered and scaled stream width
    model.p=lm(p.logit ~ WidthStd , data=data_for_model)
    summary(model.p)
    
    n.sims=250 				                                                              # the number of simulations
    n.passes=3 				                                                              # the number of passes
    
    ## each of the below are used to populate the empty data frame power.table
    simsN=sim(model.pop.threepass,n.sims) 		                                      # simulates 300 parameter sets of the abundance model
    simsP=sim(model.p, n.sims)   		                                                # simulates 300 parameter sets of the capture probability model
    p.vals=array(NA,c(n.sims,2))  		                                              # creates an empty array with two columns and 300 rows for p-values
    relgrad=array(NA,c(n.sims,2))		                                                # creates an empty array with two columns and 300 rows to monitor convergence
    
    ## create an empty df to rbind p-values for each setting later
    power.table <- data.frame(n.sites=integer(0),	                                  # populated with the number of sites
                              n.years=integer(0),			                              # populated with the number of years
                              rhat=numeric(0),			                                # populated with the percent population decline
                              power.firstPass=numeric(0),		                        # populated with the proportion of p-values < 0.05 from the first pass
                              power.cs=numeric(0),		                              # populated with the proportion of p-values < 0.05 from the three pass method
                              power.firstPass.cnvg=numeric(0),	                    # populated with value from monitoring of model convergence in the first pass
                              power.cs.cnvg=numeric(0))		                          # populated with value from monitoring of model convergence in the three pass method
    
    
    for(n.sites in c(1)){ 			                                                    # add  number of sites
      for(n.years in c(seq(from = 5, to = 30, by = 5))){ 		                        # add  number of years
        for(rhat in c(-0.05)){ 		                                                  # add  mean rate of population trends(s)
          
          ## start simulations
          for(s in 1:n.sims){  
            cat(paste('n.sites=', n.sites, 'n.years=', 
                      n.years,'rhat=', rhat, 'simulation', s, '\n'))
            mu=simsN@fixef[s,1]                                                     # overall mean abundance
            sigma.alpha=model.pop.threepass@theta[1]                                # SD of abundance among sites
            sigma.beta=model.pop.threepass@theta[2]                                 # SD of abundance among years
            
            ## stochastic factors affecting abundance
            alpha=rnorm(n.sites,0,sigma.alpha)  		                                # site random effect on abundance
            beta=rnorm(n.years,0,sigma.beta)		                                    # year random effect on abundance
            
            
            ## factors affecting capture probability
            gamma=simsP@coef[s,1]  			                                            # mean capture probability on the logit scale 
            delta=simsP@coef[s,2] 			                                            # coefficient of width on capture probability
            
            width.sd=array(rnorm(n.sites*n.years,mean=0,sd=1), 
                           dim=c(n.sites,n.years))                                  # standardized width
            sigma.eps=sd(model.p$residuals) 					                              # SD of capture probability on the logit scale 
            eps=array(rnorm(n.sites*n.years,mean=0,sd=sigma.eps),
                      dim=c(n.sites,n.years))                                       # sample random effect on capture probability
            
            ## setup a number of sites by number of years array for N, p, and lambda.
            N <- p <- lambda <- array(NA,dim=c(n.sites,n.years), 
                                      dimnames=list(paste("site",1:n.sites),
                                                    paste("year",1:n.years)))
            
            ## setup a number of sites by number of years array for y        
            y=array(NA,dim=c(n.sites,n.years,n.passes),
                    dimnames=list(paste("site",1:n.sites),
                                  paste("year",1:n.years),
                                  paste("pass",1:n.passes)))
            
            ## generate a range for the population decline (here set to a 5% range)  
            r=runif(n.sites, rhat-0.05, rhat+0.05) 		 
            trend=log(1+r)
            
            ## generate data
            for(i in 1:n.sites){
              for(j in 1:n.years){
                
                lambda[i,j]<-exp(mu+((trend[i])*(j-1)) + alpha[i] + beta[j])
                N[i,j]<-rpois(1,lambda[i,j])
                
                p[i,j]<-plogis(gamma+delta*width.sd[i,j] +
                                 eps[i,j])  
                
                y[i,j,1]<-rbinom(1,N[i,j],p[i,j])
                y[i,j,2]<-rbinom(1,N[i,j]-y[i,j,1],p[i,j])
                y[i,j,3]<-rbinom(1,N[i,j]-y[i,j,1]-y[i,j,2],p[i,j])
                
                ## replace all-0 samples (no fish in all three passes)
                if (sum(y[i,j,]) == 0){
                  i=i
                  for(j in 1:n.years){
                    lambda[i,]<-exp(mu+((trend[i])*(j-1))+beta[j]) 
                    N[i,j]<-rpois(1,lambda[i,j])
                    p[i,j]<-plogis(gamma + eps[i,j])
                    
                    
                    y[i,j,1]<-rbinom(1,N[i,j],p[i,j])
                    y[i,j,2]<-rbinom(1,N[i,j]-y[i,j,1],p[i,j])
                    y[i,j,3]<-rbinom(1,N[i,j]-y[i,j,1]-y[i,j,2],p[i,j])
                  }
                }
                
              }
            }
            
            # CStrub1<-apply(y, MARGIN=c(1,2), FUN=removal,
            #                method="CarleStrub", just.ests=TRUE)
            # 
            # CStrub2 <- reshape2::melt(CStrub1[1,,], value.name="CSEst",
            #                           varnames=c('site','year'))                    # estimated abundance from simulated data organized by number of sites and number of years
            # firstPass <- reshape2::melt(y[,,1], value.name="firstPcount",
            #                             varnames=c('site','year'))                  # count from simulated data for the first pass organized by number of sites and number of years
            # result <- merge(CStrub2, firstPass) %>%
            #   mutate(CSEst = abs(CSEst))
            
            
            #used to estimate the abundance (N) and capture probability (p) from the simulated data
            CStrub1<-apply(y, MARGIN=c(1,2), FUN=removal,
                           method="CarleStrub", just.ests=TRUE)

            CStrub2 <- reshape2::melt(CStrub1[1,,], value.name="CSEst", # change this in order to adapt the upper CI and lower CI of prediction model
                                      varnames=c('site','year'))%>%
              tibble::rownames_to_column(., "year") %>%
              separate(year, into = c("name", "year"), remove = F) %>%
              mutate(year = as.numeric(year)) %>%
              dplyr::select(-name) # estimated abundance from simulated data organized by number of sites and number of years

            firstPass <- reshape2::melt(y[,,1], value.name="firstPcount",
                                        varnames=c('site','year'))%>%
              tibble::rownames_to_column(., "year") %>%
              separate(year, into = c("name", "year"), remove = F) %>%
              mutate(year = as.numeric(year)) %>%
              dplyr::select(-name)                 # count from simulated data for the first pass organized by number of sites and number of years

            result <- merge(CStrub2, firstPass) %>%
              mutate(site = 1)
            
            
            # merged estimated abundance and count from the simulated data
            result$YearNum <- as.numeric(result$year)	
            result$YearNumStdz<-(result$YearNum)/(sd(result$YearNum))		            # standardize year 
            result <- arrange(result, site, YearNumStdz)				                    # organizing the data frame
            
            
            
            # generalized model for estimated abundance of the simulated data
            if(n.sites == 1){
              lm.sim.cs=glmer(CSEst ~ YearNumStdz + (1|year),
                              data=result, family=poisson)
            } else {
              lm.sim.cs=glmer(CSEst ~ YearNumStdz + (1|site) + (1|year),
                              data=result, family=poisson)
            }
            

            p.vals[s,2]=coef(summary(lm.sim.cs))[2,4]	                              # P-values for estimated abundance of the simulated data	                  # Check for convergence 
            relgrad[s,2]=max(abs(with(lm.sim.cs@optinfo$derivs,
                                      solve(Hessian,gradient))))                    # Check for convergence
          }
          
          ## append (rbind) settings and power
          temp<- data.frame(n.sites=n.sites, 
                            n.years=n.years,
                            rhat=rhat,
                            removal.05=length(p.vals[,2][p.vals[,2]<0.05])/n.sims,
                            removal.cnvg=length(relgrad[,2][relgrad[,2]<0.001])/n.sims)
          
          power.table <- rbind(power.table, temp)  
          
        }
      }
    }
    
    
    ## power.table stores all of the results after simulations have been completed for all combinations of number of sites, number of years, and for each rate of population decline.  In this example the table is populated with p-values < 0.05 for both single pass and three pass methods. 
    
    ##this saves the results as a .csv file called PowerSims to the working directory
    write.csv(power.table, file=paste0(getwd(),"/output/PowerSims_1reach_",sites[g],".csv"))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

e <- Sys.time()
t=e-s
print(t)