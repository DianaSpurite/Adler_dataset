#Setting up workspace and data----------------------------------------------------------------------------------------

#Project spas
#By Diana Spurite
#diana.spurite@posteo.de


ks_spas <- read.csv("ks_spas.csv")
spo_asp <- read.csv("data_df.csv")


spo_asp <- spo_asp %>%
  mutate(basalArea_genet = ks_spas$basalArea_genet)




surv        <- subset( spo_asp, !is.na( survives ) ) %>%
  subset( basalArea_genet  != 0 ) %>%
  select( quad, year, track_id,
          basalArea_genet, logsize_t0, logsize_t1,
          survives, size_t1 )


grow        <- spo_asp %>%
  subset( basalArea_genet  != 0) %>%
  subset( size_t1 != 0) %>%
  select( quad, year, track_id,
          basalArea_genet, logsize_t0, logsize_t1,
          survives, size_t1)



quad_df     <- spo_asp %>%
  group_by( species, quad, year ) %>%
  summarise( totParea = sum( basalArea_genet ) ) %>%
  ungroup



group_df    <- quad_df %>%
  group_by( species, year ) %>%
  summarise( Gcov = mean( totParea ) ) %>%
  ungroup



cover_df     <-  left_join( quad_df, group_df ) %>%
  mutate( year = year + 1 ) %>%
  mutate( year = as.integer( year ) ) %>%
  drop_na()



recr_df     <- spo_asp %>%
  group_by( species, quad, year ) %>%
  summarise( NRquad   = sum( recruit, na.rm=T ) ) %>%
  ungroup



recr        <- left_join( cover_df, recr_df ) %>%
  drop_na



write.csv( surv, "data/spo_asp/survival_df.csv" )
write.csv( grow, "data/spo_asp/growth_df.csv" )
write.csv( recr, "data/spo_asp/recruitment_df.csv" )




#_______________________________________________________________________________
#### 7. Plotting the raw data ####


par( mfrow = c( 1, 2 ), mar = c( 3.5, 3.5, 1, 0.2 ), mgp = c( 2, 0.7, 0 ), cex = 0.8 )

hist( spo_asp$basalArea_genet, main = "Histogram of size at time t0", xlab = "Size at time t0" )
hist( spo_asp$size_t1, main = "Histogram of size at time t1", xlab = "Size at time t1" )



par( mfrow = c( 1, 2 ), mar = c( 3.5, 3.5, 1, 0.2 ), mgp = c( 2, 0.7, 0 ), cex = 0.8 )

hist( spo_asp$logsize_t0, main = "Histogram of log-transformed size at time t0", xlab = "log( Size at time t0 )" )
hist( log( spo_asp$size_t1 ), main = "Histogram of log-transformed size at time t1", xlab = "log( Size at time t1 )" )



h    <- ( max( spo_asp$logsize_t0, na.rm = T ) - min( spo_asp$logsize_t0, na.rm = T ) ) / 200
lwr  <- min( spo_asp$logsize_t0, na.rm = T ) + ( h * c( 0:( 200 - 1 ) ) )
upr  <- lwr + h
mid  <- lwr + ( 1/2 * h )

binned_prop <- function( lwr_x, upr_x, response ){
  
  id  <- which( spo_asp$logsize_t0 > lwr_x & spo_asp$logsize_t0 < upr_x )
  tmp <- spo_asp[id,]
  
  if( response == 'prob' ){   return( sum( tmp$survives, na.rm = T ) / nrow( tmp ) ) }
  if( response == 'n_size' ){ return( nrow( tmp ) ) }
  
}

y_binned <- Map( binned_prop, lwr, upr, 'prob' ) %>% unlist
x_binned <- mid
y_n_size <- Map( binned_prop, lwr, upr, 'n_size' ) %>% unlist

surv_binned <- data.frame( xx  = x_binned,
                           yy  = y_binned,
                           nn  = y_n_size ) %>%
  setNames( c( 'logsize_t0', 'survives', 'n_size' ) )


ggplot( data  = surv_binned, aes( x = logsize_t0, y = survives ) ) +
  geom_point( alpha = 1,
              pch   = 16,
              size  = 1,
              color = 'red' ) +
  theme_bw( ) +
  theme( axis.text = element_text( size = 8 ),
         title     = element_text( size = 10 ) ) +
  labs( x = expression( 'log(size)'[t0] ),
        y = expression( 'Survival to time t1' ) )



ggplot(data  = grow, aes( x = logsize_t0, y = log( size_t1 ) ) ) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 0.7,
              color = 'red' ) +
  theme_bw( ) +
  theme( axis.text     = element_text( size   = 8 ),
         title         = element_text( size   = 10 ) ) +
  labs( x = expression( 'log( size )'[t0] ),
        y = expression( 'log( size )'[t1] ) )






#_______________________________________________________________________________
#### 8. Fitting vital rate models for the mean IPM ####

grow_df      <- grow %>%
  mutate( logsize.t0 = log( basalArea_genet ),
          logsize.t1 = log( size_t1 ) )

surv_df      <- surv %>%
  mutate( logsize = log( basalArea_genet ) )



# Growth model

gr_mod_mean <- lm( logsize.t1 ~ logsize.t0, data = grow_df)



grow_df$pred <- predict( gr_mod_mean, type = "response" )

grow_line <- ggplot( grow_df, aes( x = logsize.t0, y = logsize.t1 ) ) +
  geom_point( ) +
  geom_abline( aes( intercept = coef( gr_mod_mean )[1],
                    slope     = coef( gr_mod_mean )[2] ),
               color     = 'red',
               lwd       = 2 )

grow_pred <- ggplot( grow_df, aes( x = logsize.t1, y = pred ) ) +
  geom_point( ) +
  geom_abline( aes( intercept = 0,
                    slope = 1 ),
               color = "red",
               lwd = 2 )


grow_line + grow_pred + plot_layout( )



x         <- fitted( gr_mod_mean )
y         <- resid( gr_mod_mean )^2
gr_var_m  <- nls( y ~ a * exp( b * x ), start = list( a = 1, b = 0 ) )


# Survival model

su_mod_mean <- glm( survives ~ logsize, data = surv_df, family = "binomial" )



surv_x <- seq( min( surv_df$logsize ), max( surv_df$logsize ), length.out = 100)
surv_pred <- boot::inv.logit( coef( su_mod_mean )[1] + coef( su_mod_mean )[2] * surv_x )

surv_pred_df <- data.frame( logsize = surv_x, survives = surv_pred )

surv_line <- ggplot( ) +
  geom_jitter( data = surv_df, aes( x        = logsize,
                                    y        = survives ),
               alpha    = 0.25,
               width    = 0,
               height   = 0.25 ) +
  geom_line( data = surv_pred_df, aes( x     = logsize,
                                       y     = survives ),
             color = 'red',
             lwd   = 2 )

surv_bin <- ggplot( ) +
  geom_point( data = surv_binned, aes( x     = logsize_t0,
                                       y     = survives ) ) +
  geom_line( data = surv_pred_df, aes( x     = logsize,
                                       y     = survives ),
             color = 'red',
             lwd   = 2 )


surv_line + surv_bin + plot_layout( )



# Recruitment model

rec_mod_mean <- MASS::glm.nb( NRquad ~ 1, data = recr_df )



recr_df        <- recr_df %>%
  mutate( pred_mod_mean = predict( rec_mod_mean, type = "response" ) )

rec_sums_df_m <- recr_df %>%
  summarize( NRquad    = sum( NRquad ),
             pred_mod_mean = sum( pred_mod_mean ) )



indiv_m <- surv_df %>%
  summarize( n_adults = n( ) )



repr_pc_m <- indiv_m %>%
  bind_cols( rec_sums_df_m ) %>%
  mutate( repr_pc_mean   = pred_mod_mean / n_adults ) %>%
  mutate( repr_pc_obs    = NRquad / n_adults ) %>%
  drop_na



repr_pc_m


# Exporting parameter estimates

grow_fe <- data.frame( coefficient = names( coef( gr_mod_mean ) ),
                       value       = coef( gr_mod_mean ) )
var_m   <- data.frame( coefficient = names( coef( gr_var_m ) ),
                       value       = coef( gr_var_m ) )

grow_out <- Reduce( function(...) rbind(...), list( grow_fe, var_m ) ) %>%
  mutate( coefficient = as.character( coefficient ) ) %>%
  mutate( coefficient = replace( coefficient,
                                 grepl( "Intercept", coefficient ),
                                 "b0" ) )

write.csv( grow_out, "data/spo_asp/grow_pars.csv", row.names = F )



surv_fe <- data.frame( coefficient = names( coef( su_mod_mean ) ),
                       value       = coef( su_mod_mean ) )

surv_out<- Reduce( function(...) rbind(...), list( surv_fe ) ) %>%
  mutate( coefficient = as.character( coefficient ) ) %>%
  mutate( coefficient = replace( coefficient,
                                 grepl( "Intercept", coefficient ),
                                 "b0" ) )

write.csv( surv_out, "data/spo_asp/surv_pars.csv", row.names = F )



recSize <- spo_asp %>% subset( recruit == 1)

others  <- data.frame( coefficient = c( "rec_siz", "rec_sd",
                                        "max_siz", "min_siz",
                                        "fecu_b0" ),
                       value       = c( mean( log( recSize$size_t0 ) ),
                                        sd( log( recSize$size_t0 ) ),
                                        grow_df$logsize.t0 %>% max,
                                        grow_df$logsize.t0 %>% min,
                                        repr_pc_m$repr_pc_mean ) )

write.csv( others, "data/spo_asp/other_pars.csv", row.names = F )




#_______________________________________________________________________________
#### 9. Building the IPM from scratch ####

extr_value <- function(x, field){ subset( x, coefficient == field )$value }

pars_mean  <- list(  surv_b0 = extr_value( surv_out, "b0" ),
                     surv_b1 = extr_value( surv_out, "logsize" ),
                     grow_b0 = extr_value( grow_out, "b0" ),
                     grow_b1 = extr_value( grow_out, "logsize.t0" ),
                     a       = extr_value( grow_out, "a" ),
                     b       = extr_value( grow_out, "b" ),
                     fecu_b0 = extr_value( others, "fecu_b0" ),
                     recr_sz = extr_value( others, "rec_siz" ),
                     recr_sd = extr_value( others, "rec_sd" ),
                     L       = extr_value( others, "min_siz" ),
                     U       = extr_value( others, "max_siz" ),
                     mat_siz = 200 )

pars       <- pars_mean

write.csv( pars_mean, "data/spo_asp/pars_mean.csv", row.names = F )



lm( logsize.t1 ~ logsize.t0, data = grow_df )
nls( y ~ a * exp( b * x ), start = list( a = 1, b = 0 ) )



# Function describing standard deviation of growth model
grow_sd <- function( x, pars ) {
  pars$a * ( exp( pars$b * x ) ) %>% sqrt
}

# Function describing growth from size x to size y
gxy <- function( x, y, pars ) {
  return( dnorm( y,  mean = pars$grow_b0 + pars$grow_b1*x,
                 sd   = grow_sd( x, pars ) ) )
}



glm( survives  ~ logsize, data = surv_df, family = 'binomial' )



inv_logit <- function( x ) { exp( x ) / ( 1 + exp( x ) ) }

sx <- function( x, pars ) {
  return( inv_logit( pars$surv_b0 + pars$surv_b1 * x ) )
}



pxy <- function( x, y, pars ) {
  return( sx( x, pars ) * gxy( x, y, pars ) )
}



fy <- function( y, pars, h ){
  n_recr  <- pars$fecu_b0
  recr_y  <- dnorm( y, pars$recr_sz, pars$recr_sd ) * h
  recr_y  <- recr_y / sum( recr_y )
  f       <- n_recr * recr_y
  return( f )
  
}



kernel <- function( pars ) {
  
  n   <- pars$mat_siz                      # number of bins over which to integrate
  L   <- pars$L                            # lower limit of integration
  U   <- pars$U                            # upper limit of integration
  h   <- ( U - L ) / n                     # bin size
  b   <- L + c( 0:n ) * h                  # lower boundaries of bins
  y   <- 0.5 * ( b[1:n] + b[2:( n + 1 )] ) # midpoints of bins
  
  # Fertility matrix
  Fmat        <- matrix( 0, n, n )
  Fmat[]      <- matrix( fy( y, pars, h ), n, n )
  
  # Survival vector
  Smat   <- c( )
  Smat   <- sx( y, pars )
  
  # Growth matrix
  Gmat   <- matrix( 0, n, n )
  Gmat[] <- t( outer( y, y, gxy, pars ) ) * h
  
  # Growth/survival transition matrix
  Tmat   <- matrix( 0, n, n )
  
  # Correct for eviction of offspring
  for( i in 1:( n / 2 ) ) {
    Gmat[1,i] <- Gmat[1,i] + 1 - sum( Gmat[,i] )
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Correct eviction of large adults
  for( i in ( n / 2 + 1 ):n ) {
    Gmat[n,i] <- Gmat[n,i] + 1 - sum( Gmat[,i] )
    Tmat[,i]  <- Gmat[,i] * Smat[i]
  }
  
  # Full Kernel is simply a summation of fertility and transition matrices
  k_yx <- Fmat + Tmat
  
  return( list( k_yx    = k_yx,
                Fmat    = Fmat,
                Tmat    = Tmat,
                Gmat    = Gmat,
                meshpts = y ) )
  
}



lambda_ipm <- function( i ) {
  return( Re( eigen( kernel( i )$k_yx )$value[1] ) )
}

lam_mean <- lambda_ipm( pars_mean )
lam_mean




#_______________________________________________________________________________
#### 10. Building the IPM with ipmr ####


proto_ipm_p <- init_ipm( sim_gen   = "simple",
                         di_dd     = "di",
                         det_stoch = "det" ) %>%
  define_kernel(
    
    name      = "P",
    
    family    = "CC",
    
    formula   = s * g,
    
    s         = plogis( surv_b0 +
                          surv_b1 * size_1 ),
    
    g         = dnorm( size_2, mu_g, grow_sig ),
    mu_g      = grow_b0 + grow_b1 * size_1,
    grow_sig  = sqrt( a * exp( b * size_1 ) ),
    
    data_list = pars_mean,
    states    = list( c( 'size' ) ),
    
    evict_cor = TRUE,
    evict_fun = truncated_distributions( fun    = 'norm',
                                         target = 'g' )
    
  ) %>%
  define_impl(
    
    make_impl_args_list(
      kernel_names = c( "P" ),
      int_rule     = rep( "midpoint", 1 ),
      state_start  = rep( "size", 1 ),
      state_end    = rep( "size", 1 )
    )
    
  ) %>%
  define_domains(
    
    size = c( pars_mean$L,
              pars_mean$U,
              pars_mean$mat_siz
    )
    
  ) %>%
  define_pop_state(
    n_size     = rep( 1 / 200, 200 )
  )

ipmr_p <- make_ipm( proto_ipm  = proto_ipm_p,
                    iterations = 200 )

lambda( ipmr_p )



plot( ipmr_p )



proto_ipm_pf <- init_ipm( sim_gen   = "simple",
                          di_dd     = "di",
                          det_stoch = "det" ) %>%
  
  define_kernel(
    name      = "P",
    family    = "CC",
    formula   = s * g,
    s         = plogis( surv_b0 +
                          surv_b1 * size_1),
    g         = dnorm( size_2, mu_g, grow_sig ),
    mu_g      = grow_b0 + grow_b1 * size_1,
    grow_sig  = sqrt( a * exp( b * size_1 ) ),
    data_list = pars_mean,
    states    = list( c( 'size' ) ),
    evict_cor = TRUE,
    evict_fun = truncated_distributions( fun    = 'norm',
                                         target = 'g' )
  ) %>%
  
  define_kernel(
    name      = 'F',
    family    = 'CC',
    formula   = fecu_b0 * r_d,
    r_d       = dnorm( size_2, recr_sz, recr_sd ),
    data_list = pars_mean,
    states    = list( c( 'size' ) ),
    evict_cor = TRUE,
    evict_fun = truncated_distributions( "norm", "r_d" )
  ) %>%
  
  define_impl(
    make_impl_args_list(
      kernel_names = c( "P", "F" ),
      int_rule     = rep( "midpoint", 2 ),
      state_start  = rep( "size", 2 ),
      state_end    = rep( "size", 2 )
    )
  ) %>%
  
  define_domains(
    size = c(pars_mean$L,
             pars_mean$U,
             pars_mean$mat_siz
    )
  ) %>%
  
  define_pop_state(
    n_size = rep( 1 / 200, 200 )
  )

ipmr_pf <- make_ipm( proto_ipm = proto_ipm_pf,
                     iterations = 200 )
lam_mean_ipmr <- lambda( ipmr_pf )
lam_mean_ipmr



##   lambda
## 1.099794



plot( ipmr_pf )
