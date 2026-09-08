# Local country-level HIV incidence extraction, fitting and prediction.
hiv_incidence_data <- function(path) {
  x <- as.data.frame(suppressMessages(readxl::read_excel(path, sheet="Data", skip=1, guess_max=100000)))
  x <- x[x$Type=="Country" & x$Sex=="Both" & x$Indicator==
    "Estimated incidence rate (new HIV infection per 1,000 uninfected population)",]
  x$year <- as.integer(x$Year)
  x$value <- suppressWarnings(as.numeric(x$Value))
  x$censored <- grepl("^<",x$Value)
  stopifnot(all(x$Value[x$censored]=="<0.01"), !anyDuplicated(x[c("ISO3","year","Age")]))
  a <- x[x$Age=="Age 15-19" & (is.finite(x$value) | x$censored),]
  c <- x[x$Age=="Age 0-14",]
  d <- data.frame(iso3=a$ISO3, country_name=a$`Country/Region`, region=a$`UNICEF Region`,
    year=a$year, adolescent_rate=a$value, adolescent_censored=a$censored,
    adolescent_source_value=as.character(a$Value), adolescent_lower=as.character(a$Lower),
    adolescent_upper=as.character(a$Upper))
  j <- match(paste(d$iso3,d$year),paste(c$ISO3,c$year))
  d$child_rate <- c$value[j]; d$child_censored <- !is.na(j) & c$censored[j]
  d$child_source_value <- as.character(c$Value[j])
  d$child_lower <- as.character(c$Lower[j]); d$child_upper <- as.character(c$Upper[j])
  d <- d[order(d$iso3,d$year),]; rownames(d)<-NULL
  stopifnot(!anyNA(d$region),all(d$adolescent_rate[!d$adolescent_censored]>0),
            all(d$child_rate[is.finite(d$child_rate)]>0))
  d
}

hiv_stan_data <- function(d, held_out=character()) {
  countries <- unique(d$iso3); regions <- sort(unique(d$region))
  country <- match(d$iso3,countries)
  obs <- which((is.finite(d$child_rate) | d$child_censored) & !d$iso3 %in% held_out)
  make_prev <- function(g,t) {
    prev <- c(0L,seq_len(length(g)-1L)); prev[c(TRUE,diff(g)!=0)]<-0L
    gap <- rep(1,length(g)); hit<-which(prev>0);gap[hit]<-t[hit]-t[prev[hit]]
    list(prev=prev,gap=gap)
  }
  pa<-make_prev(country,d$year);pc<-make_prev(country[obs],d$year[obs])
  ca<-integer(nrow(d)); ca[d$adolescent_censored]<-seq_len(sum(d$adolescent_censored))
  cc<-integer(length(obs)); cc[d$child_censored[obs]]<-seq_len(sum(d$child_censored[obs]))
  # Fixed, modest natural cubic time basis, identical across validation folds.
  basis<-splines::ns((d$year-2012)/12, knots=c(-.5,0,.5),Boundary.knots=c(-1,1))
  basis<-sweep(basis,2,colMeans(basis))
  list(N=nrow(d),C=length(countries),R=length(regions),K=ncol(basis),country=country,
    region=match(d$region[match(countries,d$iso3)],regions),prev_a=pa$prev,gap_a=pa$gap,
    time_basis=unclass(basis),log_a=ifelse(d$adolescent_censored,log(.01),log(d$adolescent_rate)),
    NA_cens=max(ca),cens_a=ca,M=length(obs),child_row=obs,prev_c=pc$prev,gap_c=pc$gap,
    log_c=ifelse(d$child_censored[obs],log(.01),log(d$child_rate[obs])),
    NC_cens=max(cc),cens_c=cc)
}

hiv_fit_diagnostics <- function(fit) {
  parameters<-c("alpha_a","alpha_c","beta_between","beta_within","time_a","time_c",
    "region_a_sd","region_c_sd","country_a_sd","country_c_sd","sigma_a","sigma_c","rho_a","rho_c")
  sm<-rstan::summary(fit,pars=parameters)$summary
  sampler<-rstan::get_sampler_params(fit,inc_warmup=FALSE)
  data.frame(max_rhat=max(sm[,"Rhat"]),min_ess=min(sm[,"n_eff"]),
    divergences=sum(vapply(sampler,function(x)sum(x[,"divergent__"]),numeric(1))),
    treedepth_hits=sum(vapply(sampler,function(x)sum(x[,"treedepth__"]>=12),numeric(1))))
}

hiv_sample <- function(model,sd,seed,chains=4L,iter=2000L) {
  init<-function() list(cens_a_uniform=rep(.5,sd$NA_cens),cens_c_uniform=rep(.5,sd$NC_cens),
    rho_a=.9,rho_c=.9,sigma_a=.5,sigma_c=.5)
  rstan::sampling(model,data=sd,chains=chains,iter=iter,warmup=iter%/%2,
    cores=min(chains,2L),seed=seed,init=init,refresh=100,
    pars=c("mu_a","lp_a","lp_c","country_n"),include=FALSE,
    control=list(adapt_delta=.95,max_treedepth=12))
}

hiv_child_draws <- function(fit,sd,d,n_draws=200L,seed=20260908L) {
  set.seed(seed)
  post<-rstan::extract(fit,pars=c("mu_c","y","sigma_c","rho_c"),permuted=TRUE)
  pick<-sample(seq_along(post$sigma_c),n_draws,replace=FALSE)
  answer<-matrix(NA_real_,n_draws,nrow(d))
  # Draw entire country trajectories conditional on the retained child values.
  # Entirely unobserved countries retain posterior draws of their random effect.
  for (country in unique(d$iso3)) {
    idx<-which(d$iso3==country); obs_j<-which(sd$child_row %in% idx)
    oi<-match(sd$child_row[obs_j],idx); mi<-setdiff(seq_along(idx),oi)
    dist<-abs(outer(d$year[idx],d$year[idx],"-"))
    for (b in seq_len(n_draws)) {
      s<-pick[b]; mu<-post$mu_c[s,idx]
      if(length(oi)) answer[b,idx[oi]]<-post$y[s,obs_j]
      if(length(mi)) {
        V<-post$sigma_c[s]^2 * post$rho_c[s]^dist
        if(length(oi)) {
          W<-V[mi,oi,drop=FALSE] %*% solve(V[oi,oi,drop=FALSE])
          cm<-mu[mi]+drop(W %*% (post$y[s,obs_j]-mu[oi]))
          CV<-V[mi,mi,drop=FALSE]-W %*% V[oi,mi,drop=FALSE]
        } else {cm<-mu[mi]; CV<-V[mi,mi,drop=FALSE]}
        CV<-(CV+t(CV))/2
        answer[b,idx[mi]]<-cm+drop(t(chol(CV+diag(1e-10,length(mi)))) %*% rnorm(length(mi)))
      }
    }
  }
  stopifnot(all(is.finite(answer)),all(answer[,which(d$child_censored & seq_len(nrow(d)) %in% sd$child_row),drop=FALSE] < log(.01)+1e-8))
  answer
}
