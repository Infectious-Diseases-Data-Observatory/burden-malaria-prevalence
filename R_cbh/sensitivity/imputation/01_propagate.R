#!/usr/bin/env Rscript
# Multiple-imputation check for the imputed-covariate sensitivity version.
# For each of M imputed datasets (regional mice imputations, national health-
# expenditure draws and child HIV incidence draws), refit the seven age-band
# models with smoothing parameters fixed at the point-imputation fit and pool
# the PfPR log hazard ratios with Rubin's rules. Pattern follows
# R_cbh/analysis/04_propagate_incidence.R. Aggregate outputs only.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
source("R_cbh/sensitivity/model.R")
source("R_cbh/primary/settings.R")
source("R_cbh/primary/specification.R")
source("R_cbh/reporting/labels.R")
library(mgcv); library(data.table); library(ggplot2)
args <- commandArgs(trailingOnly=TRUE); stopifnot(all(args %in% "--force"))
settings <- cbh_primary_settings("regional_imputed"); stopifnot(isTRUE(settings$imputed))
out <- file.path(settings$out,"multiple_imputation"); private <- file.path(settings$private,"multiple_imputation")
for(p in c(out,private)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
write_csv <- function(d,name) cbh_atomic_csv(d,file.path(out,name))
prepared <- readRDS(settings$data); d0 <- prepared$data; scaling <- prepared$scaling
manifest <- cbh_read_csv(file.path(settings$out,"fit_manifest.csv"))
stopifnot(nrow(manifest)==7L,all(manifest$prepared_data_md5==cbh_file_hash(settings$data)),
  identical(unname(vapply(manifest$model_file,cbh_file_hash,"")),manifest$md5))
regional_imp <- readRDS(file.path(dirname(settings$imputed_overlay),"regional_imputations.rds"))
national_imp <- readRDS(file.path(dirname(settings$imputed_overlay),"national_imputations.rds"))
hiv_draws <- readRDS(sub("child_incidence_country_year_extended[.]csv$","child_incidence_draws_extended.rds",settings$hiv_panel))
hiv_panel <- cbh_read_csv(settings$hiv_panel)
stopifnot(identical(unique(hiv_panel$imputation_signature),hiv_draws$signature))
M <- regional_imp$m; stopifnot(M==ncol(national_imp$log_health_expenditure_pc))
set.seed(20260918); hiv_pick <- sort(sample(seq_len(nrow(hiv_draws$log_incidence)),M))
spec <- cbh_primary_regional_spec(settings); form <- cbh_primary_regional_formula(settings)
ages <- cbh_config()$age_bands$age_band
knot_table <- cbh_read_csv(settings$knots)
# Row lookups, computed once.
rk <- match(paste(d0$survey,d0$regkey),paste(regional_imp$keys$survey,regional_imp$keys$regkey)); stopifnot(!anyNA(rk))
nk <- match(paste(d0$country,d0$entry_year),paste(national_imp$keys$iso3,national_imp$keys$year))
hk <- match(paste(d0$country,d0$entry_year),paste(hiv_draws$keys$iso3,hiv_draws$keys$year)); stopifnot(!anyNA(hk))
hiv_status <- hiv_panel$hiv_incidence_status[match(paste(d0$country,d0$entry_year),paste(hiv_panel$iso3,hiv_panel$year))]
vary_hiv <- hiv_status!="reported_numeric"
imputed_flags <- regional_imp$imputed[rk,,drop=FALSE]
varied <- data.frame(component=c("regional covariates (any)","health expenditure","child HIV incidence (imputed or censored series)"),
  records=c(sum(rowSums(as.matrix(imputed_flags))>0),sum(!is.na(nk)),sum(vary_hiv)))
write_csv(varied,"records_varied_by_component.csv")
impute_dataset <- function(k) {
  d <- d0
  for(v in spec$regional) { ii <- imputed_flags[[v]]; if(any(ii)) d[[v]][ii] <- regional_imp$imputations[[k]][[v]][rk[ii]] }
  ii <- !is.na(nk); d$log_health_expenditure_pc[ii] <- national_imp$log_health_expenditure_pc[nk[ii],k]
  d$log_hiv_incidence[vary_hiv] <- hiv_draws$log_incidence[hiv_pick[k],hk[vary_hiv]]
  for(i in seq_len(nrow(scaling))) { v <- scaling$variable[i]; d[[paste0("z_",v)]] <- (d[[v]]-scaling$mean[i])/scaling$sd[i] }
  stopifnot(!anyNA(d[all.vars(form)])); d
}
code_hash <- vapply(c("R_cbh/sensitivity/imputation/01_propagate.R","R_cbh/primary/specification.R"),cbh_file_hash,"")
curves <- diagnostics <- list()
for(k in seq_len(M)) {
  d <- impute_dataset(k)
  for(i in seq_along(ages)) {
    base <- readRDS(manifest$model_file[i]); f0 <- base$fit
    di <- droplevels(d[d$age_band==ages[i],,drop=FALSE])
    kt <- knot_table[knot_table$age_band==ages[i],]
    knots <- setNames(lapply(c("pfpr_pct","calendar_year"),function(v){z<-kt[kt$variable==v,];z$value[order(z$index)]}),c("pfpr_pct","calendar_year"))
    stopifnot(identical(f0$smooth[[1]]$xp,knots$pfpr_pct))
    id <- sprintf("imputation_%02d_age_%d",k,i); path <- file.path(private,paste0(id,".rds"))
    sig <- cbh_hash(list(manifest$md5[i],regional_imp$seed,national_imp$seed,hiv_draws$signature,hiv_pick[k],k,code_hash,deparse(form)))
    cache <- if(file.exists(path) && !"--force" %in% args) readRDS(path) else NULL
    if(is.null(cache) || !identical(cache$signature,sig)) {
      message("Imputation ",k,"/",M,", age ",ages[i],": ",nrow(di)," records")
      set.seed(settings$seed)
      timing <- system.time(f <- bam(form,data=di,knots=knots,family=binomial(link="cloglog"),method="fREML",
        discrete=TRUE,gamma=settings$gamma,select=FALSE,nthreads=settings$nthreads,gc.level=1,na.action=na.fail,
        sp=f0$sp,coef=coef(f0),control=gam.control(maxit=100)))
      stopifnot(isTRUE(f$converged),all(is.finite(coef(f))),all(is.finite(f$Vp)),f$rank==length(coef(f)))
      cc <- cbh_sensitivity_curves(f,di,id)
      cache <- list(signature=sig,curves=cc,converged=f$converged,iterations=f$iter,elapsed_seconds=unname(timing["elapsed"]),
        coefficients=coef(f),sp=f$sp)
      cbh_atomic_rds(cache,path); rm(f)
    }
    cc <- cache$curves; cc$imputation <- k; curves[[id]] <- cc
    diagnostics[[id]] <- data.frame(imputation=k,age_band=ages[i],converged=cache$converged,iterations=cache$iterations,
      elapsed_seconds=cache$elapsed_seconds)
    rm(base,f0,di); gc(FALSE)
  }
  write_csv(do.call(rbind,diagnostics),"fit_diagnostics.csv")
  rm(d); gc(FALSE)
}
cv <- rbindlist(curves)
# Rubin pooling at every grid point: within = mean variance, between = variance of estimates.
pooled <- cv[,.(estimate=mean(log_hazard_ratio),within=mean(standard_error^2),between=var(log_hazard_ratio),
  pfpr_p025=pfpr_p025[1],pfpr_p975=pfpr_p975[1],within_central_support=within_central_support[1]),by=.(age_band,pfpr_pct)]
pooled[,total:=within+(1+1/M)*between]
pooled[,lambda:=ifelse(total>0,(1+1/M)*between/total,0)]
pooled[,df:=ifelse(lambda>0,(M-1)/lambda^2,Inf)]
pooled[,t:=ifelse(is.finite(df),qt(.975,df),qnorm(.975))]
pooled[,`:=`(standard_error=sqrt(total),lower_95=estimate-t*sqrt(total),upper_95=estimate+t*sqrt(total))]
pooled[,between_share:=ifelse(total>0,(1+1/M)*between/total,0)]
write_csv(as.data.frame(pooled),"pooled_pfpr_curves.csv")
point <- cbh_read_csv(file.path(settings$out,"pfpr_curves.csv")); point <- point[point$series=="map_full",]
contrast <- function(x,label,sign) {
  z <- pooled[pfpr_pct==x]; p <- point[point$pfpr_pct==x,]; p <- p[match(z$age_band,p$age_band),]
  data.frame(age_band=z$age_band,contrast=label,
    point_hr=exp(sign*p$log_hazard_ratio),point_lower_95=exp(sign*p$log_hazard_ratio-1.96*p$standard_error),point_upper_95=exp(sign*p$log_hazard_ratio+1.96*p$standard_error),
    pooled_hr=exp(sign*z$estimate),pooled_lower_95=exp(sign*z$estimate-z$t*z$standard_error),pooled_upper_95=exp(sign*z$estimate+z$t*z$standard_error),
    between_imputation_share=z$between_share,pooled_df=z$df)
}
contrasts <- rbind(contrast(40,"40% to 20%",-1),contrast(0,"20% to 0%",1))
write_csv(contrasts,"pooled_contrasts.csv")
# Overlay: point-imputation fit with its conditional interval against the pooled MI interval.
pl <- as.data.frame(pooled[within_central_support==TRUE]); pl$age_label <- factor(pl$age_band,levels=ages,labels=ifelse(ages=="<1","<1 month",paste(ages,"months")))
pt <- point[point$within_central_support,]; pt$age_label <- factor(pt$age_band,levels=ages,labels=levels(pl$age_label))
p <- ggplot()+geom_hline(yintercept=0,colour="grey60",linewidth=.4)+geom_vline(xintercept=20,colour="grey80",linewidth=.4)+
  geom_ribbon(data=pl,aes(pfpr_pct,ymin=lower_95,ymax=upper_95),fill="#C34D26",alpha=.18)+
  geom_ribbon(data=pt,aes(pfpr_pct,ymin=lower_95,ymax=upper_95),fill="#215E91",alpha=.18)+
  geom_line(data=pl,aes(pfpr_pct,estimate,colour="Multiple imputation (pooled, M = 10)"),linewidth=1)+
  geom_line(data=pt,aes(pfpr_pct,log_hazard_ratio,colour="Point imputation"),linewidth=1,linetype=2)+
  facet_wrap(~age_label,ncol=4)+scale_colour_manual(values=c(`Point imputation`="#215E91",`Multiple imputation (pooled, M = 10)`="#C34D26"),name=NULL)+
  scale_x_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100))+
  labs(x="PfPR[2–10] (%)",y="Log hazard ratio\nrelative to PfPR = 20%")+
  theme_minimal(base_size=18)+theme(panel.grid.minor=element_blank(),legend.position="bottom",strip.text=element_text(face="bold"))
ggsave(file.path(out,"pooled_vs_point_curves.png"),p,width=14,height=8,dpi=220,device=ragg::agg_png,bg="white")
fmt <- function(x) sprintf("%.2f (%.2f–%.2f)",x[[1]],x[[2]],x[[3]])
rows <- unlist(lapply(c("40% to 20%","20% to 0%"),function(c) { z <- contrasts[contrasts$contrast==c,]
  c(paste0("### PfPR ",c),"","| Age (months) | Point imputation HR (95% interval) | Pooled MI HR (95% interval) | Between-imputation share of variance |","|---|---:|---:|---:|",
    sprintf("| %s | %s | %s | %.1f%% |",z$age_band,fmt(z[c("point_hr","point_lower_95","point_upper_95")]),
      fmt(z[c("pooled_hr","pooled_lower_95","pooled_upper_95")]),100*z$between_imputation_share),"") }))
writeLines(c("# Multiple-imputation check: imputed-covariate sensitivity","",
  sprintf("The point fit in `%s` uses one imputed dataset (mean of the mice imputations, GAM fitted means and posterior-median HIV incidence). This check refits all seven age-band models in each of %d imputed datasets with the smoothing parameters fixed at the point fit, and pools the PfPR log hazard ratios with Rubin's rules (within-imputation variance from the conditional covariance, between-imputation variance across the %d fits; Barnard–Rubin degrees of freedom).",settings$id,M,M),"",
  "Varied components per imputation: the mice imputations of the regional covariates, the health-expenditure GAM draws, and one posterior draw of the child HIV incidence series for every country-year whose value is imputed or censored (Liberia and Sao Tome and Principe among them). Observed values never change. Records affected:","",
  sprintf("- %s: %s records",varied$component,format(varied$records,big.mark=",")),"",
  "![Pooled versus point curves](pooled_vs_point_curves.png)","",rows,
  "The between-imputation share is (1 + 1/M) B / T. A small share means the covariate imputation adds little to the sampling uncertainty already in the point fit. Smoothing parameters are fixed, so smoothing uncertainty is not part of either interval. Survey design and exposure uncertainty remain unpropagated.","",
  "Reproduce: `Rscript R_cbh/sensitivity/imputation/01_propagate.R` after `Rscript R_cbh/primary/run_regional.R --imputed`."),
  file.path(out,"REPORT.md"))
paths <- c(settings$data,manifest$model_file,file.path(dirname(settings$imputed_overlay),c("regional_imputations.rds","national_imputations.rds")),
  settings$hiv_panel,"R_cbh/sensitivity/imputation/01_propagate.R")
write_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),"provenance.csv")
print(contrasts[,c("age_band","contrast","point_hr","pooled_hr","pooled_lower_95","pooled_upper_95","between_imputation_share")])
message("Multiple-imputation check complete: ",out)
