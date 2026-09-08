#!/usr/bin/env Rscript
# Ten coherent external posterior draws; conditional mortality fits with fixed
# smoothing parameters from the median-incidence model. Other X remain complete cases.
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
stopifnot(requireNamespace("mgcv",quietly=TRUE))
spec<-cbh_trial_spec()
private<-file.path("data/derived_cbh/models",spec$id)
out<-file.path("results/cbh",spec$id,"multiple_imputation")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
prepared<-readRDS(file.path(private,"complete_case_dataset.rds"))
base<-readRDS(file.path(private,"fit.rds"))
stopifnot(identical(prepared$signature,base$signature),isTRUE(base$fit$converged))
external<-readRDS(spec$incidence_draws)
panel<-cbh_read_csv(spec$incidence_panel)
stopifnot(identical(unique(panel$imputation_signature),external$signature))
d<-prepared$data
j<-match(paste(d$country,d$entry_year),paste(external$keys$iso3,external$keys$year))
stopifnot(!anyNA(j),!anyNA(external$log_incidence[,j[!duplicated(j)],drop=FALSE]))
set.seed(20260908)
selected<-sort(sample(seq_len(nrow(external$log_incidence)),10))
writeLines(c("Ten coherent external incidence draws; fixed mortality smoothing parameters; Rubin-style pointwise pooling.",
  paste("Draw indices:",paste(selected,collapse=", ")),
  paste("Base fit signature:",base$signature),paste("Incidence signature:",external$signature)),
  file.path(out,"specification.txt"))
scaling<-prepared$scaling[prepared$scaling$variable=="log_hiv_incidence",]
age<-levels(d$age_band)
templates<-d[match(age,d$age_band),]
curves<-contrasts<-diagnostics<-list()
for(k in seq_along(selected)) {
  draw<-selected[k]
  path<-file.path(private,sprintf("incidence_draw_%03d.rds",draw))
  sig<-cbh_hash(list(base$signature,external$signature,draw,base$fit$sp,
                    cbh_file_hash("R_cbh/analysis/04_propagate_incidence.R")))
  cache<-if(file.exists(path)) readRDS(path) else NULL
  if(is.null(cache)||!identical(cache$signature,sig)) {
    message("Mortality uncertainty fit ",k,"/",length(selected),"; incidence draw ",draw)
    d$log_hiv_incidence<-external$log_incidence[draw,j]
    d$hiv_incidence_per1000<-exp(d$log_hiv_incidence)
    d$z_log_hiv_incidence<-(d$log_hiv_incidence-scaling$mean)/scaling$sd
    cache<-cbh_trial_fit(d,spec,trace=FALSE,start=coef(base$fit),smoothing=base$fit$sp)
    cache$signature<-sig;cache$draw<-draw
    cache$base_input_signature<-base$signature;cache$external_signature<-external$signature
    cbh_atomic_rds(cache,path)
  }
  f<-cache$fit
  stopifnot(isTRUE(f$converged),all(is.finite(coef(f))),all(is.finite(f$Vp)))
  diagnostics[[k]]<-data.frame(draw=draw,converged=f$converged,iterations=f$iter,
                               elapsed_seconds=cache$elapsed_seconds)
  z<-cbh_trial_contrasts(f,d);z$draw<-draw;contrasts[[k]]<-z
  curve<-lapply(seq_along(age),function(i) {
    support<-quantile(d$pfpr_pct[d$age_band==age[i]],c(.025,.975),names=FALSE)
    grid<-sort(unique(c(seq(support[1],support[2],length.out=151),20,40)))
    nd<-templates[rep(i,length(grid)),];ref<-nd
    nd$pfpr_pct<-grid;ref$pfpr_pct<-20
    L<-predict(f,nd,type="lpmatrix",discrete=FALSE)-predict(f,ref,type="lpmatrix",discrete=FALSE)
    vv<-rowSums((L%*%f$Vp)*L)
    stopifnot(min(vv)>-1e-8)
    data.frame(draw=draw,age_band=age[i],pfpr_pct=grid,
      estimate=drop(L%*%coef(f)),variance=pmax(vv,0))
  })
  curves[[k]]<-do.call(rbind,curve)
  cbh_atomic_csv(do.call(rbind,diagnostics),file.path(out,"fit_diagnostics.csv"))
  rm(cache,f);gc(FALSE)
}
pool<-function(x) {
  m<-nrow(x);u<-mean(x$variance);b<-var(x$estimate);v<-u+(1+1/m)*b
  df<-if(b>1e-20) (m-1)*(1+u/((1+1/m)*b))^2 else Inf
  est<-mean(x$estimate);half<-qt(.975,df)*sqrt(v)
  data.frame(estimate=est,within_variance=u,between_variance=b,total_se=sqrt(v),
             df=df,lower_95=est-half,upper_95=est+half,
             monte_carlo_se=sqrt(b/m),imputations=m)
}
cc<-do.call(rbind,curves)
pc<-do.call(rbind,lapply(split(cc,paste(cc$age_band,cc$pfpr_pct)),function(x)
  cbind(x[1,c("age_band","pfpr_pct")],pool(x))))
rownames(pc)<-NULL
cbh_atomic_csv(pc,file.path(out,"pooled_pfpr_curves.csv"))
zz<-do.call(rbind,contrasts)
zz$estimate<-log(zz$hazard_ratio)
zz$variance<-((log(zz$upper_95)-log(zz$lower_95))/(2*1.96))^2
pooled<-do.call(rbind,lapply(split(zz,zz$age_band),function(x)cbind(age_band=x$age_band[1],pool(x))))
pooled$hazard_ratio<-exp(pooled$estimate)
pooled$hr_lower_95<-exp(pooled$lower_95);pooled$hr_upper_95<-exp(pooled$upper_95)
pooled<-pooled[match(age,pooled$age_band),];rownames(pooled)<-NULL
cbh_atomic_csv(pooled,file.path(out,"pooled_pfpr_40_to_20.csv"))
cbh_atomic_csv(zz,file.path(out,"individual_pfpr_40_to_20.csv"))
library(ggplot2)
pc$age_band<-factor(pc$age_band,levels=age,
  labels=ifelse(age=="<1","<1 month",paste0(age," months")))
p<-ggplot(pc,aes(pfpr_pct,estimate))+
  geom_hline(yintercept=0,linetype="dashed",colour="grey60")+
  geom_vline(xintercept=20,colour="grey85")+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95),fill="#176B87",alpha=.18)+
  geom_line(colour="#176B87",linewidth=.8)+facet_wrap(~age_band,ncol=4)+
  labs(title="PfPR splines adjusted for child HIV incidence",
    subtitle="Ten country-year incidence imputations | Central 95% of observed PfPR values",
    x="PfPR at band entry (%)",y="Log mortality hazard ratio relative to 20% PfPR",
    caption="Pointwise intervals pool conditional model uncertainty and between-imputation variation.\nMortality smoothing parameters fixed at the median-incidence fit; survey-design and source-estimate uncertainty are not included.")+
  theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0),
    strip.text=element_text(face="bold",hjust=0))
ggsave(file.path(out,"pooled_pfpr_splines.png"),p,device=ragg::agg_png,width=14,height=8,dpi=160,bg="white")
print(pooled)
