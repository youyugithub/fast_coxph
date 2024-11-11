# fast_coxph
```
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List fast_coxph(
    arma::vec time_,
    arma::uvec event_,
    arma::mat xx_){
  
  arma::uvec idx_sort=arma::sort_index(time_);
  
  arma::vec time=time_.elem(idx_sort);
  arma::uvec event=event_.elem(idx_sort);
  arma::mat xx=xx_.rows(idx_sort);
  
  List result_temp(10);
  
  int nind=xx.n_rows;
  int ndim=xx.n_cols;
  arma::vec SS(ndim);
  arma::mat II(ndim,ndim);
  double denominator;
  arma::vec numerator1(ndim);
  arma::mat numerator2(ndim,ndim);
  int iter,ii,ll,ttidx;
  
  arma::vec beta(ndim,arma::fill::zeros);
  arma::vec beta_old(ndim,arma::fill::zeros);
  
  arma::vec expr(nind);
  arma::field<arma::vec> exprx(nind);
  arma::field<arma::mat> exprx2(nind);
  arma::field<arma::vec> F_xx(nind);
  for(ii=0;ii<nind;ii++){
    F_xx(ii)=xx.row(ii).t();
  }
  
  arma::vec unique_time=arma::unique(time.elem(arma::find(event)));
  unique_time=arma::sort(unique_time);
  int n_unique_time=unique_time.n_elem;
  
  arma::uvec ttidx2ii_begin(n_unique_time);
  arma::uvec ttidx2ii_end(n_unique_time);
  arma::uvec temp_uvec;
  for(ttidx=0;ttidx<n_unique_time-1;ttidx++){
    temp_uvec=arma::find(time>=unique_time(ttidx)&&time<unique_time(ttidx+1));
    ttidx2ii_begin(ttidx)=arma::min(temp_uvec);
    ttidx2ii_end(ttidx)=arma::max(temp_uvec);
  }
  temp_uvec=arma::find(time>=unique_time(ttidx));
  ttidx2ii_begin(n_unique_time-1)=arma::min(temp_uvec);
  ttidx2ii_end(n_unique_time-1)=arma::max(temp_uvec);
  
  arma::field<arma::vec> F_SSx(n_unique_time);
  arma::vec SSx(ndim);
  arma::vec nevent(n_unique_time,arma::fill::zeros);
  for(ttidx=0;ttidx<n_unique_time;ttidx++){
    SSx.fill(0.0);
    for(ii=ttidx2ii_begin(ttidx);ii<=ttidx2ii_end(ttidx);ii++){
      if(!event(ii))continue;
      nevent(ttidx)+=1.0;
      SSx+=F_xx(ii);
    }
    F_SSx(ttidx)=SSx;
  }
  
  for(iter=0;iter<300;iter++){
    SS.fill(0.0);
    II.fill(0.0);
    denominator=0.0;
    numerator1.fill(0.0);
    numerator2.fill(0.0);
    expr=arma::exp(xx*beta);
    for(ii=0;ii<nind;ii++){
      exprx(ii)=expr(ii)*F_xx(ii);
      exprx2(ii)=expr(ii)*F_xx(ii)*F_xx(ii).t();
    }
    for(ttidx=n_unique_time-1;ttidx>=0;ttidx--){
      for(ii=ttidx2ii_begin(ttidx);ii<=ttidx2ii_end(ttidx);ii++){
        denominator=denominator+expr(ii);
        numerator1=numerator1+exprx(ii);
        numerator2=numerator2+exprx2(ii);
      }
      SS=SS+F_SSx(ttidx)-nevent(ttidx)*(numerator1/denominator);
      II=II+nevent(ttidx)*(numerator2/denominator-
        (numerator1/denominator)*(numerator1/denominator).t());
    }
    beta_old=beta;
    beta=beta+arma::solve(II,SS);
    if(arma::sum(arma::abs(beta-beta_old))<1e-8)break;
  }
  
  expr=arma::exp(xx*beta);
  arma::vec hazard0(n_unique_time);
  arma::vec hazard0_pulled(n_unique_time);
  denominator=0.0;
  for(ttidx=n_unique_time-1;ttidx>=0;ttidx--){
    for(ii=ttidx2ii_begin(ttidx);ii<=ttidx2ii_end(ttidx);ii++){
      denominator+=expr(ii);
    }
    hazard0(ttidx)=1.0/denominator;
    hazard0_pulled(ttidx)=nevent(ttidx)/denominator;
  }

  // arma::mat IIbb(ndim,ndim);
  // arma::mat IIbh(ndim,n_unique_time);
  // arma::vec IIhh(n_unique_time);
  // IIhh=nevent%sum_expr%sum_expr;
  // 
  // IIbb.fill(0.0);
  // for(ttidx=0;ttidx<n_unique_time;ttidx++){
  //   for(ii=0;ii<nind;ii++){
  //     if(time(ii)<unique_time(ttidx))continue;
  //     IIbb=IIbb+hazard0(ttidx)*exprx2(ii);
  //   }
  // }
  // 
  // IIbh.fill(0.0);
  // for(ttidx=0;ttidx<n_unique_time;ttidx++){
  //   for(ii=0;ii<nind;ii++){
  //     if(time(ii)<unique_time(ttidx))continue;
  //     IIbh.col(ttidx)=IIbh.col(ttidx)-nevent(ttidx)*exprx(ii);
  //   }
  // }
   
  List result(10);
  // result(0)=beta;
  // result(1)=II;
  // result(2)=unique_time;
  // result(3)=hazard0_pulled;
  // result(4)=nevent;
  // result(5)=IIhh;
  // result(6)=IIbb;
  // result(7)=IIbh;
  result(0)=beta;
  result(1)=hazard0_pulled;
  result(2)=unique_time;
  return(result);
}
```

# fast_coxph with information
```
#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List fast_coxph2(
  arma::vec time_,
  arma::uvec event_,
  arma::mat xx_){
  
  arma::uvec idx_sort=arma::sort_index(time_);
  
  arma::vec time=time_.elem(idx_sort);
  arma::uvec event=event_.elem(idx_sort);
  arma::mat xx=xx_.rows(idx_sort);
  
  List result_temp(10);
  
  int nind=xx.n_rows;
  int ndim=xx.n_cols;
  arma::vec SS(ndim);
  arma::mat II(ndim,ndim);
  double denominator;
  arma::vec numerator1(ndim);
  arma::mat numerator2(ndim,ndim);
  int iter,ii,ll,ttidx;
  
  arma::vec beta(ndim,arma::fill::zeros);
  arma::vec beta_old(ndim,arma::fill::zeros);
  
  arma::vec expr(nind);
  arma::field<arma::vec> exprx(nind);
  arma::field<arma::mat> exprx2(nind);
  arma::field<arma::vec> F_xx(nind);
  for(ii=0;ii<nind;ii++){
    F_xx(ii)=xx.row(ii).t();
  }
  
  arma::vec unique_time=arma::unique(time.elem(arma::find(event)));
  unique_time=arma::sort(unique_time);
  int n_unique_time=unique_time.n_elem;
  
  arma::uvec ttidx2ii_begin(n_unique_time);
  arma::uvec ttidx2ii_end(n_unique_time);
  arma::uvec temp_uvec;
  for(ttidx=0;ttidx<n_unique_time-1;ttidx++){
    temp_uvec=arma::find(time>=unique_time(ttidx)&&time<unique_time(ttidx+1));
    ttidx2ii_begin(ttidx)=arma::min(temp_uvec);
    ttidx2ii_end(ttidx)=arma::max(temp_uvec);
  }
  temp_uvec=arma::find(time>=unique_time(ttidx));
  ttidx2ii_begin(n_unique_time-1)=arma::min(temp_uvec);
  ttidx2ii_end(n_unique_time-1)=arma::max(temp_uvec);
  
  arma::field<arma::vec> F_SSx(n_unique_time);
  arma::vec SSx(ndim);
  arma::vec nevent(n_unique_time,arma::fill::zeros);
  for(ttidx=0;ttidx<n_unique_time;ttidx++){
    SSx.fill(0.0);
    for(ii=ttidx2ii_begin(ttidx);ii<=ttidx2ii_end(ttidx);ii++){
      if(!event(ii))continue;
      nevent(ttidx)+=1.0;
      SSx+=F_xx(ii);
    }
    F_SSx(ttidx)=SSx;
  }
  
  for(iter=0;iter<300;iter++){
    SS.fill(0.0);
    II.fill(0.0);
    denominator=0.0;
    numerator1.fill(0.0);
    numerator2.fill(0.0);
    expr=arma::exp(xx*beta);
    for(ii=0;ii<nind;ii++){
      exprx(ii)=expr(ii)*F_xx(ii);
      exprx2(ii)=expr(ii)*F_xx(ii)*F_xx(ii).t();
    }
    for(ttidx=n_unique_time-1;ttidx>=0;ttidx--){
      for(ii=ttidx2ii_begin(ttidx);ii<=ttidx2ii_end(ttidx);ii++){
        denominator=denominator+expr(ii);
        numerator1=numerator1+exprx(ii);
        numerator2=numerator2+exprx2(ii);
      }
      SS=SS+F_SSx(ttidx)-nevent(ttidx)*(numerator1/denominator);
      II=II+nevent(ttidx)*(numerator2/denominator-
                             (numerator1/denominator)*(numerator1/denominator).t());
    }
    beta_old=beta;
    beta=beta+arma::solve(II,SS);
    if(arma::sum(arma::abs(beta-beta_old))<1e-8)break;
  }
  
  expr=arma::exp(xx*beta);
  arma::vec hazard0(n_unique_time);
  arma::vec hazard0_pulled(n_unique_time);
  denominator=0.0;
  arma::vec sum_expr(n_unique_time);
  for(ttidx=n_unique_time-1;ttidx>=0;ttidx--){
    for(ii=ttidx2ii_begin(ttidx);ii<=ttidx2ii_end(ttidx);ii++){
      denominator+=expr(ii);
    }
    sum_expr(ttidx)=denominator;
    hazard0(ttidx)=1.0/denominator;
    hazard0_pulled(ttidx)=nevent(ttidx)/denominator;
  }
  
  
  arma::mat IIbb(ndim,ndim);
  arma::mat IIbh(ndim,n_unique_time);
  arma::vec IIhh(n_unique_time);
  IIhh=nevent%sum_expr%sum_expr;

  IIbb.fill(0.0);
  for(ttidx=0;ttidx<n_unique_time;ttidx++){
    for(ii=0;ii<nind;ii++){
      if(time(ii)<unique_time(ttidx))continue;
      IIbb=IIbb+hazard0(ttidx)*exprx2(ii);
    }
  }

  IIbh.fill(0.0);
  for(ttidx=0;ttidx<n_unique_time;ttidx++){
    for(ii=0;ii<nind;ii++){
      if(time(ii)<unique_time(ttidx))continue;
      IIbh.col(ttidx)=IIbh.col(ttidx)-nevent(ttidx)*exprx(ii);
    }
  }

  List result(10);
  result(0)=beta;
  result(1)=II;
  result(2)=unique_time;
  result(3)=hazard0_pulled;
  
  result(5)=IIhh;
  result(6)=IIbb;
  result(7)=IIbh;
  
  return(result);
}
```
