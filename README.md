---
title: "Measuring Volatility Spillovers and Connectedness"
output: 
  html_document: 
    keep_md: yes
    code_folding: hide
    mathjax: default
  html_notebook:
    code_folding: hide
    mathjax: default
---

### NOTE:
The dashboard which this vignette supports can be found here: (insert link)

## Secton 1: Introduction

This vignette combines and replicates the work done in [Dieblod and Yilmaz 2009](http://www.ssc.upenn.edu/~fdiebold/papers/paper75/DY2final.pdf), [DY 2010](http://www.ssc.upenn.edu/~fdiebold/papers/paper99/DirectionofSpillovers_Mar10.pdf), and [DY 2011](http://www.ssc.upenn.edu/~fdiebold/papers/paper106/dieboldandyilmaz2011.pdf), as well as scaling the methods for quick application and replication. The DY papers describe how to quantify in an index how much error in a forecast can be attributed to shocks to its variables. These "spillovers" are found using standard time-series techniques for finding an impulse-response function and forecast error variance decompositions. In addition, DY 2011 show how the effects of these shocks can be represented as a network when one relaxes the assumptions that vector autoregressions have normally distributed, i.i.d. error terms and that shocks to each variable must be orthongonal to the other shocks. In the spirit of DY 2011, stocks of large banks and AIG will be the subjects of analysis. Like all the DY publications, the analysis will focus on connections and spillovers of volatility.

This vignette will proceed as follows: Section 2) describe data source and summarize data; Section 3) to motivate the concepts of spillovers and connectedness a full sample analysis of the data will be done; Section 4) using a rolling window over the time-series data, apply the same concepts from the full sample analysis to derive index levels of spillover/conncectedness and network behavior over time.

## Section 2: Data

The data is pulled from Yahoo! Finance using the getSymbols() function from the quantmod package, which includes Open, Close, High, and Low prices for a user specified set of stocks/indexes over a chosen period of time measured in days. The data is similar to that in [DY 2010](http://www.ssc.upenn.edu/~fdiebold/papers/paper99/DirectionofSpillovers_Mar10.pdf), so measures of stock return volatility used in the analysis and app are the same as in that publication. First daily variance for stock $i$ at time $t$ is estimated using high and low prices: 

\begin{equation}
  \label{eq:variance}
  \sigma^2_{it} = 0.361[ln(P^{max}_{it}) - ln(P^{min}_{it}))]^2 
\end{equation} 

Since voltatilities are skewed, it is common practice to use log-volatilites which closely approximate a normally distribution. However, to control for instances when volatility is zero, $sinh^{-1}$ is used instead of taking the log ($sinh^{-1}(x) = log(2x)$).

\begin{equation} 
  \sigma_{it} = sinh^{-1}(\sqrt{252*\sigma^2_{it}})
  \label{eq:stdev}
\end{equation}

Using the quantmod getSymbols function, read in the set of stock prices for a period of time. Then calculate their daily volatilities.


```
##  [1] "AIG"  "BAC"  "BK"   "BNS"  "C"    "CS"   "DB"   "GS"   "HSBC" "JPM" 
## [11] "MS"   "TD"   "WFC"
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Name </th>
   <th style="text-align:right;"> count </th>
   <th style="text-align:right;"> mean </th>
   <th style="text-align:right;"> med </th>
   <th style="text-align:right;"> sd </th>
   <th style="text-align:right;"> min </th>
   <th style="text-align:right;"> max </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> BK </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 25.36850 </td>
   <td style="text-align:right;"> 17.97149 </td>
   <td style="text-align:right;"> 26.32555 </td>
   <td style="text-align:right;"> 3.623530 </td>
   <td style="text-align:right;"> 472.3500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JPM </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 25.17900 </td>
   <td style="text-align:right;"> 18.16364 </td>
   <td style="text-align:right;"> 22.94510 </td>
   <td style="text-align:right;"> 4.099518 </td>
   <td style="text-align:right;"> 240.1290 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GS </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 25.28016 </td>
   <td style="text-align:right;"> 19.02999 </td>
   <td style="text-align:right;"> 22.28633 </td>
   <td style="text-align:right;"> 4.768647 </td>
   <td style="text-align:right;"> 314.3188 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 31.09633 </td>
   <td style="text-align:right;"> 20.49637 </td>
   <td style="text-align:right;"> 36.07751 </td>
   <td style="text-align:right;"> 4.012153 </td>
   <td style="text-align:right;"> 567.5506 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MS </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 32.82048 </td>
   <td style="text-align:right;"> 23.69845 </td>
   <td style="text-align:right;"> 35.86672 </td>
   <td style="text-align:right;"> 4.616174 </td>
   <td style="text-align:right;"> 713.4594 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BNS </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 18.19439 </td>
   <td style="text-align:right;"> 14.28432 </td>
   <td style="text-align:right;"> 14.65502 </td>
   <td style="text-align:right;"> 3.388271 </td>
   <td style="text-align:right;"> 193.0767 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AIG </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 34.42879 </td>
   <td style="text-align:right;"> 20.38823 </td>
   <td style="text-align:right;"> 52.74940 </td>
   <td style="text-align:right;"> 3.021434 </td>
   <td style="text-align:right;"> 1344.8577 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TD </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 17.53805 </td>
   <td style="text-align:right;"> 13.88249 </td>
   <td style="text-align:right;"> 13.39231 </td>
   <td style="text-align:right;"> 2.498563 </td>
   <td style="text-align:right;"> 164.3744 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HSBC </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 13.78878 </td>
   <td style="text-align:right;"> 10.52876 </td>
   <td style="text-align:right;"> 11.24531 </td>
   <td style="text-align:right;"> 2.627882 </td>
   <td style="text-align:right;"> 128.7405 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DB </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 22.18650 </td>
   <td style="text-align:right;"> 16.88544 </td>
   <td style="text-align:right;"> 18.17595 </td>
   <td style="text-align:right;"> 3.667744 </td>
   <td style="text-align:right;"> 194.1361 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BAC </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 29.70548 </td>
   <td style="text-align:right;"> 20.62477 </td>
   <td style="text-align:right;"> 31.35738 </td>
   <td style="text-align:right;"> 3.742814 </td>
   <td style="text-align:right;"> 458.1310 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CS </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 20.35036 </td>
   <td style="text-align:right;"> 15.31075 </td>
   <td style="text-align:right;"> 17.76882 </td>
   <td style="text-align:right;"> 3.588789 </td>
   <td style="text-align:right;"> 205.2473 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WFC </td>
   <td style="text-align:right;"> 2768 </td>
   <td style="text-align:right;"> 25.54616 </td>
   <td style="text-align:right;"> 16.98521 </td>
   <td style="text-align:right;"> 26.87410 </td>
   <td style="text-align:right;"> 4.068103 </td>
   <td style="text-align:right;"> 245.8170 </td>
  </tr>
</tbody>
</table>


```r
chart.TimeSeries(vol.data,lwd=2,auto.grid=F,ylab="Annualized Log Volatility",xlab="Time",
                     main="Log Volatility",lty=1,
                     legend.loc="topright")
```

![](README_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


## Section 3: Full Sample Analysis
## Section 3.1: Basic Time-Series: VAR(p), MA(q), Impulse-Response, FEVD

In time-series analysis (and most regressions generally) the most interesting information is found in the distribution of the error terms, especially the variance. In a multivariate model the distibution of error terms embeds the assumptions made regarding correlations of dependent variables. In financial econometrics forecasting stock price returns entails two components - systematic and idiosyncratic shown here:

\begin{align}
  \sigma_{it} = \beta \sigma_{M,t-1} + \epsilon_{it}
  \label{eq:beta}
\end{align}

Where $R_{it}$ is the return on stock $i$ at time $t$ and $\Sigma_M$ is the market return (i.e. the return on the S&P 500, Dow Jones, etc.), $\beta$ is a measure of the systematic portion of a stock's return, or the strength/measure of its relationship with the broader market, and $\epsilon$ is the idiosyncratic portion. The interesting question is how responsive is the future volatility $\sigma_{i,t+1}$, at time $t+1$, to a unit shock of $\epsilon_{it}$? And in a system of several different returns, how does a shock to the idiosyncratic term $\epsilon_{jt}$ for variable $j$ at time $t$ effect $\sigma_{i,t+1}$ for variable $i$ at time $t+1$? In other words, we are interested in finding whether and how idiosyncratic risks of one stock effcts others. These effects are known as spillovers or connectedness.

In order to derive a measure of these spillovers, the calculations will proceed in four steps: 1) Estimate a VAR for a set of $N$ variables of $\sigma_i$ for $i = 1,..,N$ variables which will replace $\sigma_M$; 2) convert the VAR into a moving average representation of order q (MA(q)); 3) derive the impulse-response function for each variable using the MA(q) for each variable $j$ on variable $i$; 4) calculate FEVD matrix for the entire set of variables using the impulse-response function. For textbook treatment of how this is done algebraically, see [Zivot](https://faculty.washington.edu/ezivot/econ584/notes/varModels.pdf), [Cochrane Chapter 7](http://econ.lse.ac.uk/staff/wdenhaan/teach/cochrane.pdf), and [Pesaran and Shin 1997](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.153.5596&rep=rep1&type=pdf) Section 2 for estimating both correlated an uncorrelated shocks.

## Section 3.2: Full Sample Spillover Table
### Step 1. Vector Autoregression (VAR)

Using the VAR and SVAR functions from the [vars](https://cran.r-project.org/web/packages/vars/vars.pdf) package using the set of volatilities estimate a vector autoregression for the volatilities. The SVAR function estimates a Structural VAR which imposes a lower triangular matrix of coefficients (with 1s along the diagonal). This is necessary in order not to violate that assumption that the error terms are normally distributed and i.i.d., which will matter when impulse response functions are calculated for orthogonal and generalized impulses. Thus, the VAR is "non-orthogonal" and the SVAR is "orthogonal". From these models we will need the residual errors (i.e. the idiosyncratic portion) for the calculations that follow.

```r
vol_var = VAR(vol.data,p=3,type="none")
amat <- diag(ncol(vol.data))
amat[lower.tri(amat)] <- NA
vol_svar = SVAR(vol_var,Amat = amat,estmethod = "direct")
### extract residuals of the VAR
res_t <- residuals(vol_var)
svar_ecov <- vol_svar$Sigma.U
```

### Step 2. Moving Average Representation
Next, to estimate the MA(p) representation of the VAR, input the VAR to the Phi function, specifying how many lags (p) we want - in this case 10. Out of this estimate we will need the MA coefficients, along with the residual errors from the VAR, for the next steps.

```r
MA_lag <- 10
theta_temp <- Phi(vol_var,nstep = MA_lag)
svar_theta_temp <- Phi(vol_svar,nstep = MA_lag)
### extract MA coefficients
theta.list <- alply(theta_temp,3)
svar_theta_list <- alply(svar_theta_temp,3)
```

### Step 3 Impulse Response Function (IRF)

The function fevd.matrix implements the IRF and FEVD calculations, taking the VAR residuals and MA coefficients as inputs. For the IRF first calculate the covariance of the residuals, then find the lower-triangular Cholesky-decomposition of the residuals' covariance matrix. The values of this matrix are the "shocks" to the system. 

### Step 4 Forecast Error Variance Decomposition (FEVD)

Second, mulitply each lag of MA cofficients by the Cholesky matrix. The FEVD is estimated by picking out the MA coefficients for variable $j$ and multiplying those by the shocks from variable $i$ as shown in the formula below. The results are shown in the following table.

\begin{equation}
  \theta^o_{ij} = \frac{\sum^n_{l=0} (e'_i A_l P e_j)^2}{\sum^n_{l=0} (e'_i A_l \Sigma A'_l e_i)}
  \label{eq:fevd}
\end{equation}

$\theta^o_{ij}$ is the orthogonal FEVD, $A_l$ is the MA coefficient matrix at lag $l$ out of $n$ lags, $\Sigma$ is the residual covariance matrix, $P$ is the lower triangular matrix, and $e_i$ and $e_j$ are basis vectors with unity at index $i$ and $j$, respectively.

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> BK </th>
   <th style="text-align:right;"> JPM </th>
   <th style="text-align:right;"> GS </th>
   <th style="text-align:right;"> C </th>
   <th style="text-align:right;"> MS </th>
   <th style="text-align:right;"> BNS </th>
   <th style="text-align:right;"> AIG </th>
   <th style="text-align:right;"> TD </th>
   <th style="text-align:right;"> HSBC </th>
   <th style="text-align:right;"> DB </th>
   <th style="text-align:right;"> BAC </th>
   <th style="text-align:right;"> CS </th>
   <th style="text-align:right;"> WFC </th>
   <th style="text-align:right;"> From </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> BK </td>
   <td style="text-align:right;"> 0.58 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JPM </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 0.36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GS </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.67 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.13 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.33 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.48 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 0.52 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MS </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.14 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BNS </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.78 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AIG </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.12 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.21 </td>
   <td style="text-align:right;"> 0.69 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TD </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.85 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HSBC </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.84 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DB </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BAC </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.13 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.49 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 0.51 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CS </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.75 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WFC </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 0.08 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> To </td>
   <td style="text-align:right;"> 0.15 </td>
   <td style="text-align:right;"> 0.44 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 0.13 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.72 </td>
   <td style="text-align:right;"> 0.14 </td>
   <td style="text-align:right;"> 0.74 </td>
   <td style="text-align:right;"> 0.39 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 0.43 </td>
   <td style="text-align:right;"> 0.88 </td>
   <td style="text-align:right;"> 0.37 </td>
  </tr>
</tbody>
</table>

This is the Spillover Table. The row sums (excluding the diagonal terms) measure how much volatility the stock adds to the system. The column sums measure how much volatility each stock receives from the system. Hence, each cell in the table is a measure of how much volatility from column $i$ is given to row $j$. One can estimate the net pairwise spillover by subtracting two spillovers. For instance, the spillover from Goldman Sachs (GS) to Morgan Stanley (MS) is 0.03, and from MS to GS is 0.02. Therefore, the net pairwise spillover from GS to MS is 0.01.

## Section 3.3: Full Sample Connectedness Table

Connectedness requires similar calculations, except that the shocks from IRF are no longer assumed orthogonal. To achieve this, the lower triangular matrix $P$ for $\theta^o_{ij}$ in the numerator is replaced by the covariance matrix of the residuals. The generalized FEVD equation is then,

\begin{equation}
  \theta^g_{ij} = \frac{\sigma^{-1}_{ii} \sum^n_{l=0} (e'_i A_l \Sigma e_j)^2}{\sum^n_{l=0} (e'_i A_l \Sigma A'_l e_i)}
  \label{eq:gfevd}
\end{equation}

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> BK </th>
   <th style="text-align:right;"> JPM </th>
   <th style="text-align:right;"> GS </th>
   <th style="text-align:right;"> C </th>
   <th style="text-align:right;"> MS </th>
   <th style="text-align:right;"> BNS </th>
   <th style="text-align:right;"> AIG </th>
   <th style="text-align:right;"> TD </th>
   <th style="text-align:right;"> HSBC </th>
   <th style="text-align:right;"> DB </th>
   <th style="text-align:right;"> BAC </th>
   <th style="text-align:right;"> CS </th>
   <th style="text-align:right;"> WFC </th>
   <th style="text-align:right;"> From </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> BK </td>
   <td style="text-align:right;"> 12.04 </td>
   <td style="text-align:right;"> 6.06 </td>
   <td style="text-align:right;"> 4.86 </td>
   <td style="text-align:right;"> 11.15 </td>
   <td style="text-align:right;"> 11.66 </td>
   <td style="text-align:right;"> 1.75 </td>
   <td style="text-align:right;"> 27.05 </td>
   <td style="text-align:right;"> 1.76 </td>
   <td style="text-align:right;"> 0.94 </td>
   <td style="text-align:right;"> 3.02 </td>
   <td style="text-align:right;"> 9.56 </td>
   <td style="text-align:right;"> 3.42 </td>
   <td style="text-align:right;"> 6.71 </td>
   <td style="text-align:right;"> 87.96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> JPM </td>
   <td style="text-align:right;"> 6.27 </td>
   <td style="text-align:right;"> 9.66 </td>
   <td style="text-align:right;"> 5.38 </td>
   <td style="text-align:right;"> 16.36 </td>
   <td style="text-align:right;"> 9.74 </td>
   <td style="text-align:right;"> 2.23 </td>
   <td style="text-align:right;"> 18.28 </td>
   <td style="text-align:right;"> 2.18 </td>
   <td style="text-align:right;"> 1.22 </td>
   <td style="text-align:right;"> 3.46 </td>
   <td style="text-align:right;"> 12.89 </td>
   <td style="text-align:right;"> 3.38 </td>
   <td style="text-align:right;"> 8.95 </td>
   <td style="text-align:right;"> 90.34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GS </td>
   <td style="text-align:right;"> 6.62 </td>
   <td style="text-align:right;"> 6.29 </td>
   <td style="text-align:right;"> 10.14 </td>
   <td style="text-align:right;"> 11.35 </td>
   <td style="text-align:right;"> 14.29 </td>
   <td style="text-align:right;"> 2.05 </td>
   <td style="text-align:right;"> 23.93 </td>
   <td style="text-align:right;"> 2.04 </td>
   <td style="text-align:right;"> 1.12 </td>
   <td style="text-align:right;"> 3.69 </td>
   <td style="text-align:right;"> 8.96 </td>
   <td style="text-align:right;"> 3.65 </td>
   <td style="text-align:right;"> 5.86 </td>
   <td style="text-align:right;"> 89.86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C </td>
   <td style="text-align:right;"> 4.06 </td>
   <td style="text-align:right;"> 5.19 </td>
   <td style="text-align:right;"> 3.54 </td>
   <td style="text-align:right;"> 30.65 </td>
   <td style="text-align:right;"> 6.46 </td>
   <td style="text-align:right;"> 1.68 </td>
   <td style="text-align:right;"> 18.58 </td>
   <td style="text-align:right;"> 1.62 </td>
   <td style="text-align:right;"> 0.82 </td>
   <td style="text-align:right;"> 2.67 </td>
   <td style="text-align:right;"> 14.35 </td>
   <td style="text-align:right;"> 2.67 </td>
   <td style="text-align:right;"> 7.71 </td>
   <td style="text-align:right;"> 69.35 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MS </td>
   <td style="text-align:right;"> 6.65 </td>
   <td style="text-align:right;"> 5.15 </td>
   <td style="text-align:right;"> 5.95 </td>
   <td style="text-align:right;"> 9.31 </td>
   <td style="text-align:right;"> 20.68 </td>
   <td style="text-align:right;"> 1.51 </td>
   <td style="text-align:right;"> 28.48 </td>
   <td style="text-align:right;"> 1.46 </td>
   <td style="text-align:right;"> 0.85 </td>
   <td style="text-align:right;"> 3.05 </td>
   <td style="text-align:right;"> 8.44 </td>
   <td style="text-align:right;"> 3.16 </td>
   <td style="text-align:right;"> 5.31 </td>
   <td style="text-align:right;"> 79.32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BNS </td>
   <td style="text-align:right;"> 7.46 </td>
   <td style="text-align:right;"> 6.71 </td>
   <td style="text-align:right;"> 5.65 </td>
   <td style="text-align:right;"> 13.93 </td>
   <td style="text-align:right;"> 10.83 </td>
   <td style="text-align:right;"> 6.34 </td>
   <td style="text-align:right;"> 15.56 </td>
   <td style="text-align:right;"> 4.40 </td>
   <td style="text-align:right;"> 1.35 </td>
   <td style="text-align:right;"> 4.80 </td>
   <td style="text-align:right;"> 11.08 </td>
   <td style="text-align:right;"> 4.60 </td>
   <td style="text-align:right;"> 7.32 </td>
   <td style="text-align:right;"> 93.66 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AIG </td>
   <td style="text-align:right;"> 1.79 </td>
   <td style="text-align:right;"> 1.42 </td>
   <td style="text-align:right;"> 0.99 </td>
   <td style="text-align:right;"> 4.80 </td>
   <td style="text-align:right;"> 1.49 </td>
   <td style="text-align:right;"> 0.57 </td>
   <td style="text-align:right;"> 78.93 </td>
   <td style="text-align:right;"> 0.44 </td>
   <td style="text-align:right;"> 0.36 </td>
   <td style="text-align:right;"> 0.72 </td>
   <td style="text-align:right;"> 4.44 </td>
   <td style="text-align:right;"> 0.69 </td>
   <td style="text-align:right;"> 3.37 </td>
   <td style="text-align:right;"> 21.07 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TD </td>
   <td style="text-align:right;"> 6.88 </td>
   <td style="text-align:right;"> 7.03 </td>
   <td style="text-align:right;"> 5.71 </td>
   <td style="text-align:right;"> 14.71 </td>
   <td style="text-align:right;"> 9.63 </td>
   <td style="text-align:right;"> 4.65 </td>
   <td style="text-align:right;"> 15.87 </td>
   <td style="text-align:right;"> 5.89 </td>
   <td style="text-align:right;"> 1.44 </td>
   <td style="text-align:right;"> 4.76 </td>
   <td style="text-align:right;"> 11.30 </td>
   <td style="text-align:right;"> 4.39 </td>
   <td style="text-align:right;"> 7.73 </td>
   <td style="text-align:right;"> 94.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HSBC </td>
   <td style="text-align:right;"> 5.85 </td>
   <td style="text-align:right;"> 6.13 </td>
   <td style="text-align:right;"> 5.27 </td>
   <td style="text-align:right;"> 15.25 </td>
   <td style="text-align:right;"> 9.18 </td>
   <td style="text-align:right;"> 3.09 </td>
   <td style="text-align:right;"> 19.40 </td>
   <td style="text-align:right;"> 2.65 </td>
   <td style="text-align:right;"> 2.98 </td>
   <td style="text-align:right;"> 5.71 </td>
   <td style="text-align:right;"> 12.13 </td>
   <td style="text-align:right;"> 4.81 </td>
   <td style="text-align:right;"> 7.56 </td>
   <td style="text-align:right;"> 97.02 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DB </td>
   <td style="text-align:right;"> 6.43 </td>
   <td style="text-align:right;"> 5.79 </td>
   <td style="text-align:right;"> 5.69 </td>
   <td style="text-align:right;"> 12.33 </td>
   <td style="text-align:right;"> 11.59 </td>
   <td style="text-align:right;"> 2.79 </td>
   <td style="text-align:right;"> 17.19 </td>
   <td style="text-align:right;"> 2.51 </td>
   <td style="text-align:right;"> 1.80 </td>
   <td style="text-align:right;"> 9.89 </td>
   <td style="text-align:right;"> 12.11 </td>
   <td style="text-align:right;"> 5.65 </td>
   <td style="text-align:right;"> 6.22 </td>
   <td style="text-align:right;"> 90.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BAC </td>
   <td style="text-align:right;"> 4.88 </td>
   <td style="text-align:right;"> 6.04 </td>
   <td style="text-align:right;"> 3.71 </td>
   <td style="text-align:right;"> 20.26 </td>
   <td style="text-align:right;"> 7.79 </td>
   <td style="text-align:right;"> 1.78 </td>
   <td style="text-align:right;"> 14.09 </td>
   <td style="text-align:right;"> 1.60 </td>
   <td style="text-align:right;"> 1.01 </td>
   <td style="text-align:right;"> 3.28 </td>
   <td style="text-align:right;"> 22.54 </td>
   <td style="text-align:right;"> 3.13 </td>
   <td style="text-align:right;"> 9.89 </td>
   <td style="text-align:right;"> 77.46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CS </td>
   <td style="text-align:right;"> 7.55 </td>
   <td style="text-align:right;"> 5.90 </td>
   <td style="text-align:right;"> 5.37 </td>
   <td style="text-align:right;"> 13.00 </td>
   <td style="text-align:right;"> 11.67 </td>
   <td style="text-align:right;"> 2.54 </td>
   <td style="text-align:right;"> 19.33 </td>
   <td style="text-align:right;"> 2.44 </td>
   <td style="text-align:right;"> 1.48 </td>
   <td style="text-align:right;"> 5.80 </td>
   <td style="text-align:right;"> 11.30 </td>
   <td style="text-align:right;"> 7.23 </td>
   <td style="text-align:right;"> 6.40 </td>
   <td style="text-align:right;"> 92.77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WFC </td>
   <td style="text-align:right;"> 5.90 </td>
   <td style="text-align:right;"> 6.41 </td>
   <td style="text-align:right;"> 3.98 </td>
   <td style="text-align:right;"> 18.22 </td>
   <td style="text-align:right;"> 7.44 </td>
   <td style="text-align:right;"> 2.04 </td>
   <td style="text-align:right;"> 18.32 </td>
   <td style="text-align:right;"> 1.87 </td>
   <td style="text-align:right;"> 1.12 </td>
   <td style="text-align:right;"> 2.89 </td>
   <td style="text-align:right;"> 15.17 </td>
   <td style="text-align:right;"> 3.02 </td>
   <td style="text-align:right;"> 13.62 </td>
   <td style="text-align:right;"> 86.38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> To </td>
   <td style="text-align:right;"> 70.32 </td>
   <td style="text-align:right;"> 68.11 </td>
   <td style="text-align:right;"> 56.08 </td>
   <td style="text-align:right;"> 160.68 </td>
   <td style="text-align:right;"> 111.76 </td>
   <td style="text-align:right;"> 26.66 </td>
   <td style="text-align:right;"> 236.11 </td>
   <td style="text-align:right;"> 24.98 </td>
   <td style="text-align:right;"> 13.52 </td>
   <td style="text-align:right;"> 43.86 </td>
   <td style="text-align:right;"> 131.72 </td>
   <td style="text-align:right;"> 42.57 </td>
   <td style="text-align:right;"> 83.04 </td>
   <td style="text-align:right;"> 82.26 </td>
  </tr>
</tbody>
</table>

The interpretation of the table as well as the calculation for net pairwise connectedness is not the same as that of the spillover. In DY 2011 the connectedness table can be interpreted as a network in which each entry represents a weight of directional connections between the nodes. The row sums are the from-degrees, column sums are to-degrees. The total connectedness (bottom-right entry in the table) is the average degree of the network (to or from since the row sums equal column sums). Since generalized FEVD matrices do not have row/columns sums equal to one, the entries in the connectedness have been normalized to that end. 

## Section 3.4: Full Sample FEVD as a Network

DY 2011 demonstrates that variance decompositions define weighted, directed networks. Thus the Connectedness table can be transformed into a network as shown below. (To see how the network evolves over time, use the app here: (insert link))

```
## + 13/13 vertices, named, from 7a96360:
##  [1] AIG  C    BAC  BK   JPM  GS   MS   BNS  TD   HSBC DB   CS   WFC
```

![](README_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## Section 4: Rolling Window - Spillover Vs. Connected

To apply the methods shown above to a rolling window, it is neccessary to only define the size of the rolling window, the size of the increment to move the window, which will determine the how many times to loop through the VAR function. The function vector_autoreg takes in the an .xts data frame, window size, increment, AR lag order, and MA lag order. It also calls on function "make.table" and "normalize" which can be found in the function script on the GitHub page.


```r
vector_autoreg <- function(theData.xts,window=100,iter=10,ar_lag=3,ma_lag=10){
  ### create a sequence to use as cut offs for the rolling window
  ### estimation. iter determines how far to move the window forward
  
  mySeq <- seq(1,nrow(theData.xts)-window-1,by=iter)
  min_iter = min(diff(index(theData.xts[mySeq])),na.rm=T)
  date_names <- as.character(seq.Date(from=min(index(theData.xts)),by=min_iter,length.out = length(mySeq)))
  
  model_list <- vector("list",length=length(date_names))
  names(model_list) <- date_names
  ### initialize a list to hold index values
  fevd.list <- vector("list",length=length(date_names))
  ### names of the list come from the index of .xts data object
  names(fevd.list) <- date_names
  
  spillover.table.list <- vector("list",length=length(date_names))
  names(spillover.table.list) <- date_names
  
  vol.conn.index <- vector("list",length=length(date_names))
  ### names of the list come from the index of .xts data object
  names(vol.conn.index) <- date_names
  
  network.list <- vector("list",length=length(date_names))
  names(network.list) = date_names
  net.df.list <- network.list
  nvars <- ncol(theData.xts)
  
  for(d in 1:length(mySeq)){
    ### estimate VAR(p) and SVAR(p)
    var.vol.temp <- VAR(theData.xts[mySeq[d]:(mySeq[d]+window-1)],p=ar_lag,type="none")
    amat <- diag(nvars)
    amat[lower.tri(amat)] <- NA
    svar.vol.temp <- SVAR(var.vol.temp,estmethod = "direct",Amat=amat)
    
    model_list[[date_names[d]]] <- list(orthogonal=svar.vol.temp,
                                        non_orthogonal = var.vol.temp)
    ### get MA coefficients
    theta.temp <- Phi(var.vol.temp,nstep=ma_lag)
    svar.theta <- Phi(svar.vol.temp,nstep=ma_lag)
    ### convert array to listn.ahead
    theta.list <- alply(theta.temp,3)
    svar.theta.list <- alply(svar.theta,3)
    ################## SPILLOVER #########################
    svar_ecov <-svar.vol.temp$Sigma.U
    error.variance.total <- fevd.matrix(svar_ecov,svar.theta.list)*100
    spillover.table.list[[date_names[d]]] <- make.table(error.variance.total,colnames(theData.xts),date_names[d])
    
    
    Qinv <- t(chol(svar_ecov))
    irf <- llply(theta.list,function(x) x%*%Qinv)
    num = llply(irf,cond.sum)
    ### calculate own variance shares (i.e. do not subtract the diagonal of the matrix)
    den = llply(irf,get.trace)
    sum.num = sum(unlist(num))
    sum.den = sum(unlist(den))
    S = sum.num*100/sum.den
    fevd.list[[date_names[d]]] <- S
    ########################################################
    
    ################## CONNECTEDNESS #######################
    e.cov <- cov(residuals(var.vol.temp))
    ### create directional network table
    Sconn = directional.matrix(residuals(var.vol.temp),theta.list)
    ### to normalize or not to normalize
    D = normalize(Sconn)*100
    ### calculate total connectedness index level
    vol.conn.index[[date_names[d]]] = (sum(apply(D,1,sum)  - diag(D)))/ncol(D)
    ########################################################
    
    ################## Networks #######################
    rownames(D) <- colnames(theData.xts)
    colnames(D) <- rownames(D)
    
    df.list <- make.table(D,rownames(D),date_names[d])
    
    net.df.list[[date_names[d]]] <- df.list
    ### Create network Graph
    
    ### calculate net pairwise directional connectedness, Diebold pg 4.
    ### the calculation t(D) - D is the net pairwise equivalent of 
    ### the gross net = to - from, the calculation just above.
    abs.D = abs(D)
    net.mat = t(abs.D) - abs.D
    net.melt1 <- melt(net.mat)
    net.melt <- net.melt1 %>% filter(value != 0 )

    ### calculate each percentile of the net pairwise connectedness values
    ### and choose only the top 10%
    net.quant <- quantile(net.melt$value,prob=seq(0,1,by=0.01))
    new.net1 <- net.melt[net.melt$value > net.quant[90],]
    new.net <- new.net1[new.net1[,1] != new.net1[,2],]
    
    ### create igraph graph
    net.network <- graph.data.frame(new.net,direct=T)
    ### set graph nodes
    V(net.network)
    ### set edge colors
    E(net.network)$color <- ifelse(E(net.network)$value >= net.quant[99],"black",
                                   ifelse((E(net.network)$value < net.quant[99] & E(net.network)$value >= net.quant[95]),"red",
                                          ifelse((E(net.network)$value < net.quant[95] & E(net.network)$value >= net.quant[90]),"orange","grey")))
    
    ### set node size
    V(net.network)$size <- degree(net.network)/.5
    
    network.list[[date_names[d]]] <- net.network
    #####################################################################################
    remove(var.vol.temp,theta.temp,theta.list,e.cov,Sconn,
           num,den,sum.num,sum.den,S)
  }
  
  list(models = model_list,
       spillover_table = spillover.table.list,
       spill = fevd.list,
       conn=vol.conn.index,
       nets = network.list,
       net_df = net.df.list,
       dates = date_names)
  
}

rolling_var <- vector_autoreg(theData.xts = vol.data)
```
# Section 4.1: Spillover & Connectedness Indexes


```r
############## CONNECTEDNESS ################
vol.conn.index.out <- rolling_var[["conn"]]
### unlist output data
vol.conn.index.df <- data.frame(unlist(vol.conn.index.out))
colnames(vol.conn.index.df) = c("Index")
### .xts object
vol.conn.index.df$Date <- as.POSIXct(rownames(vol.conn.index.df),format="%Y-%m-%d")
vol.conn.index.xts <- xts(vol.conn.index.df[,-2],order.by=vol.conn.index.df[,2])

############## SPILLOVER ###################
fevd.list.out <- rolling_var[["spill"]]

fevd.df <- data.frame(unlist(fevd.list.out))
colnames(fevd.df) = c("Index")
fevd.df$Date <- as.POSIXct(rownames(fevd.df),format="%Y-%m-%d")
fevd.xts <- xts(fevd.df[,-2],order.by=fevd.df[,2])

### compare the two index measures
indice = merge.all(vol.conn.index.xts,fevd.xts)
colnames(indice) = c("Connectedness","Spillover")

chart.TimeSeries(indice,lwd=2,auto.grid=F,ylab="Index",xlab="Time",
                      main="Comparing Spillover Index levels",lty=1,
                      legend.loc="topright")
```

![](README_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

The indexes show the percentage of the 10-day ahead forecast error which can be explained by spillovers/connectedness. Mathematically, for each matrix calculated at each point in time the sum of all the off-diagonal entries are divided by the sum of the entire matrix. Surpisingly, the indexes are very correlated, hence the use of orthogonal shocks yields a similar amount of forecast error as correlated shocks in this instance.

## Section 4.2: Total Directional Connectedness and Spillover

Using the time-series of connectedness and spillover tables it is also possible to see how much volatility each stock gave or received, and which ones were net givers/receivers of volatility over time. The red lines are the Spillover measures and the black lines are Connectedness.

These results show the consequences of the assumptions behind the SVAR and VAR models. During the run up to the 2008 financial crisis and in its wake, AIG underwrote several option contracts for various asset-backed securities with the large banks. However, the SVAR does not pick up on that connection, unlike the "non-orthogonal" VAR, and instead estimates a lower volatility spillover of AIG stock price to the other prices. A likely reason is that SVARs are sensitive to the ordering of variables used to satisfy the normality assumption, and, perhaps, a different ordering would deliver a different result. In addition, these results will change for the SVAR when the parameters to the function vector_autoreg() are altered. 


```r
net.df.list <- rolling_var[["net_df"]]
spill.df.list <- rolling_var[["spillover_table"]]

net.cols <- names(net.df.list)

to.selected.net.df.list <- llply(net.df.list,function(df){df[["To"]]})
from.selected.net.df.list <- llply(net.df.list,function(df){df[["From"]]})
net.selected.net.df.list <- llply(net.df.list,function(df){df[["Net"]]})

to.selected.net.df <- as.data.frame(do.call(rbind,lapply(to.selected.net.df.list, function(x) t(data.frame(x)))))
from.selected.net.df <- as.data.frame(do.call(rbind,lapply(from.selected.net.df.list, function(x) t(data.frame(x)))))
net.selected.net.df <- as.data.frame(do.call(rbind,lapply(net.selected.net.df.list, function(x) t(data.frame(x)))))

to.selected.net.df <- to.selected.net.df %>% mutate(Date = as.POSIXct(net.cols,format="%Y-%m-%d"))
from.selected.net.df <- from.selected.net.df %>% mutate(Date = as.POSIXct(net.cols,format="%Y-%m-%d"))
net.selected.net.df <- net.selected.net.df %>% mutate(Date = as.POSIXct(net.cols,format="%Y-%m-%d"))

spill.cols <- names(spill.df.list)

to.selected.spill.df.list <- llply(spill.df.list,function(df){df[["To"]]})
from.selected.spill.df.list <- llply(spill.df.list,function(df){df[["From"]]})
net.selected.spill.df.list <- llply(spill.df.list,function(df){df[["Net"]]})

to.selected.spill.df <- as.data.frame(do.call(rbind,lapply(to.selected.spill.df.list, function(x) t(data.frame(x)))))
from.selected.spill.df <- as.data.frame(do.call(rbind,lapply(from.selected.spill.df.list, function(x) t(data.frame(x)))))
net.selected.spill.df <- as.data.frame(do.call(rbind,lapply(net.selected.spill.df.list, function(x) t(data.frame(x)))))

to.selected.spill.df <- to.selected.spill.df %>% mutate(Date = as.POSIXct(spill.cols,format="%Y-%m-%d"))
from.selected.spill.df <- from.selected.spill.df %>% mutate(Date = as.POSIXct(spill.cols,format="%Y-%m-%d"))
net.selected.spill.df <- net.selected.spill.df %>% mutate(Date = as.POSIXct(spill.cols,format="%Y-%m-%d"))

to.net.melt <- melt(to.selected.net.df,id="Date",value.name = "Connectedness")
to.spill.melt <- melt(to.selected.spill.df,id="Date",value.name = "Spillover")
to.tbl <- merge(to.net.melt,to.spill.melt,by=c("Date","variable"))

from.net.melt <- melt(from.selected.net.df,id="Date",value.name = "Connectedness")
from.spill.melt <- melt(from.selected.spill.df,id="Date",value.name = "Spillover")
from.tbl <- merge(from.net.melt,from.spill.melt,by=c("Date","variable"))

net.net.melt <- melt(net.selected.net.df,id="Date",value.name = "Connectedness")
net.spill.melt <- melt(net.selected.spill.df,id="Date",value.name = "Spillover")
net.tbl <- merge(net.net.melt,net.spill.melt,by=c("Date","variable"))
```

### Giving


```r
ggplot(to.tbl) + 
  geom_line(aes(x=Date,y=Connectedness),colour="black") +
  geom_line(aes(x=Date,y=Spillover),colour="red") +
  facet_wrap(~variable,scales="free_y") +
  theme(legend.position = "bottom")
```

![](README_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

### Receiving


```r
ggplot(from.tbl) + 
  geom_line(aes(x=Date,y=Connectedness),colour="black") +
  geom_line(aes(x=Date,y=Spillover),colour="red") +
  facet_wrap(~variable,scales="free_y") +
  theme(legend.position = "bottom")
```

![](README_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

### Net


```r
ggplot(net.tbl) + 
  geom_line(aes(x=Date,y=Connectedness),colour="black") +
  geom_line(aes(x=Date,y=Spillover),colour="red") +
  facet_wrap(~variable,scales="free_y") +
  theme(legend.position = "bottom")
```

![](README_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

## Section 4.3: Rolling Window Net Pairwise Spillovers & Connectedness

Using the connectedness measure, it is also possible to compare two stocks and see which one was the net giver or receiver. For instance if we compare the Morgan Stanley (MS) and Goldman Sachs (GS), it appears that Morgan Stanley "gave" volatility to Goldman Sachs for most of the time considered.


```r
sender <- "MS"
receiver <- "GS"

net_stuff <- rolling_var[["net_df"]]
tables <- llply(net_stuff,function(stuff) stuff[["table"]])

pairwise_data <- ldply(tables,function(table){
  sender_index <- which(colnames(table) == sender)
  receiver_index <- which(colnames(table) == receiver)
  
  sender_shock <- table[receiver_index,sender_index]
  receiver_shock <- table[sender_index,receiver_index]
  
  net_shock <- sender_shock - receiver_shock
})
colnames(pairwise_data) <- c("Time","Connectedness")
pairwise_data <- pairwise_data %>% mutate(Time = as.POSIXct(Time,format="%Y-%m-%d"))

ggplot(pairwise_data,aes(x=Time,y=Connectedness,group=1)) + 
  geom_line() +
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  ggtitle(paste(sender,"to", receiver)) + 
  xlab("Time") + 
  ylab("Connectedness")
```

![](README_files/figure-html/unnamed-chunk-14-1.png)<!-- -->
