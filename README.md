Measuring Volatility Spillovers and Connectedness
================
Rhamey Sayed

### NOTE:

The dashboard which this vignette supports can be found here: (insert link)

Secton 1: Introduction
----------------------

This vignette combines and replicates the work done in [Dieblod and Yilmaz 2009](http://www.ssc.upenn.edu/~fdiebold/papers/paper75/DY2final.pdf), [DY 2010](http://www.ssc.upenn.edu/~fdiebold/papers/paper99/DirectionofSpillovers_Mar10.pdf), and [DY 2011](http://www.ssc.upenn.edu/~fdiebold/papers/paper106/dieboldandyilmaz2011.pdf), as well as scaling the methods for quick application and replication. The DY papers describe how to quantify in an index how much error in a forecast can be attributed to shocks to its variables. These "spillovers" are found using standard time-series techniques for finding an impulse-response function and forecast error variance decompositions. In addition, DY 2011 show how the effects of these shocks can be represented as a network when one relaxes the assumptions that vector autoregressions have normally distributed, i.i.d. error terms and that shocks to each variable must be orthongonal to the other shocks. In the spirit of DY 2011, stocks of large banks and AIG will be the subjects of analysis. Like all the DY publications, the analysis will focus on connections and spillovers of volatility.

This vignette will proceed as follows: Section 2) describe data source and summarize data; Section 3) to motivate the concepts of spillovers and connectedness a full sample analysis of the data will be done; Section 4) using a rolling window over the time-series data, apply the same concepts from the full sample analysis to derive index levels of spillover/conncectedness and network behavior over time.

Section 2: Data
---------------

The data is pulled from Yahoo! Finance using the getSymbols() function from the quantmod package, which includes Open, Close, High, and Low prices for a user specified set of stocks/indexes over a chosen period of time measured in days. The data is similar to that in [DY 2010](http://www.ssc.upenn.edu/~fdiebold/papers/paper99/DirectionofSpillovers_Mar10.pdf), so measures of stock return volatility used in the analysis and app are the same as in that publication. First daily variance for stock *i* at time *t* is estimated using high and low prices:

To be able to analyze the volatility of closed-end mutual funds the minimum price is observed at *t* − 1. Otherwise, returns and volatilities for these funds would be zero.

Since voltatilities are skewed, it is common practice to use log-volatilites which closely approximate a normally distribution. However, to control for instances when volatility is zero, *s**i**n**h*<sup>−1</sup> is used instead of taking the log (*s**i**n**h*<sup>−1</sup>(*x*)=*l**o**g*(2*x*)).

Using the quantmod getSymbols function, read in the set of stock prices for a period of time. Then calculate their daily volatilities.

    ##  [1] "AIG"  "BAC"  "BK"   "BNS"  "C"    "CS"   "DB"   "GS"   "HSBC" "JPM" 
    ## [11] "MS"   "TD"   "WFC"

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
count
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
med
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
min
</th>
<th style="text-align:right;">
max
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BK
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2414055
</td>
<td style="text-align:right;">
0.1739952
</td>
<td style="text-align:right;">
0.2355072
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.348819
</td>
</tr>
<tr>
<td style="text-align:left;">
JPM
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2412604
</td>
<td style="text-align:right;">
0.1746596
</td>
<td style="text-align:right;">
0.2310831
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.992054
</td>
</tr>
<tr>
<td style="text-align:left;">
GS
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2413190
</td>
<td style="text-align:right;">
0.1802763
</td>
<td style="text-align:right;">
0.2276783
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.215029
</td>
</tr>
<tr>
<td style="text-align:left;">
C
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2922059
</td>
<td style="text-align:right;">
0.1889913
</td>
<td style="text-align:right;">
0.3121762
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.643312
</td>
</tr>
<tr>
<td style="text-align:left;">
MS
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.3061973
</td>
<td style="text-align:right;">
0.2229228
</td>
<td style="text-align:right;">
0.2916601
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.954439
</td>
</tr>
<tr>
<td style="text-align:left;">
BNS
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.1815797
</td>
<td style="text-align:right;">
0.1339244
</td>
<td style="text-align:right;">
0.1753930
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.553637
</td>
</tr>
<tr>
<td style="text-align:left;">
AIG
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.3065092
</td>
<td style="text-align:right;">
0.1883449
</td>
<td style="text-align:right;">
0.3500704
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3.567131
</td>
</tr>
<tr>
<td style="text-align:left;">
TD
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.1742599
</td>
<td style="text-align:right;">
0.1295371
</td>
<td style="text-align:right;">
0.1651353
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.646948
</td>
</tr>
<tr>
<td style="text-align:left;">
HSBC
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.1610396
</td>
<td style="text-align:right;">
0.1141410
</td>
<td style="text-align:right;">
0.1648465
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.679011
</td>
</tr>
<tr>
<td style="text-align:left;">
DB
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2509589
</td>
<td style="text-align:right;">
0.1842325
</td>
<td style="text-align:right;">
0.2367057
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.873443
</td>
</tr>
<tr>
<td style="text-align:left;">
BAC
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2811012
</td>
<td style="text-align:right;">
0.1944293
</td>
<td style="text-align:right;">
0.2905163
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.469604
</td>
</tr>
<tr>
<td style="text-align:left;">
CS
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2287849
</td>
<td style="text-align:right;">
0.1654587
</td>
<td style="text-align:right;">
0.2200668
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.848902
</td>
</tr>
<tr>
<td style="text-align:left;">
WFC
</td>
<td style="text-align:right;">
2767
</td>
<td style="text-align:right;">
0.2406566
</td>
<td style="text-align:right;">
0.1612006
</td>
<td style="text-align:right;">
0.2528949
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.152670
</td>
</tr>
</tbody>
</table>
``` r
chart.TimeSeries(vol.data,lwd=2,auto.grid=F,ylab="Annualized Log Volatility",xlab="Time",
                     main="Log Volatility",lty=1,
                     legend.loc="topright")
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

Section 3: Full Sample Analysis
-------------------------------

Section 3.1: Basic Time-Series: VAR(p), MA(q), Impulse-Response, FEVD
---------------------------------------------------------------------

In time-series analysis (and most regressions generally) the most interesting information is found in the distribution of the error terms, especially the variance. In a multivariate model the distibution of error terms embeds the assumptions made regarding correlations of dependent variables. In financial econometrics forecasting stock price returns entails two components - systematic and idiosyncratic shown here:

Where *R*<sub>*i**t*</sub> is the return on stock *i* at time *t* and *Σ*<sub>*M*</sub> is the market return (i.e. the return on the S&P 500, Dow Jones, etc.), *β* is a measure of the systematic portion of a stock's return, or the strength/measure of its relationship with the broader market, and *ϵ* is the idiosyncratic portion. The interesting question is how responsive is the future volatility *σ*<sub>*i*, *t* + 1</sub>, at time *t* + 1, to a unit shock of *ϵ*<sub>*i**t*</sub>? And in a system of several different returns, how does a shock to the idiosyncratic term *ϵ*<sub>*j**t*</sub> for variable *j* at time *t* effect *σ*<sub>*i*, *t* + 1</sub> for variable *i* at time *t* + 1? In other words, we are interested in finding whether and how idiosyncratic risks of one stock effcts others. These effects are known as spillovers or connectedness.

In order to derive a measure of these spillovers, the calculations will proceed in four steps: 1) Estimate a VAR for a set of *N* variables of *σ*<sub>*i*</sub> for *i* = 1, ..,*N* variables which will replace *σ*<sub>*M*</sub>; 2) convert the VAR into a moving average representation of order q (MA(q)); 3) derive the impulse-response function for each variable using the MA(q) for each variable *j* on variable *i*; 4) calculate FEVD matrix for the entire set of variables using the impulse-response function. For textbook treatment of how this is done algebraically, see [Zivot](https://faculty.washington.edu/ezivot/econ584/notes/varModels.pdf), [Cochrane Chapter 7](http://econ.lse.ac.uk/staff/wdenhaan/teach/cochrane.pdf), and [Pesaran and Shin 1997](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.153.5596&rep=rep1&type=pdf) Section 2 for estimating both correlated an uncorrelated shocks.

Section 3.2: Full Sample Spillover Table
----------------------------------------

### Step 1. Vector Autoregression (VAR)

Using the VAR and SVAR functions from the [vars](https://cran.r-project.org/web/packages/vars/vars.pdf) package using the set of volatilities estimate a vector autoregression for the volatilities. The SVAR function estimates a Structural VAR which imposes a lower triangular matrix of coefficients (with 1s along the diagonal). This is necessary in order not to violate that assumption that the error terms are normally distributed and i.i.d., which will matter when impulse response functions are calculated for orthogonal and generalized impulses. Thus, the VAR is "non-orthogonal" and the SVAR is "orthogonal". From these models we will need the residual errors (i.e. the idiosyncratic portion) for the calculations that follow.

``` r
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

``` r
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

Second, mulitply each lag of MA cofficients by the Cholesky matrix. The FEVD is estimated by picking out the MA coefficients for variable *j* and multiplying those by the shocks from variable *i* as shown in the formula below. The results are shown in the following table.

*θ*<sub>*i**j*</sub><sup>*o*</sup> is the orthogonal FEVD, *A*<sub>*l*</sub> is the MA coefficient matrix at lag *l* out of *n* lags, *Σ* is the residual covariance matrix, *P* is the lower triangular matrix, and *e*<sub>*i*</sub> and *e*<sub>*j*</sub> are basis vectors with unity at index *i* and *j*, respectively.

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
BK
</th>
<th style="text-align:right;">
JPM
</th>
<th style="text-align:right;">
GS
</th>
<th style="text-align:right;">
C
</th>
<th style="text-align:right;">
MS
</th>
<th style="text-align:right;">
BNS
</th>
<th style="text-align:right;">
AIG
</th>
<th style="text-align:right;">
TD
</th>
<th style="text-align:right;">
HSBC
</th>
<th style="text-align:right;">
DB
</th>
<th style="text-align:right;">
BAC
</th>
<th style="text-align:right;">
CS
</th>
<th style="text-align:right;">
WFC
</th>
<th style="text-align:right;">
From
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BK
</td>
<td style="text-align:right;">
69
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
31
</td>
</tr>
<tr>
<td style="text-align:left;">
JPM
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
32
</td>
</tr>
<tr>
<td style="text-align:left;">
GS
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
34
</td>
</tr>
<tr>
<td style="text-align:left;">
C
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
45
</td>
</tr>
<tr>
<td style="text-align:left;">
MS
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
59
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
41
</td>
</tr>
<tr>
<td style="text-align:left;">
BNS
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
79
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
AIG
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
45
</td>
</tr>
<tr>
<td style="text-align:left;">
TD
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
25
</td>
</tr>
<tr>
<td style="text-align:left;">
HSBC
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
26
</td>
</tr>
<tr>
<td style="text-align:left;">
DB
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
38
</td>
</tr>
<tr>
<td style="text-align:left;">
BAC
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
53
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
47
</td>
</tr>
<tr>
<td style="text-align:left;">
CS
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
65
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
35
</td>
</tr>
<tr>
<td style="text-align:left;">
WFC
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
33
</td>
</tr>
<tr>
<td style="text-align:left;">
To
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
86
</td>
<td style="text-align:right;">
35
</td>
</tr>
</tbody>
</table>
This is the Spillover Table. The row sums (excluding the diagonal terms) measure how much volatility the stock adds to the system. The column sums measure how much volatility each stock receives from the system. Hence, each cell in the table is a measure of how much volatility from column *i* is given to row *j*. One can estimate the net pairwise spillover by subtracting two spillovers. For instance, the spillover from Goldman Sachs (GS) to Morgan Stanley (MS) is 2, and from MS to GS is 5. Therefore, the net pairwise spillover from GS to MS is -3.

Section 3.3: Full Sample Connectedness Table
--------------------------------------------

Connectedness requires similar calculations, except that the shocks from IRF are no longer assumed orthogonal. To achieve this, the lower triangular matrix *P* for *θ*<sub>*i**j*</sub><sup>*o*</sup> in the numerator is replaced by the covariance matrix of the residuals. The generalized FEVD equation is then,

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
BK
</th>
<th style="text-align:right;">
JPM
</th>
<th style="text-align:right;">
GS
</th>
<th style="text-align:right;">
C
</th>
<th style="text-align:right;">
MS
</th>
<th style="text-align:right;">
BNS
</th>
<th style="text-align:right;">
AIG
</th>
<th style="text-align:right;">
TD
</th>
<th style="text-align:right;">
HSBC
</th>
<th style="text-align:right;">
DB
</th>
<th style="text-align:right;">
BAC
</th>
<th style="text-align:right;">
CS
</th>
<th style="text-align:right;">
WFC
</th>
<th style="text-align:right;">
From
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BK
</td>
<td style="text-align:right;">
10.62
</td>
<td style="text-align:right;">
7.42
</td>
<td style="text-align:right;">
6.61
</td>
<td style="text-align:right;">
12.10
</td>
<td style="text-align:right;">
12.58
</td>
<td style="text-align:right;">
3.14
</td>
<td style="text-align:right;">
10.46
</td>
<td style="text-align:right;">
2.93
</td>
<td style="text-align:right;">
2.85
</td>
<td style="text-align:right;">
6.96
</td>
<td style="text-align:right;">
10.28
</td>
<td style="text-align:right;">
5.92
</td>
<td style="text-align:right;">
8.13
</td>
<td style="text-align:right;">
89.38
</td>
</tr>
<tr>
<td style="text-align:left;">
JPM
</td>
<td style="text-align:right;">
7.29
</td>
<td style="text-align:right;">
9.50
</td>
<td style="text-align:right;">
6.72
</td>
<td style="text-align:right;">
12.97
</td>
<td style="text-align:right;">
12.18
</td>
<td style="text-align:right;">
3.17
</td>
<td style="text-align:right;">
9.64
</td>
<td style="text-align:right;">
2.94
</td>
<td style="text-align:right;">
3.00
</td>
<td style="text-align:right;">
6.89
</td>
<td style="text-align:right;">
11.11
</td>
<td style="text-align:right;">
5.75
</td>
<td style="text-align:right;">
8.84
</td>
<td style="text-align:right;">
90.50
</td>
</tr>
<tr>
<td style="text-align:left;">
GS
</td>
<td style="text-align:right;">
7.17
</td>
<td style="text-align:right;">
7.38
</td>
<td style="text-align:right;">
10.58
</td>
<td style="text-align:right;">
11.32
</td>
<td style="text-align:right;">
14.41
</td>
<td style="text-align:right;">
3.16
</td>
<td style="text-align:right;">
9.18
</td>
<td style="text-align:right;">
2.90
</td>
<td style="text-align:right;">
3.06
</td>
<td style="text-align:right;">
7.73
</td>
<td style="text-align:right;">
9.43
</td>
<td style="text-align:right;">
6.57
</td>
<td style="text-align:right;">
7.13
</td>
<td style="text-align:right;">
89.42
</td>
</tr>
<tr>
<td style="text-align:left;">
C
</td>
<td style="text-align:right;">
6.44
</td>
<td style="text-align:right;">
7.09
</td>
<td style="text-align:right;">
5.87
</td>
<td style="text-align:right;">
17.90
</td>
<td style="text-align:right;">
11.19
</td>
<td style="text-align:right;">
2.92
</td>
<td style="text-align:right;">
10.67
</td>
<td style="text-align:right;">
2.78
</td>
<td style="text-align:right;">
2.90
</td>
<td style="text-align:right;">
6.47
</td>
<td style="text-align:right;">
11.96
</td>
<td style="text-align:right;">
5.52
</td>
<td style="text-align:right;">
8.29
</td>
<td style="text-align:right;">
82.10
</td>
</tr>
<tr>
<td style="text-align:left;">
MS
</td>
<td style="text-align:right;">
7.27
</td>
<td style="text-align:right;">
7.18
</td>
<td style="text-align:right;">
7.63
</td>
<td style="text-align:right;">
11.69
</td>
<td style="text-align:right;">
16.81
</td>
<td style="text-align:right;">
2.99
</td>
<td style="text-align:right;">
9.84
</td>
<td style="text-align:right;">
2.73
</td>
<td style="text-align:right;">
2.93
</td>
<td style="text-align:right;">
7.39
</td>
<td style="text-align:right;">
10.06
</td>
<td style="text-align:right;">
6.31
</td>
<td style="text-align:right;">
7.18
</td>
<td style="text-align:right;">
83.19
</td>
</tr>
<tr>
<td style="text-align:left;">
BNS
</td>
<td style="text-align:right;">
6.91
</td>
<td style="text-align:right;">
6.93
</td>
<td style="text-align:right;">
6.25
</td>
<td style="text-align:right;">
11.63
</td>
<td style="text-align:right;">
11.37
</td>
<td style="text-align:right;">
7.02
</td>
<td style="text-align:right;">
9.28
</td>
<td style="text-align:right;">
5.05
</td>
<td style="text-align:right;">
3.36
</td>
<td style="text-align:right;">
7.88
</td>
<td style="text-align:right;">
9.63
</td>
<td style="text-align:right;">
7.00
</td>
<td style="text-align:right;">
7.70
</td>
<td style="text-align:right;">
92.98
</td>
</tr>
<tr>
<td style="text-align:left;">
AIG
</td>
<td style="text-align:right;">
5.85
</td>
<td style="text-align:right;">
5.95
</td>
<td style="text-align:right;">
4.62
</td>
<td style="text-align:right;">
12.24
</td>
<td style="text-align:right;">
9.81
</td>
<td style="text-align:right;">
2.81
</td>
<td style="text-align:right;">
25.34
</td>
<td style="text-align:right;">
2.56
</td>
<td style="text-align:right;">
2.83
</td>
<td style="text-align:right;">
5.50
</td>
<td style="text-align:right;">
10.08
</td>
<td style="text-align:right;">
4.67
</td>
<td style="text-align:right;">
7.73
</td>
<td style="text-align:right;">
74.66
</td>
</tr>
<tr>
<td style="text-align:left;">
TD
</td>
<td style="text-align:right;">
6.90
</td>
<td style="text-align:right;">
7.08
</td>
<td style="text-align:right;">
6.27
</td>
<td style="text-align:right;">
12.17
</td>
<td style="text-align:right;">
11.44
</td>
<td style="text-align:right;">
5.32
</td>
<td style="text-align:right;">
9.33
</td>
<td style="text-align:right;">
5.95
</td>
<td style="text-align:right;">
3.23
</td>
<td style="text-align:right;">
7.97
</td>
<td style="text-align:right;">
9.74
</td>
<td style="text-align:right;">
6.73
</td>
<td style="text-align:right;">
7.87
</td>
<td style="text-align:right;">
94.05
</td>
</tr>
<tr>
<td style="text-align:left;">
HSBC
</td>
<td style="text-align:right;">
6.39
</td>
<td style="text-align:right;">
6.55
</td>
<td style="text-align:right;">
6.07
</td>
<td style="text-align:right;">
12.21
</td>
<td style="text-align:right;">
11.00
</td>
<td style="text-align:right;">
3.59
</td>
<td style="text-align:right;">
9.23
</td>
<td style="text-align:right;">
3.08
</td>
<td style="text-align:right;">
7.24
</td>
<td style="text-align:right;">
9.47
</td>
<td style="text-align:right;">
10.02
</td>
<td style="text-align:right;">
7.93
</td>
<td style="text-align:right;">
7.22
</td>
<td style="text-align:right;">
92.76
</td>
</tr>
<tr>
<td style="text-align:left;">
DB
</td>
<td style="text-align:right;">
6.40
</td>
<td style="text-align:right;">
6.35
</td>
<td style="text-align:right;">
6.15
</td>
<td style="text-align:right;">
10.48
</td>
<td style="text-align:right;">
11.16
</td>
<td style="text-align:right;">
3.42
</td>
<td style="text-align:right;">
8.13
</td>
<td style="text-align:right;">
3.12
</td>
<td style="text-align:right;">
4.04
</td>
<td style="text-align:right;">
15.52
</td>
<td style="text-align:right;">
9.64
</td>
<td style="text-align:right;">
9.26
</td>
<td style="text-align:right;">
6.35
</td>
<td style="text-align:right;">
84.48
</td>
</tr>
<tr>
<td style="text-align:left;">
BAC
</td>
<td style="text-align:right;">
6.66
</td>
<td style="text-align:right;">
7.33
</td>
<td style="text-align:right;">
5.61
</td>
<td style="text-align:right;">
14.31
</td>
<td style="text-align:right;">
10.85
</td>
<td style="text-align:right;">
2.95
</td>
<td style="text-align:right;">
9.99
</td>
<td style="text-align:right;">
2.75
</td>
<td style="text-align:right;">
2.99
</td>
<td style="text-align:right;">
6.87
</td>
<td style="text-align:right;">
15.41
</td>
<td style="text-align:right;">
5.52
</td>
<td style="text-align:right;">
8.76
</td>
<td style="text-align:right;">
84.59
</td>
</tr>
<tr>
<td style="text-align:left;">
CS
</td>
<td style="text-align:right;">
6.47
</td>
<td style="text-align:right;">
6.37
</td>
<td style="text-align:right;">
6.35
</td>
<td style="text-align:right;">
10.93
</td>
<td style="text-align:right;">
11.52
</td>
<td style="text-align:right;">
3.42
</td>
<td style="text-align:right;">
8.54
</td>
<td style="text-align:right;">
3.07
</td>
<td style="text-align:right;">
3.79
</td>
<td style="text-align:right;">
10.82
</td>
<td style="text-align:right;">
9.48
</td>
<td style="text-align:right;">
12.70
</td>
<td style="text-align:right;">
6.55
</td>
<td style="text-align:right;">
87.30
</td>
</tr>
<tr>
<td style="text-align:left;">
WFC
</td>
<td style="text-align:right;">
7.29
</td>
<td style="text-align:right;">
7.61
</td>
<td style="text-align:right;">
5.73
</td>
<td style="text-align:right;">
13.81
</td>
<td style="text-align:right;">
10.91
</td>
<td style="text-align:right;">
3.07
</td>
<td style="text-align:right;">
10.07
</td>
<td style="text-align:right;">
2.93
</td>
<td style="text-align:right;">
3.01
</td>
<td style="text-align:right;">
6.19
</td>
<td style="text-align:right;">
12.07
</td>
<td style="text-align:right;">
5.22
</td>
<td style="text-align:right;">
12.08
</td>
<td style="text-align:right;">
87.92
</td>
</tr>
<tr>
<td style="text-align:left;">
To
</td>
<td style="text-align:right;">
81.04
</td>
<td style="text-align:right;">
83.24
</td>
<td style="text-align:right;">
73.89
</td>
<td style="text-align:right;">
145.84
</td>
<td style="text-align:right;">
138.41
</td>
<td style="text-align:right;">
39.96
</td>
<td style="text-align:right;">
114.36
</td>
<td style="text-align:right;">
36.84
</td>
<td style="text-align:right;">
37.98
</td>
<td style="text-align:right;">
90.14
</td>
<td style="text-align:right;">
123.50
</td>
<td style="text-align:right;">
76.40
</td>
<td style="text-align:right;">
91.74
</td>
<td style="text-align:right;">
87.18
</td>
</tr>
</tbody>
</table>
The interpretation of the table as well as the calculation for net pairwise connectedness is not the same as that of the spillover. In DY 2011 the connectedness table can be interpreted as a network in which each entry represents a weight of directional connections between the nodes. The row sums are the from-degrees, column sums are to-degrees. The total connectedness (bottom-right entry in the table) is the average degree of the network (to or from since the row sums equal column sums). Since generalized FEVD matrices do not have row/columns sums equal to one, the entries in the connectedness have been normalized to that end.

Section 3.4: Full Sample FEVD as a Network
------------------------------------------

DY 2011 demonstrates that variance decompositions define weighted, directed networks. Thus the Connectedness table can be transformed into a network as shown below. (To see how the network evolves over time, use the app here: (insert link))

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

Section 4: Rolling Window - Spillover Vs. Connected
---------------------------------------------------

To apply the methods shown above to a rolling window, it is neccessary to only define the size of the rolling window, the size of the increment to move the window, which will determine the how many times to loop through the VAR function. The function vector\_autoreg takes in the an .xts data frame, window size, increment, AR lag order, and MA lag order. It also calls on function "make.table" and "normalize" which can be found in the function script on the GitHub page.

``` r
vector_autoreg <- function(theData.xts,window=100,iter=10,ar_lag=3,ma_lag=10){
  ### create a sequence to use as cut offs for the rolling window
  ### estimation. iter determines how far to move the window forward
  
  mySeq <- seq(1,nrow(theData.xts)-window-1,by=iter)
  min_iter = min(diff(index(theData.xts[mySeq])),na.rm=T)
  date_names <- as.character(seq.Date(from=min(index(theData.xts)),by=min_iter,length.out = length(mySeq)))
  
  model_list <- vector("list",length=length(date_names))
  names(model_list) <- date_names
  centrality <- model_list
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
    colnames(net.melt)[3] <- "weight"

    ### calculate each percentile of the net pairwise connectedness values
    ### and choose only the top 10%
    net.quant <- quantile(net.melt$weight,prob=seq(0,1,by=0.01))
    new.net1 <- net.melt[net.melt$weight > net.quant[90],]
    new.net <- new.net1[new.net1[,1] != new.net1[,2],]
    
    ### create igraph graph
    total_network <- graph.data.frame(net.melt,direct=T)
    E(total_network)$weight <- as.numeric(net.melt$weight)
    
    net_clo <- centr_clo(total_network)$centralization
    net_betw <- centr_betw(total_network,directed = TRUE,normalized = F)$centralization
    net_deg <- centr_degree(total_network,mode="total")$centralization
    
    centrality[[date_names[d]]] <- list(clo = net_clo,betw = net_betw,deg = net_deg)
    
    net.network <- graph.data.frame(new.net,direct=T)
    E(net.network)$weight <- as.numeric(new.net$weight)
    ### set graph nodes
    #V(net.network)
    ### set edge colors
    E(net.network)$color <- ifelse(E(net.network)$weight >= net.quant[99],"black",
                                   ifelse(E(net.network)$weight < net.quant[99] & E(net.network)$weight >= net.quant[95],"red",
                                          ifelse(E(net.network)$weight < net.quant[95] & E(net.network)$weight >= net.quant[90],"orange","blue")))
    
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
       centrality = centrality,
       dates = date_names)
  
}

rolling_var <- vector_autoreg(theData.xts = vol.data)
```

Section 4.1: Spillover & Connectedness Indexes
==============================================

``` r
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

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

The indexes show the percentage of the 10-day ahead forecast error which can be explained by spillovers/connectedness. Mathematically, for each matrix calculated at each point in time the sum of all the off-diagonal entries are divided by the sum of the entire matrix. Surpisingly, the indexes are very correlated, hence the use of orthogonal shocks yields a similar amount of forecast error as correlated shocks in this instance.

Section 4.2: Total Directional Connectedness and Spillover
----------------------------------------------------------

Using the time-series of connectedness and spillover tables it is also possible to see how much volatility each stock gave or received, and which ones were net givers/receivers of volatility over time. The red lines are the Spillover measures and the black lines are Connectedness.

These results show the consequences of the assumptions behind the SVAR and VAR models. During the run up to the 2008 financial crisis and in its wake, AIG underwrote several option contracts for various asset-backed securities with the large banks. However, the SVAR does not pick up on that connection, unlike the "non-orthogonal" VAR, and instead estimates a lower volatility spillover of AIG stock price to the other prices. A likely reason is that SVARs are sensitive to the ordering of variables used to satisfy the normality assumption, and, perhaps, a different ordering would deliver a different result. In addition, these results will change for the SVAR when the parameters to the function vector\_autoreg() are altered.

``` r
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

``` r
ggplot(to.tbl) + 
  geom_line(aes(x=Date,y=Connectedness),colour="black") +
  geom_line(aes(x=Date,y=Spillover),colour="red") +
  facet_wrap(~variable,scales="free_y") +
  theme(legend.position = "bottom")
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

### Receiving

``` r
ggplot(from.tbl) + 
  geom_line(aes(x=Date,y=Connectedness),colour="black") +
  geom_line(aes(x=Date,y=Spillover),colour="red") +
  facet_wrap(~variable,scales="free_y") +
  theme(legend.position = "bottom")
```

![](README_files/figure-markdown_github/unnamed-chunk-12-1.png)

### Net

``` r
ggplot(net.tbl) + 
  geom_line(aes(x=Date,y=Connectedness),colour="black") +
  geom_line(aes(x=Date,y=Spillover),colour="red") +
  facet_wrap(~variable,scales="free_y") +
  theme(legend.position = "bottom")
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

Section 4.3: Rolling Window Net Pairwise Spillovers & Connectedness
-------------------------------------------------------------------

Using the connectedness measure, it is also possible to compare two stocks and see which one was the net giver or receiver. For instance if we compare the Morgan Stanley (MS) and Goldman Sachs (GS), it appears that Morgan Stanley "gave" volatility to Goldman Sachs for most of the time considered.

``` r
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

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)

Section 5: Monetary Transmission
--------------------------------

Research from [Niepmann & Schmidt-Eisenlohr](https://voxeu.org/article/when-dollar-appreciates-us-corporate-credit-tightens) shows that when the US dollar appreciates business loan terms tighten. The main reason for this is that banks often package and sell their business loans in capital markets to Collateralized Loan Obligations (CLO) insurance companies, mutual funds, hedge funds, and sometimes other banks, which can have exposure to foreign exchange markets. Hence, when the dollar appreciates these entities will buy fewer business loans, and this lack of demand causes banks to tighten loan terms so as to offer better returns in the capital markets. This development is significant for several reasons, but perhaps the most salient point is that accroding to economic theory currency appreciation is associated with asset higher returns (i.e. interest rates), as this graph from [Krugman, Obstfeld, and Melitz](https://www.amazon.com/International-Economics-Theory-Policy-Pearson/dp/0133423646) demonstrates:

![International Economics: Theory and Policy](C:\Users\Rhamey\Desktop\financeR\ConnectedNess\krugman_obstfeld.png)

Using the methodology above and few proxies for relevant financial assets we can get an idea of how connected loan demand and foreign exchange are. According to theory when interest rates the dollar will appreciate and lower returns, hence a proxy investment for dollar exchange rate, relative to EMs and developed markets, and money market returns is necessary. To measure the transmission of emerging markets to CLOs it is necessary need measure the returns for buyers and sellers of CLOs. Tthe price of oil is an [important factor](https://www.bankofcanada.ca/wp-content/uploads/2010/05/wp10-5.pdf) driving exchange rates, thus the price of WTI will be necessary. Finally, Citibank, JP Morgan, and Wells Fargo are among the largest sellers (buyer in Citibank's case) as discussed in a [note from S&P](https://www.spglobal.com/marketintelligence/en/news-insights/latest-news-headlines/leveraged-loan-news/those-700b-in-us-clos-who-holds-them-what-risk-they-pose).

Luckily there is an ETF for everything. Consider the following:

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Ticker
</th>
<th style="text-align:left;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BSCK
</td>
<td style="text-align:left;">
Money Market Proxy (ETF)
</td>
</tr>
<tr>
<td style="text-align:left;">
VCSH
</td>
<td style="text-align:left;">
Money Market Proxy (ETF)
</td>
</tr>
<tr>
<td style="text-align:left;">
USO
</td>
<td style="text-align:left;">
Light, sweet crude (ETF)
</td>
</tr>
<tr>
<td style="text-align:left;">
LEMB
</td>
<td style="text-align:left;">
Emerging Market Currencies (excl. China) (ETF)
</td>
</tr>
<tr>
<td style="text-align:left;">
FXCH
</td>
<td style="text-align:left;">
Chinese/USD Yuan (ETF)
</td>
</tr>
<tr>
<td style="text-align:left;">
UUP
</td>
<td style="text-align:left;">
USD Index relative to developed economies (ETF)
</td>
</tr>
<tr>
<td style="text-align:left;">
FFRHX
</td>
<td style="text-align:left;">
Fidelity Floating Rate High Income Fund (Closed-End Fund)
</td>
</tr>
<tr>
<td style="text-align:left;">
OXLC
</td>
<td style="text-align:left;">
Oxford Lane Corporation (Closed End Fund)
</td>
</tr>
<tr>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
Citibank
</td>
</tr>
<tr>
<td style="text-align:left;">
WFC
</td>
<td style="text-align:left;">
Wells Fargo
</td>
</tr>
<tr>
<td style="text-align:left;">
JPM
</td>
<td style="text-align:left;">
JP Morgan Chase & Co.
</td>
</tr>
</tbody>
</table>
    ##  [1] "BSCK"  "C"     "FFRHX" "FXCH"  "JPM"   "LEMB"  "OXLC"  "USO"  
    ##  [9] "UUP"   "VCSH"  "WFC"

Price History (Log Prices)
--------------------------

First, the black line in the graph that falls precipitously is USO (oil), the other black line is UUP (USD exhange rate with other developed economies), and the green line at the bottom is Oxford Lane Capital, a purchaser of CLO securities. Of course the dollar and oil are negatively correlated, but it should be surprising that the prices of oil and OXLC move together. Moreover, it appears that as the dollar appreciates against the yen (FXHC, top blue line) and other emerging markets (LEMB, pink line) OXLC also falls in price.

``` r
chart.TimeSeries(all_data[["log_price"]],lwd=2,auto.grid=F,ylab="Log Prices",xlab="Time",
                     main="Log Prices",lty=1,
                     legend.loc="topright")
```

![](README_files/figure-markdown_github/unnamed-chunk-17-1.png)

Cumulative Returns
------------------

``` r
chart.TimeSeries(all_data[["cum_returns"]],lwd=2,auto.grid=F,ylab="Cumulative Returns",xlab="Time",
                     main="Returns (%)",lty=1,
                     legend.loc="topright")
```

![](README_files/figure-markdown_github/unnamed-chunk-18-1.png)

Daily Volatility
----------------

``` r
chart.TimeSeries(all_data[["return_vol"]],lwd=2,auto.grid=F,ylab="Annualized Log Volatility",xlab="Time",
                     main="Log Volatility",lty=1,
                     legend.loc="topright")
```

![](README_files/figure-markdown_github/unnamed-chunk-19-1.png)

Daily Returns
-------------

``` r
chart.TimeSeries(all_data[["returns"]],lwd=2,auto.grid=F,ylab="Return (%)",xlab="Time",
                     main="Daily Return Percentage",lty=1,
                     legend.loc="topright")
```

![](README_files/figure-markdown_github/unnamed-chunk-20-1.png)

### Network

``` r
vol_var = VAR(vol.data,p=3,type="none")
amat <- diag(ncol(vol.data))
amat[lower.tri(amat)] <- NA
### extract residuals of the VAR
res_t <- residuals(vol_var)
MA_lag <- 10
theta_temp <- Phi(vol_var,nstep = MA_lag)
### extract MA coefficients
theta.list <- alply(theta_temp,3)
dir.mat <- directional.matrix(res_t,theta.list)
rownames(dir.mat) <- colnames(vol.data)
colnames(dir.mat) <- rownames(dir.mat)

D <- normalize(dir.mat)*100
df.list <- make.table(D,rownames(D),start_date)
df.list <- round(df.list[["table"]],2)

kable(df.list) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
UUP
</th>
<th style="text-align:right;">
JPM
</th>
<th style="text-align:right;">
OXLC
</th>
<th style="text-align:right;">
FXCH
</th>
<th style="text-align:right;">
BSCK
</th>
<th style="text-align:right;">
LEMB
</th>
<th style="text-align:right;">
C
</th>
<th style="text-align:right;">
VCSH
</th>
<th style="text-align:right;">
USO
</th>
<th style="text-align:right;">
FFRHX
</th>
<th style="text-align:right;">
WFC
</th>
<th style="text-align:right;">
From
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
UUP
</td>
<td style="text-align:right;">
13.19
</td>
<td style="text-align:right;">
9.29
</td>
<td style="text-align:right;">
15.04
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
3.09
</td>
<td style="text-align:right;">
13.53
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
34.35
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
10.06
</td>
<td style="text-align:right;">
86.81
</td>
</tr>
<tr>
<td style="text-align:left;">
JPM
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
31.49
</td>
<td style="text-align:right;">
6.26
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
30.96
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
12.14
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
17.76
</td>
<td style="text-align:right;">
68.51
</td>
</tr>
<tr>
<td style="text-align:left;">
OXLC
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
4.66
</td>
<td style="text-align:right;">
68.21
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
8.54
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
13.31
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
4.17
</td>
<td style="text-align:right;">
31.79
</td>
</tr>
<tr>
<td style="text-align:left;">
FXCH
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
2.98
</td>
<td style="text-align:right;">
3.07
</td>
<td style="text-align:right;">
74.25
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.41
</td>
<td style="text-align:right;">
3.27
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
11.62
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
2.76
</td>
<td style="text-align:right;">
25.75
</td>
</tr>
<tr>
<td style="text-align:left;">
BSCK
</td>
<td style="text-align:right;">
1.39
</td>
<td style="text-align:right;">
13.35
</td>
<td style="text-align:right;">
16.21
</td>
<td style="text-align:right;">
1.39
</td>
<td style="text-align:right;">
10.50
</td>
<td style="text-align:right;">
6.78
</td>
<td style="text-align:right;">
17.52
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
23.84
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
8.64
</td>
<td style="text-align:right;">
89.50
</td>
</tr>
<tr>
<td style="text-align:left;">
LEMB
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
12.92
</td>
<td style="text-align:right;">
10.19
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
19.99
</td>
<td style="text-align:right;">
17.15
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
28.15
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
8.35
</td>
<td style="text-align:right;">
80.01
</td>
</tr>
<tr>
<td style="text-align:left;">
C
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
22.52
</td>
<td style="text-align:right;">
7.37
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
39.33
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
13.27
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
16.10
</td>
<td style="text-align:right;">
60.67
</td>
</tr>
<tr>
<td style="text-align:left;">
VCSH
</td>
<td style="text-align:right;">
1.40
</td>
<td style="text-align:right;">
11.65
</td>
<td style="text-align:right;">
14.80
</td>
<td style="text-align:right;">
1.21
</td>
<td style="text-align:right;">
1.79
</td>
<td style="text-align:right;">
5.32
</td>
<td style="text-align:right;">
17.54
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
34.22
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
11.03
</td>
<td style="text-align:right;">
98.96
</td>
</tr>
<tr>
<td style="text-align:left;">
USO
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
6.36
</td>
<td style="text-align:right;">
6.18
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
1.21
</td>
<td style="text-align:right;">
9.97
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
69.66
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
5.93
</td>
<td style="text-align:right;">
30.34
</td>
</tr>
<tr>
<td style="text-align:left;">
FFRHX
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
13.90
</td>
<td style="text-align:right;">
18.04
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
1.05
</td>
<td style="text-align:right;">
22.11
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
31.89
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
11.34
</td>
<td style="text-align:right;">
99.42
</td>
</tr>
<tr>
<td style="text-align:left;">
WFC
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
21.54
</td>
<td style="text-align:right;">
6.48
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
28.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
13.26
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
29.45
</td>
<td style="text-align:right;">
70.55
</td>
</tr>
<tr>
<td style="text-align:left;">
To
</td>
<td style="text-align:right;">
6.03
</td>
<td style="text-align:right;">
119.16
</td>
<td style="text-align:right;">
103.63
</td>
<td style="text-align:right;">
6.08
</td>
<td style="text-align:right;">
3.82
</td>
<td style="text-align:right;">
22.04
</td>
<td style="text-align:right;">
168.60
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
216.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
96.13
</td>
<td style="text-align:right;">
67.48
</td>
</tr>
</tbody>
</table>
``` r
net.mat = t(D) - D
net.melt1 <- melt(net.mat)
net.melt <- net.melt1 %>% filter(value != 0 )
colnames(net.melt)[3] <- "weight"

### calculate each percentile of the net pairwise connectedness values
### and choose only the top 10%
net.quant <- quantile(net.melt$weight,prob=seq(0,1,by=0.01))
new.net1 <- net.melt[net.melt$weight >= net.quant[75],]
new.net <- new.net1[new.net1[,1] != new.net1[,2],]

### create igraph graph
net.network <- graph.data.frame(new.net,direct=T)
#net.network <- set_edge_attr(net.network, "weight", value= new.net$weight)
E(net.network)$weight <- as.numeric(new.net$weight)
### set graph nodes
#V(net.network)
### set edge colors
E(net.network)$color <- ifelse(E(net.network)$weight >= net.quant[99],"black",
                             ifelse(E(net.network)$weight < net.quant[99] & E(net.network)$weight >= net.quant[95],"red",
                                    ifelse(E(net.network)$weight < net.quant[95] & E(net.network)$weight >= net.quant[90],"orange","blue")))

### set node size
V(net.network)$size <- degree(net.network)/.5

plot(net.network,layout=layout.circle(net.network),
     main="Firm Connectedness \n (intra day return volatility)",
     xlab = "black: 1st percentile \n red: 5th percentile \n orange: 10th percentile \n blue: 25th percentile \n node size: number of edeges connected to the node")
```

![](README_files/figure-markdown_github/unnamed-chunk-22-1.png)

### Volatility Index

Financial stocks are highly correlated, but it should come as a surprise that these variables are as connected as the index shows. However, the Connectedness Index is also very volatile here, falling from just over 80% to under 70% between 2015 - 2017, before jumping up above 80 again.

``` r
############## CONNECTEDNESS ################
vol.conn.index.out <- monetary_transmission[["conn"]]
### unlist output data
vol.conn.index.df <- data.frame(unlist(vol.conn.index.out))
colnames(vol.conn.index.df) = c("Index")
### .xts object
vol.conn.index.df$Date <- as.POSIXct(rownames(vol.conn.index.df),format="%Y-%m-%d")
vol.conn.index.xts <- xts(vol.conn.index.df[,-2],order.by=vol.conn.index.df[,2])

############## SPILLOVER ###################
fevd.list.out <- monetary_transmission[["spill"]]

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

![](README_files/figure-markdown_github/unnamed-chunk-23-1.png)

### Net Volatility Connectedness

The Net Volatility Connectedness and Spillovers show why the choice of assuming uncorrelated errors is not always wise. The prices of VCSH, FFRHX, and BSCK show hardly any signs of life, but the FEVD of these ETFs is high and volatile. In almost every other instance orthogonal errors assign an asset as giver (receiver) of volatility when the correlated errors assign the same asset as a receiver (giver). After considering these graphs, it is all the more surprising that both the Connectedness and Spillover indexes track each other as well as they do.

``` r
net.df.list <- monetary_transmission[["net_df"]]
spill.df.list <- monetary_transmission[["spillover_table"]]

net.cols <- names(net.df.list)
net.selected.net.df.list <- llply(net.df.list,function(df){df[["Net"]]})
net.selected.net.df <- as.data.frame(do.call(rbind,lapply(net.selected.net.df.list, function(x) t(data.frame(x)))))
net.selected.net.df <- net.selected.net.df %>% mutate(Date = as.POSIXct(net.cols,format="%Y-%m-%d"))

spill.cols <- names(spill.df.list)
net.selected.spill.df.list <- llply(spill.df.list,function(df){df[["Net"]]})
net.selected.spill.df <- as.data.frame(do.call(rbind,lapply(net.selected.spill.df.list, function(x) t(data.frame(x)))))
net.selected.spill.df <- net.selected.spill.df %>% mutate(Date = as.POSIXct(spill.cols,format="%Y-%m-%d"))

net.net.melt <- melt(net.selected.net.df,id="Date",value.name = "Connectedness")
net.spill.melt <- melt(net.selected.spill.df,id="Date",value.name = "Spillover")
net.tbl <- merge(net.net.melt,net.spill.melt,by=c("Date","variable"))
    
ggplot(net.tbl) + 
  geom_line(aes(x=Date,y=Connectedness),colour="black") +
  geom_line(aes(x=Date,y=Spillover),colour="red") +
  facet_wrap(~variable,scales="free_y") +
  theme(legend.position = "bottom")
```

![](README_files/figure-markdown_github/unnamed-chunk-24-1.png)
