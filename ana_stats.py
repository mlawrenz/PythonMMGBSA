import pandas, numpy, scipy
import random
import statsmodels.formula.api as sm
import pylab
import os
import sys
import operator
import optparse
from scipy import stats


def WLS(xdata, ydata, xerr):
    ws=pandas.DataFrame({'x':xdata, 'y':ydata})
    weights=pandas.Series(xerr)
    fit=sm.wls('y ~ x', data=ws, weights=1/weights).fit()
    Int, x=fit.pvalues
    residuals=fit.resid
    rval=fit.rsquared
    residuals=[abs(i) for i in residuals]
    newerr=numpy.sqrt(sum(residuals)/(len(residuals)-2))
    return round(rval, 2), fit.predict(), round(newerr,2)

def OLS(xdata, ydata):
    ws=pandas.DataFrame({'x':xdata, 'y':ydata})
    fit=sm.ols('y ~ x', data=ws).fit()
    Int, x=fit.pvalues
    residuals=fit.resid
    rval=fit.rsquared
    residuals=[abs(i) for i in residuals]
    newerr=numpy.sqrt(sum(residuals)/(len(residuals)-2))
    return round(rval, 2), fit.predict(), round(newerr,2)

def stats_results(xdata, ydata):
    slope, intercept, r_val, p_val, std_err=stats.linregress(xdata, ydata)
    print "R^2=%s" % round(r_val**2, 2)
    print "p=%s" % round(p_val, 4)
    print "coeff error =%s" % round(std_err, 2)
    prediction=slope*numpy.array(xdata)+intercept
    residuals=[(i-j)**2 for (i,j) in zip(ydata, prediction)]
    newerr=numpy.sqrt(sum(residuals)/(len(residuals)-2))
    print "residual std. err: %s" % round(newerr, 2)
    return r_val, prediction, newerr


def cum_avg(data):
    data=pandas.DataFrame({'data':data})
    means=pandas.expanding_mean(data)
    stds=pandas.expanding_std(data)
    return numpy.array([i[0] for i in means.values]), numpy.array([i[0] for i in stds.values])


def roll_avg(data, step):
    data=pandas.DataFrame({'data':data})
    out=pandas.rolling_mean(data, 2)
    return numpy.array([i[0] for i in out.values])


def bootstrap( b, n ):
  """
  Randomly divide n data points into b blocks.
  """
  s = [ random.randint( 0, n-1 ) for t in xrange(0, b) ]
  return s

def randomize(list):
    newlist=[]
    for i in list:
        ind=numpy.random.random_integers(0,len(list)-1)
        newlist.append(list[ind])
        list.pop(ind)
    return newlist

