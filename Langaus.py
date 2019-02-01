#!/usr/bin/env python

import ROOT as ro
import numpy as np
import ipdb

class LanGaus:
	def __init__(self, histo):
		self.params = np.zeros(4, 'f8')
		self.conv_steps = 1000
		self.sigma_conv = 5
		self.mpshift = -0.22278298
		self.paramErrors = np.zeros(4, 'f8')
		self.paramsLimitsLow = np.zeros(4, 'f8')
		self.paramsLimitsHigh = np.zeros(4, 'f8')
		self.paramsFitErrors = np.zeros(4, 'f8')
		self.histo = histo
		self.chi2 = 0
		self.ndf = 0
		# self.fit_range = [max(self.histo.GetMean()*0.3, 0), min(self.histo.GetXaxis().GetXmax(), self.histo.GetMean()*3)]
		self.fit_range = {'min': max(self.histo.GetMean()*0.2, 0), 'max': min(self.histo.GetXaxis().GetXmax(), self.histo.GetMean()*3)}
		# self.fit_range = [0, min(self.histo.GetXaxis().GetXmax(), self.histo.GetMean()*5)]
		self.fit = None
		self.EstimateParameters()
		self.EstimateParametersLimits()

	def LangausFunc(self, x, params):
		mpc = params[1] - self.mpshift * params[0]
		xlow, xhigh = [x[0] + self.sigma_conv * i * params[3] for i in [-1, 1]]
		step = (xhigh - xlow) / self.conv_steps
		sums = 0
		for i in xrange(1, int(np.ceil(self.conv_steps/2.0 + 1))):
			xx = xlow + (i - 0.5) * step
			fland = ro.TMath.Landau(xx, mpc, params[0]) / params[0]
			sums += fland * ro.TMath.Gaus(x[0], xx, params[3])
			xx = xhigh - (i - 0.5) * step
			fland = ro.TMath.Landau(xx, mpc, params[0]) / params[0]
			sums += fland * ro.TMath.Gaus(x[0], xx, params[3])
		return params[2] * step * sums / (np.sqrt(2 * np.pi, dtype='f8') * params[3])

	def EstimateParameters(self):
		self.params[0] = self.histo.GetRMS()/5.7
		self.params[1] = self.histo.GetBinCenter(self.histo.GetMaximumBin())
		self.params[2] = self.histo.Integral() * 500  # 10
		self.params[3] = self.histo.GetRMS() / 3.4

	def EstimateParametersLimits(self):
		self.paramsLimitsLow[0] = self.histo.GetRMS()/15.0
		self.paramsLimitsLow[1] = max(0, self.histo.GetBinCenter(self.histo.GetMaximumBin()) - 3*self.params[0])
		self.paramsLimitsLow[2] = self.histo.Integral() / 5000.0
		self.paramsLimitsLow[3] = self.histo.GetRMS() / 5.0
		self.paramsLimitsHigh[0] = self.histo.GetRMS()/4.0
		self.paramsLimitsHigh[1] = min(self.histo.GetXaxis().GetXmax(), self.histo.GetBinCenter(self.histo.GetMaximumBin()) + 3*self.params[0])
		self.paramsLimitsHigh[2] = self.histo.Integral() * 5000
		self.paramsLimitsHigh[3] = self.histo.GetRMS()

	def LanGausFit(self, nconv=1000, doLikelihood=False, xmin=-10000000, xmax=-10000000):
		fit_name = 'fit_{n}'.format(n=self.histo.GetName())
		self.conv_steps = nconv

		fit_old = ro.gROOT.GetListOfFunctions().FindObject(fit_name)
		if fit_old:
			fit_old.Delete()
			del fit_old
		self.fit_range['min'] = xmin if xmin > -10000000 else self.histo.GetBinLowEdge(self.histo.FindFirstBinAbove(0))
		self.fit_range['max'] = xmax if xmax > -10000000 else self.histo.GetBinLowEdge(self.histo.FindLastBinAbove(0) + 1)
		self.fit = ro.TF1(fit_name, self.LangausFunc, self.histo.GetXaxis().GetXmin(), self.histo.GetXaxis().GetXmax(), 4)
		self.fit.SetNpx(1000)
		self.fit.SetParameters(self.params)
		self.fit.SetParNames('Width', 'MP', 'Area', 'GSigma')
		options = 'QB0ML' if doLikelihood else 'QB0M'
		ro.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')
		for i in xrange(len(self.params)):
			self.fit.SetParLimits(i, self.paramsLimitsLow[i], self.paramsLimitsHigh[i])
		self.histo.Fit(fit_name, options, '', self.fit_range['min'], self.fit_range['max'])
		self.fit.GetParameters(self.params)
		for i in xrange(len(self.params)):
			self.paramsFitErrors[i] = self.fit.GetParError(i)
		self.chi2 = self.fit.GetChisquare()
		self.ndf = self.fit.GetNDF()

if __name__ == '__main__':
	mu = 1000
	sigma = 20
	gaus_sigma = 5
	rand = ro.TRandom3()
	rand.SetSeed(0)
	datad = np.array([rand.Landau(mu, sigma) + rand.Gaus(0, gaus_sigma) for i in xrange(1000)])
	datah = ro.TH1F('datah', 'datah', 200, 0, 4000)
	for point in datad:
		datah.Fill(point)
	lg = LanGaus(histo=datah)
	lg.LanGausFit(nconv=100)
