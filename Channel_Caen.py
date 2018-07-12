#!/usr/bin/env python
import visa
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ipdb
from pykeyboard import PyKeyboard
from ConfigParser import ConfigParser
import subprocess as subp
import struct
import ROOT as ro
import shutil
# from DataAcquisition import DataAcquisition


class Channel_Caen:
	def __init__(self, ch=3, ch_type='signal', verb=False):
		self.ch = ch
		self.type = ch_type
		self.verb = verb

		self.base_line_u_adcs = 0
		self.sigma_adcs = 10
		self.dc_offset_percent = 0
		self.thr_counts = 0
		self.edge = -1
		self.name = str(ch_type)

		self.offseted_adc_to_volts_cal = {'p0': 0, 'p1': np.divide(2.15, 2.0**14 - 1.0, dtype='float64')}

	def Set_Channel(self, settings):
		if self.type == 'signal_ch':
			self.edge = -int(settings.bias/abs(settings.bias)) if settings.bias != 0 else -1
			self.Calculate_DC_Offset_Percentage(settings)
			# Channels 3, 6 and 7 were calibrated with dc voltage and multimeter. The calibration files are on 20180419ch{X}/Runs
			if self.ch == 3:
				# With negative bias, signals are positive. With positive bias, signals are negative
				self.offseted_adc_to_volts_cal['p0'] = 0.035089408942074955 if settings.bias < 0 else -0.02328757136517118
				self.offseted_adc_to_volts_cal['p1'] = 0.00013089621340339722 if settings.bias < 0 else 0.00013076091653412987
			elif self.ch == 6:
				self.offseted_adc_to_volts_cal['p0'] = 0.04183464530415922 if settings.bias < 0 else -0.04426166525468815
				self.offseted_adc_to_volts_cal['p1'] = 0.0001294935440001848 if settings.bias < 0 else 0.00012934155437522722
			elif self.ch == 7:
				self.offseted_adc_to_volts_cal['p0'] = 0.025512742406801153 if settings.bias < 0 else -0.024895654896994378
				self.offseted_adc_to_volts_cal['p1'] = 0.00012875546036321804 if settings.bias < 0 else 0.00013155381396351944
			else:
				self.offseted_adc_to_volts_cal['p0'] = 0
				self.offseted_adc_to_volts_cal['p1'] = settings.sigRes
		elif self.type == 'trigger_ch':
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.trig_base_line, settings.sigRes)
			self.thr_counts = settings.trig_thr_counts
			self.Calculate_DC_Offset_Percentage(settings)
			self.edge = -1
		elif self.type == 'veto':
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.ac_base_line, settings.sigRes)
			self.thr_counts = settings.ac_thr_counts
			self.Calculate_DC_Offset_Percentage(settings)
			self.edge = -1

	def Calculate_DC_Offset_Percentage(self, settings):
		if self.type == 'signal_ch':
			self.dc_offset_percent = 45 if settings.bias < 0 else -45
		else:
			limit = self.base_line_u_adcs + 10 * self.sigma_adcs
			if limit < 0:
				self.dc_offset_percent = -50
			else:
				self.dc_offset_percent = int(round(100 * (limit/float(2**settings.dig_bits - 1) - 0.5)))
				self.dc_offset_percent = 48 if self.dc_offset_percent > 48 else -48 if self.dc_offset_percent < -48 else self.dc_offset_percent

	def Calculate_Universal_ADCs(self, value_volts, sig_res):
		return np.divide(value_volts, sig_res, dtype='f8')

	def Correct_Base_Line(self, mean_adc, sigma_adc, settings):
		self.sigma_adcs = sigma_adc
		variable = mean_adc + 50 * sigma_adc
		self.base_line_u_adcs = self.Calculate_Universal_ADCs(self.ADC_to_Volts(mean_adc, settings.sigRes, settings.dig_bits), settings.sigRes)
		if variable > (2**settings.dig_bits - 1) * (0.5 + self.dc_offset_percent / 100.0):
			self.dc_offset_percent += int(round(100.0 * (variable + 1 - 2.0**settings.dig_bits) / (2.0**settings.dig_bits - 1.0)))
		else:
			self.dc_offset_percent = -50
		self.dc_offset_percent = 48 if self.dc_offset_percent > 48 else -48 if self.dc_offset_percent < -48 else self.dc_offset_percent

	def Correct_Threshold(self, sigma):
		if self.type == 'trigger_ch':
			self.thr_counts = int(round(max(10*self.sigma_adcs, self.thr_counts)))
		elif self.type == 'veto':
			self.thr_counts = int(round(max(4*self.sigma_adcs, self.thr_counts)))

	def ADC_to_Volts(self, adcs, sigres, nbits=14):
		return np.add(self.offseted_adc_to_volts_cal['p0'], np.multiply(self.offseted_adc_to_volts_cal['p1'], np.add(adcs, np.multiply(2**nbits - 1, self.dc_offset_percent/100.0 - 0.5, dtype='float64'), dtype='float64'), dtype='float64'), dtype='float64')
		# return np.multiply(sigres, np.add(adcs, np.multiply(2**nbits - 1, self.dc_offset_percent/100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8')

if __name__ == '__main__':
	print 'bla'