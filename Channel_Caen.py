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
	def __init__(self, caen_ch=0, signal_number=-1, is_trigger=False, is_veto=False):
		self.caen_ch = caen_ch
		self.signal_number = signal_number
		self.is_trigger = is_trigger and signal_number == -1
		self.is_veto = is_veto and signal_number == -1
		self.base_line_u_adcs = 0
		self.sigma_adcs = 10
		self.dc_offset_percent = 0
		self.thr_counts = 0
		self.edge = -1

		self.offseted_adc_to_volts_cal = {'p0': 0, 'p1': np.divide(2.15, 2.0**14 - 1.0, dtype='float64')}

	def Set_Channel(self, settings):
		# if it is a scintillator signal
		self.edge = -1
		self.offseted_adc_to_volts_cal['p0'] = 0
		self.offseted_adc_to_volts_cal['p1'] = settings.sigRes
		if self.signal_number == 0:
			self.thr_counts = settings.signal_thr_counts_0
		elif self.signal_number == 1:
			self.thr_counts = settings.signal_thr_counts_1
		elif self.signal_number == 2:
			self.thr_counts = settings.signal_thr_counts_2
		elif self.signal_number == 3:
			self.thr_counts = settings.signal_thr_counts_3
		elif self.is_trigger:
			self.thr_counts = settings.trig_thr_counts
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.trig_base_line, settings.sigRes)
		elif self.is_veto:
			self.thr_counts = settings.veto_thr_counts
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.veto_base_line, settings.sigRes)
		self.Calculate_DC_Offset_Percentage(settings)

	def Calculate_DC_Offset_Percentage(self, settings):
		limit = self.base_line_u_adcs + 10 * self.sigma_adcs
		if limit < 0:
			self.dc_offset_percent = -50
		else:
			self.dc_offset_percent = int(round(100 * (limit/float(2**settings.dig_bits -1) - 0.5)))
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

	def Correct_Threshold(self):
		self.thr_counts = int(round(max(10*self.sigma_adcs, self.thr_counts)))

	def ADC_to_Volts(self, adcs, nbits=14):
		return np.add(self.offseted_adc_to_volts_cal['p0'], np.multiply(self.offseted_adc_to_volts_cal['p1'], np.add(adcs, np.multiply(2**nbits - 1, self.dc_offset_percent/100.0 - 0.5, dtype='float64'), dtype='float64'), dtype='float64'), dtype='float64')

if __name__ == '__main__':
	print 'bla'