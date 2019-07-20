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
from Utils import *
# from DataAcquisition import DataAcquisition

class Channel_Caen:
	def __init__(self, caen_ch=0, signal_number=-1, is_trigger=False, is_veto=False):
		self.caen_ch = caen_ch
		self.signal_number = signal_number
		self.is_trigger = is_trigger and signal_number == -1
		self.is_veto = is_veto and signal_number == -1
		self.dc_offset_percent = 0
		self.thr_volts = -0.1
		self.polarity = -1
		self.nbits = 14
		self.sigma_adcs = 10
		self.base_line_adcs = 0

	def Get_Cal_ADC_To_V(self):
		return np.divide(2, np.subtract(np.power(2, self.nbits, dtype='f8'), 1, dtype='f8'), dtype='f8')

	def Get_Cal_V_To_ADC(self):
		return np.divide(np.subtract(np.power(2, self.nbits, dtype='f8'), 1, dtype='f8'), 2, dtype='f8')

	def ADC_to_Volts(self, adc):
		return np.subtract(np.add(np.multiply(adc, self.Get_Cal_ADC_To_V(), dtype='f8'), np.divide(self.dc_offset_percent, 50.0, dtype='f8')), 1, dtype='f8')
		# return adc * self.Get_Cal_ADC_To_V() + self.dc_offset_percent / 50.0 - 1

	def Volts_to_ADC(self, volts):
		return RoundInt(np.multiply(self.Get_Cal_V_To_ADC(), np.add(np.subtract(volts, np.divide(self.dc_offset_percent, 50.0, dtype='f8'), dtype='f8'), 1, dtype='f8'), dtype='f8'))
		# return RoundInt(self.Get_Cal_V_To_ADC() * (volts - self.dc_offset_percent / 50.0 + 1))

	def Set_Channel(self, settings):
		# if it is a signal channel
		self.nbits = settings.dig_bits
		if not (self.is_trigger or self.is_veto):
			self.polarity = settings.sig_polarity[self.signal_number]
			self.dc_offset_percent = settings.sig_offset[self.signal_number]
			self.thr_volts = settings.sig_thr[self.signal_number]
		# if it is the trigger
		elif self.is_trigger:
			self.polarity = settings.trig_polarity
			self.dc_offset_percent = settings.trig_offset
			self.thr_volts = settings.trig_thr
		else:
			self.polarity = settings.veto_polarity
			self.dc_offset_percent = settings.veto_offset
			self.thr_volts = settings.veto_thr

		self.base_line_adcs = self.Volts_to_ADC(0)

if __name__ == '__main__':
	print 'bla'