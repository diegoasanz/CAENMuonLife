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


class Settings_Caen:
	def __init__(self, infile='None', verbose=False):
		self.infile = infile
		self.verb = verbose
		self.optlink = 1
		self.node = 0
		self.vme_b_addr = 32100000
		self.wavedump_path = '/usr/local/bin'
		self.dig_bits = 14
		self.points = 2560
		self.post_trig_percent = 90
		self.num_events = 10
		self.time_calib = 300
		self.dut = 'diamond'
		self.bias = 0
		self.input_range = 2.15
		self.calib_path = ''
		self.simultaneous_conversion = False
		self.plot_waveforms = False
		self.random_test = False
		self.time_res = 2e-9
		self.do_hv_control = False
		self.pics_folder_path = ''
		self.hv_supply = ''
		self.hv_ch = 0
		self.current_limit = 0
		self.hv_ramp = 10  # in V/s
		self.hot_start = True
		self.sigCh = 0
		self.trigCh = 1
		self.trig_base_line = -0.08
		self.trig_thr_counts = 35
		self.acCh = 2
		self.ac_base_line = -0.08
		self.ac_thr_counts = 15
		self.outdir = '.'
		self.prefix = 'waves'
		self.suffix = 'default'
		self.sigRes = 0
		self.UpdateSignalResolution()
		self.filename = ''

		self.fit_signal_vcal_params = np.array([-1.21888705e-04, -8.96215025e-01], dtype='f8')
		self.fit_signal_vcal_params_errors = np.array([0.0001923, 0.00559264], dtype='f8')
		self.fit_charge_signal_params = np.array([-1.65417060e+01, -1.35735246e+05], dtype='f8')
		self.fit_charge_signal_params_errors = np.array([26.10228586,  847.0342207], dtype='f8')
		self.fit_vcal_signal_params, self.fit_vcal_signal_params_errors = None, None
		self.UpdateVcalVsSignal()

		self.struct_fmt = '@{p}H'.format(p=self.points)
		self.struct_len = struct.calcsize(self.struct_fmt)

		self.hv_struct_fmt = '@IIIff' # struct for hv file: starting event is a uint, time in seconds is a uint, nanoseconds is a uint, voltage is float32, current is float32
		self.hv_struct_len = struct.calcsize(self.hv_struct_fmt)

		self.bar = None

	def ReadInputFile(self):
		parser = ConfigParser()
		if self.infile != 'None':
			if os.path.isfile(self.infile):
				print 'Reading input file: {f} ...'.format(f=self.infile)
				parser.read(self.infile)

				if parser.has_section('OPTILINK'):
					if parser.has_option('OPTILINK', 'link'):
						self.optlink = parser.getint('OPTILINK', 'link')
					if parser.has_option('OPTILINK', 'node'):
						self.node = parser.getint('OPTILINK', 'node')
					if parser.has_option('OPTILINK', 'vme_base_address'):
						self.vme_b_addr = parser.getint('OPTILINK', 'vme_base_address')
					if parser.has_option('OPTILINK', 'wavedump_path'):
						self.wavedump_path = parser.get('OPTILINK', 'wavedump_path')

				if parser.has_section('RUN'):
					if parser.has_option('RUN', 'time'):
						self.points = int(np.ceil(parser.getfloat('RUN', 'time') * 1.0e-6 / self.time_res))
						self.struct_fmt = '@{p}H'.format(p=self.points)
						self.struct_len = struct.calcsize(self.struct_fmt)
					if parser.has_option('RUN', 'post_trigger_percent'):
						self.post_trig_percent = parser.getint('RUN', 'post_trigger_percent')
					if parser.has_option('RUN', 'num_events'):
						self.num_events = parser.getint('RUN', 'num_events')
					if parser.has_option('RUN', 'time_calib'):
						self.time_calib = parser.getfloat('RUN', 'time_calib')
					if parser.has_option('RUN', 'dut'):
						self.dut = parser.get('RUN', 'dut').lower()
					if parser.has_option('RUN', 'sample_voltage'):
						self.bias = parser.getfloat('RUN', 'sample_voltage')
					if parser.has_option('RUN', 'input_range'):
						self.input_range = parser.getfloat('RUN', 'input_range')
					if parser.has_option('RUN', 'calib_path'):
						self.calib_path = parser.get('RUN', 'calib_path')
					if parser.has_option('RUN', 'simultaneous_conversion'):
						self.simultaneous_conversion = bool(parser.getboolean('RUN', 'simultaneous_conversion'))
					if parser.has_option('RUN', 'plot_waveforms'):
						self.plot_waveforms = bool(parser.getboolean('RUN', 'plot_waveforms'))
					if parser.has_option('RUN', 'random_test'):
						self.random_test = bool(parser.getboolean('RUN', 'random_test'))

				if parser.has_section('HV'):
					if parser.has_option('HV', 'path_Pics_folder'):
						self.pics_folder_path = parser.get('HV', 'path_Pics_folder')
					if parser.has_option('HV', 'HV_supply'):
						self.hv_supply = parser.get('HV', 'HV_supply')
						if self.hv_supply != '':
							self.do_hv_control = True
					if parser.has_option('HV', 'ch'):
						self.hv_ch = parser.getint('HV', 'ch')
					if parser.has_option('HV', 'current_limit'):
						self.current_limit = abs(parser.getfloat('HV', 'current_limit'))
					if parser.has_option('HV', 'ramp'):
						self.hv_ramp = abs(parser.getfloat('HV', 'ramp'))
					# TODO: implement option in HV_control for option hot_start = False
					# if parser.has_option('HV', 'hot_start'):
					# 	self.hot_start = bool(parser.getboolean('HV', 'hot_start'))

				if parser.has_section('SIGNAL'):
					if parser.has_option('SIGNAL', 'channel'):
						self.sigCh = parser.getint('SIGNAL', 'channel')

				if parser.has_section('TRIGGER'):
					if parser.has_option('TRIGGER', 'channel'):
						self.trigCh = parser.getint('TRIGGER', 'channel')
					if parser.has_option('TRIGGER', 'base_line'):
						self.trig_base_line = parser.getfloat('TRIGGER', 'base_line')
					if parser.has_option('TRIGGER', 'thr_counts'):
						self.trig_thr_counts = parser.getint('TRIGGER', 'thr_counts')

				if parser.has_section('ANTICOINCIDENCE'):
					if parser.has_option('ANTICOINCIDENCE', 'channel'):
						self.acCh = parser.getint('ANTICOINCIDENCE', 'channel')
					if parser.has_option('ANTICOINCIDENCE', 'base_line'):
						self.ac_base_line = parser.getfloat('ANTICOINCIDENCE', 'base_line')
					if parser.has_option('ANTICOINCIDENCE', 'thr_counts'):
						self.ac_thr_counts = parser.getint('ANTICOINCIDENCE', 'thr_counts')

				if parser.has_section('OUTPUT'):
					if parser.has_option('OUTPUT', 'inDir'):
						self.outdir = parser.get('OUTPUT', 'inDir')
					if parser.has_option('OUTPUT', 'prefix'):
						self.prefix = parser.get('OUTPUT', 'prefix')
					if parser.has_option('OUTPUT', 'suffix'):
						self.suffix = parser.get('OUTPUT', 'suffix')
					else:
						self.suffix = ''
				self.UpdateSignalResolution()

	def UpdateSignalResolution(self):
		self.sigRes = np.double(np.divide(np.double(self.input_range), (np.power(2.0, 14.0, dtype='f8') - 1)))

	def Get_Calibration_Constants(self):
		if self.calib_path == '':
			pass
		else:
			if os.path.isfile(self.calib_path):
				tempf = ro.TFile(self.calib_path, 'READ')
				fit_signal_vcal = tempf.Get('TFitResult-Signal_vs_Vcal-sig_vcal_fit')
				fit_charge_signal = tempf.Get('TFitResult-Charge_vs_Signal-q_sig_fit')
				self.fit_signal_vcal_params = np.array(fit_signal_vcal.Parameters(), dtype='f8')
				self.fit_charge_signal_params = np.array(fit_charge_signal.Parameters(), dtype='f8')
				self.UpdateVcalVsSignal()

	def UpdateVcalVsSignal(self):
		self.fit_vcal_signal_params = np.array([np.divide(-self.fit_signal_vcal_params[0], self.fit_signal_vcal_params[1], dtype='f8'),
		                                        np.divide(1.0, self.fit_signal_vcal_params[1], dtype='f8')], dtype='f8')
		self.fit_vcal_signal_params_errors = np.array([np.sqrt((self.fit_signal_vcal_params_errors[0] / self.fit_signal_vcal_params[1])**2 + (self.fit_signal_vcal_params[0] * self.fit_signal_vcal_params_errors[1] / (self.fit_signal_vcal_params[1] ** 2))**2)], dtype='f8')

	def SetOutputFiles(self):
		def AddSuffix(string1):
			string1 += '_Pos' if self.bias >= 0 else '_Neg'
			string1 += '_{b}V'.format(b=abs(self.bias))
			if self.suffix != '':
				string1 += '_{s}'.format(s=self.suffix)
			return string1

		self.filename = '{d}_{p}_ccd'.format(p=self.prefix, d=self.dut)
		self.filename = AddSuffix(self.filename)

		if not os.path.isdir(self.outdir):
			os.makedirs(self.outdir)
		if not os.path.isdir('{d}/Runs'.format(d=self.outdir)):
			os.makedirs('{d}/Runs'.format(d=self.outdir))
		if not os.path.isdir('{d}/Runs/{f}'.format(d=self.outdir, f=self.filename)):
			os.makedirs('{d}/Runs/{f}'.format(d=self.outdir, f=self.filename))

	def SetupDigitiser(self, doBaseLines=False, signal=None, trigger=None, ac=None, events_written=0):
		print 'Creating digitiser CAEN V1730D configuration file... ', ; sys.stdout.flush()
		name_dest = '{d}/WaveDumpConfig_CCD_BL.txt'.format(d=self.outdir) if doBaseLines else '{d}/WaveDumpConfig_CCD.txt'.format(d=self.outdir)

		sig_polarity = 'POSITIVE' if self.bias < 0 else 'NEGATIVE'
		cont = True
		with open('default/WaveDumpConfig_CCD.txt', 'r') as source_file:
			with open(name_dest, 'w') as dest_file:
				for line in source_file:
					if cont:
						if line.startswith('OPEN PCI'):
							dest_file.write('OPEN PCI {ol} {n} {ba}\n'.format(ol=int(self.optlink), n=int(self.node), ba=int(self.vme_b_addr)))
						elif line.startswith('RECORD_LENGTH'):
							dest_file.write('RECORD_LENGTH\t{p}\n'.format(p=int(self.points)))
						elif line.startswith('POST_TRIGGER'):
							dest_file.write('POST_TRIGGER\t{pt}\n'.format(pt=int(round(self.post_trig_percent*0.9996 - 1.6384))))
						elif line.startswith('CHANNEL_TRIGGER'):
							dest_file.write('CHANNEL_TRIGGER\tDISABLED\n')
							cont = False
						else:
							dest_file.write(line)

				dest_file.write('\n# configuration for each channel [0] to [15], although it only has 8 channels ;)')
				for ch in xrange(16):
					dest_file.write('\n\n[{ch}]'.format(ch=ch))
					if ch == signal.ch or ch == trigger.ch or ch == ac.ch:
						dest_file.write('\nENABLE_INPUT\tYES')
					else:
						dest_file.write('\nENABLE_INPUT\tNO')
					if ch == signal.ch:
						dest_file.write('\nPULSE_POLARITY\t{sp}'.format(sp=sig_polarity))
						dest_file.write('\nDC_OFFSET\t{o}'.format(o=signal.dc_offset_percent))
						dest_file.write('\nCHANNEL_TRIGGER\tDISABLED')
					elif ch == self.trigCh:
						dest_file.write('\nPULSE_POLARITY\tNEGATIVE')
						dest_file.write('\nDC_OFFSET\t{o}'.format(o=trigger.dc_offset_percent))
						if doBaseLines or self.random_test:
							dest_file.write('\nCHANNEL_TRIGGER\tDISABLED')
						else:
							dest_file.write('\nCHANNEL_TRIGGER\tACQUISITION_ONLY')
							dest_file.write('\nTRIGGER_THRESHOLD\t{th}'.format(th=self.GetTriggerValueADCs(trigger)))
					elif ch == ac.ch:
						dest_file.write('\nPULSE_POLARITY\tNEGATIVE')
						dest_file.write('\nDC_OFFSET\t{o}'.format(o=ac.dc_offset_percent))
						dest_file.write('\nCHANNEL_TRIGGER\tDISABLED')
				dest_file.write('\n')
		print 'Done'

	def ADC_to_Volts(self, adcs, channel):
		return channel.ADC_to_Volts(adcs, self.sigRes, self.dig_bits)

	def GetTriggerValueADCs(self, channel):
		return int(round(channel.base_line_u_adcs - channel.thr_counts - (2.0**self.dig_bits - 1) * (channel.dc_offset_percent/100.0 - 0.5)))

	def MoveBinaryFiles(self):
		print 'Moving binary files... ', ; sys.stdout.flush()
		shutil.move('raw_wave{chs}.dat'.format(chs=self.sigCh), '{d}/Runs/{f}/{f}_signal.dat'.format(d=self.outdir, f=self.filename))
		shutil.move('raw_wave{cht}.dat'.format(cht=self.trigCh), '{d}/Runs/{f}/{f}_trigger.dat'.format(d=self.outdir, f=self.filename))
		shutil.move('raw_wave{cha}.dat'.format(cha=self.acCh), '{d}/Runs/{f}/{f}_veto.dat'.format(d=self.outdir, f=self.filename))
		self.RemoveBinaries()
		print 'Done'

	def RenameDigitiserSettings(self):
		print 'Moving digitiser settings... ', ; sys.stdout.flush()
		shutil.move('{d}/WaveDumpConfig_CCD.txt'.format(d=self.outdir), '{d}/Runs/{f}/WDConfig_CCD_{f}.txt'.format(d=self.outdir, f=self.filename))
		print 'Done'

	def RemoveBinaries(self):
		channels = [self.sigCh, self.trigCh, self.acCh]
		for ch in channels:
			if os.path.isfile('wave{c}.dat'.format(c=ch)):
				os.remove('wave{c}.dat'.format(c=ch))

	def CreateProgressBar(self, maxVal=1):
		widgets = [
			'Processed: ', progressbar.Counter(),
			' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			#  ' ', progressbar.AdaptativeTransferSpeed()
		]
		self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)


if __name__ == '__main__':
	print 'bla'
