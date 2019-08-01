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


class Settings_Caen:
	def __init__(self, infile='None'):
		self.infile = infile
		self.optlink = None
		self.node = None
		self.vme_b_addr = None
		self.use_usb = False
		self.wavedump_path = '/usr/local/bin'
		self.data_format = 'binary'
		self.dig_bits = 14
		self.points = 10240
		self.post_trig_percent = 95
		self.num_events = 10
		self.time_calib = 300
		self.test_name = 'muon_test_1'
		self.input_range = 2.0
		self.simultaneous_conversion = False
		self.plot_waveforms = False
		self.random_test = False
		self.time_res = 2e-9
		self.num_signals = 1
		self.sigCh = {}
		self.sig_thr = {}
		self.sig_polarity = {}
		self.sig_offset = {}
		self.trigCh = 7
		self.vetoCh = 2
		self.trig_polarity = -1
		self.veto_polarity = -1
		self.trig_thr = -0.1
		self.veto_thr = -0.1
		self.trig_offset = -45
		self.veto_offset = -45
		self.outdir = '.'
		self.prefix = 'waves'
		self.suffix = 'default'
		self.sigRes = 0
		self.UpdateSignalResolution()
		self.filename = ''

		self.struct_fmt = '@{p}H'.format(p=self.points)
		self.struct_len = struct.calcsize(self.struct_fmt)

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

				if parser.has_section('USB'):
					if parser.has_option('USB', 'use_usb'):
						self.use_usb = parser.getboolean('USB', 'use_usb')
					if parser.has_option('USB', 'wavedump_path'):
						self.wavedump_path = parser.get('USB', 'wavedump_path')

				if parser.has_section('RUN'):
					if parser.has_option('RUN', 'format'):
						self.data_format = parser.get('RUN', 'format').lower()
						if self.data_format not in ['binary', 'ascii']:
							ExitMessage('The format for the digitiser should be "binary" or "ascii"')
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
					if parser.has_option('RUN', 'test_name'):
						self.test_name = parser.get('RUN', 'test_name').lower()
					if parser.has_option('RUN', 'input_range'):
						self.input_range = parser.getfloat('RUN', 'input_range')
					if parser.has_option('RUN', 'simultaneous_conversion'):
						self.simultaneous_conversion = bool(parser.getboolean('RUN', 'simultaneous_conversion'))
					if parser.has_option('RUN', 'plot_waveforms'):
						self.plot_waveforms = bool(parser.getboolean('RUN', 'plot_waveforms'))
					if parser.has_option('RUN', 'random_test'):
						self.random_test = bool(parser.getboolean('RUN', 'random_test'))

				if parser.has_section('SIGNALS'):
					if parser.has_option('SIGNALS', 'number'):
						self.num_signals = parser.getint('SIGNALS', 'number')

				for sigi in xrange(self.num_signals):
					if parser.has_section('SIGNAL'+str(sigi)):
						if parser.has_option('SIGNAL'+str(sigi), 'caen_channel'):
							self.sigCh[sigi] = parser.getint('SIGNAL'+str(sigi), 'caen_channel')
						else:
							ExitMessage('Signal ' + str(sigi) + ' does not have a "caen_channel". Specify it under SIGNAL' + str(sigi), os.EX_CONFIG)
						if parser.has_option('SIGNAL'+str(sigi), 'polarity'):
							self.sig_polarity[sigi] = int(np.sign(parser.getfloat('SIGNAL'+str(sigi), 'polarity')))
						else:
							self.sig_polarity[sigi] = -1
						if parser.has_option('SIGNAL'+str(sigi), 'threshold'):
							self.sig_thr[sigi] = parser.getfloat('SIGNAL' + str(sigi), 'threshold')
						else:
							self.sig_thr[sigi] = 0.1 * self.sig_polarity[sigi]
						if parser.has_option('SIGNAL' + str(sigi), 'dc_offset'):
							self.sig_offset[sigi] = parser.getint('SIGNAL' + str(sigi), 'dc_offset')
						else:
							self.sig_offset[sigi] = self.sig_polarity[sigi] * 45

				if parser.has_section('TRIGGER'):
					if parser.has_option('TRIGGER', 'caen_channel'):
						self.trigCh = parser.getint('TRIGGER', 'caen_channel')
					if parser.has_option('TRIGGER', 'polarity'):
						self.trig_polarity = int(np.sign(parser.getfloat('TRIGGER', 'polarity')))
					else:
						self.trig_polarity = -1
					if parser.has_option('TRIGGER', 'threshold'):
						self.trig_thr = parser.getfloat('TRIGGER', 'threshold')
					else:
						self.trig_thr = self.trig_polarity * 0.1
					if parser.has_option('TRIGGER', 'dc_offset'):
						self.trig_offset = parser.getint('TRIGGER', 'dc_offset')
					else:
						self.trig_offset = self.trig_polarity * 45

				if parser.has_section('VETO'):
					if parser.has_option('VETO', 'caen_channel'):
						self.vetoCh = parser.getint('VETO', 'caen_channel')
					if parser.has_option('VETO', 'polarity'):
						self.veto_polarity = int(np.sign(parser.getfloat('VETO', 'polarity')))
					else:
						self.veto_polarity = -1
					if parser.has_option('VETO', 'thr_counts'):
						self.veto_thr = parser.getfloat('VETO', 'threshold')
					else:
						self.veto_thr = self.veto_polarity * 0.1
					if parser.has_option('VETO', 'dc_offset'):
						self.veto_offset = parser.getint('VETO', 'dc_offset')
					else:
						self.veto_offset = self.veto_polarity * 45

				if parser.has_section('OUTPUT'):
					if parser.has_option('OUTPUT', 'dir'):
						self.outdir = parser.get('OUTPUT', 'dir')
					if parser.has_option('OUTPUT', 'prefix'):
						self.prefix = parser.get('OUTPUT', 'prefix')
					if parser.has_option('OUTPUT', 'suffix'):
						self.suffix = parser.get('OUTPUT', 'suffix')
					else:
						self.suffix = ''
				self.UpdateSignalResolution()

	def UpdateSignalResolution(self):
		self.sigRes = np.divide(self.input_range, np.subtract(np.power(2, self.dig_bits, dtype='f8'), 1, dtype='f8'), dtype='f8')

	def SetOutputFiles(self):
		def AddSuffix(string1):
			if self.suffix != '':
				string1 += '_{s}'.format(s=self.suffix)
			return string1
		self.filename = '{d}_{p}'.format(p=self.prefix, d=self.test_name) if self.prefix != '' else self.test_name
		self.filename = AddSuffix(self.filename)

		if not os.path.isdir(self.outdir):
			os.makedirs(self.outdir)
		if not os.path.isdir('{d}/Runs'.format(d=self.outdir)):
			os.makedirs('{d}/Runs'.format(d=self.outdir))
		if not os.path.isdir('{d}/Runs/{f}'.format(d=self.outdir, f=self.filename)):
			os.makedirs('{d}/Runs/{f}'.format(d=self.outdir, f=self.filename))

	def SetupDigitiser(self, signal=None, trigger=None, veto=None):
		print 'Creating digitiser CAEN DT5730 configuration file... ', ; sys.stdout.flush()
		name_dest = '{d}/WaveDumpConfig_muon.txt'.format(d=self.outdir)

		cont = True
		with open('default/WaveDumpConfig_muon.txt', 'r') as source_file:
			with open(name_dest, 'w') as dest_file:
				for line in source_file:
					if cont:
						if line.startswith('OPEN PCI'):
							if not self.use_usb:
								dest_file.write('OPEN PCI {ol} {n} {ba}\n'.format(ol=int(self.optlink), n=int(self.node), ba=int(self.vme_b_addr)))
							else:
								dest_file.write('#OPEN PCI 1 0 32100000\n')
						elif line.startswith('OPEN USB'):
							if self.use_usb:
								dest_file.write('OPEN USB 0 0\n')
							else:
								dest_file.write('#OPEN USB 0 0\n')
						elif line.startswith('OUTPUT_FILE_FORMAT'):
							dest_file.write('OUTPUT_FILE_FORMAT\t{f}\n'.format(f=self.data_format.upper()))
						elif line.startswith('RECORD_LENGTH'):
							dest_file.write('RECORD_LENGTH\t{p}\n'.format(p=int(RoundInt(self.points))))
						elif line.startswith('POST_TRIGGER'):
							dest_file.write('POST_TRIGGER\t{pt}\n'.format(pt=int(round(self.post_trig_percent))))
							# dest_file.write('POST_TRIGGER\t{pt}\n'.format(pt=int(round(self.post_trig_percent*0.9996 - 1.6384))))
						elif line.startswith('CHANNEL_TRIGGER'):
							dest_file.write('CHANNEL_TRIGGER\tDISABLED\n')
							cont = False
						else:
							dest_file.write(line)

				dest_file.write('\n# configuration for each channel [0] to [15], although it only has 8 channels ;)')
				for ch in xrange(8):
					dest_file.write('\n\n[{ch}]'.format(ch=ch))
					if ch in self.sigCh.values():
						for sigi in xrange(self.num_signals):
							if self.sigCh[sigi] == ch:
								dest_file.write('\nENABLE_INPUT\tYES')
								dest_file.write('\nPULSE_POLARITY\t{sp}'.format(sp='POSITIVE' if self.sig_polarity[sigi] == 1 else 'NEGATIVE'))
								dest_file.write('\nDC_OFFSET\t{o}'.format(o=signal[sigi].dc_offset_percent))
								dest_file.write('\nCHANNEL_TRIGGER\tDISABLED')
					elif ch == self.trigCh:
						dest_file.write('\nENABLE_INPUT\tYES')
						dest_file.write('\nPULSE_POLARITY\t{sp}'.format(sp='POSITIVE' if self.trig_polarity == 1 else 'NEGATIVE'))
						dest_file.write('\nDC_OFFSET\t{o}'.format(o=trigger.dc_offset_percent))
						dest_file.write('\nCHANNEL_TRIGGER\tACQUISITION_ONLY')
						dest_file.write('\nTRIGGER_THRESHOLD\t{th}'.format(th=trigger.Volts_to_ADC(self.trig_thr)))
					elif ch == self.vetoCh:
						dest_file.write('\nENABLE_INPUT\tYES')
						dest_file.write('\nPULSE_POLARITY\t{sp}'.format(sp='POSITIVE' if self.veto_polarity == 1 else 'NEGATIVE'))
						dest_file.write('\nDC_OFFSET\t{o}'.format(o=veto.dc_offset_percent))
						dest_file.write('\nCHANNEL_TRIGGER\tDISABLED')
					else:
						dest_file.write('\nENABLE_INPUT\tNO')
				dest_file.write('\n')
		print 'Done'

	def ADC_to_Volts(self, adcs, channel):
		return channel.ADC_to_Volts(adcs, self.sigRes, self.dig_bits)

	def GetTriggerValueADCs(self, channel):
		return int(round(channel.base_line_u_adcs - channel.thr_counts - (2.0**self.dig_bits - 1) * (channel.dc_offset_percent/100.0 - 0.5)))

	def MoveBinaryFiles(self):
		print 'Moving binary files... ', ; sys.stdout.flush()
		for sigi in xrange(self.num_signals):
			shutil.move('raw_wave{chs}.dat'.format(chs=self.sigCh[sigi]), '{d}/Runs/{f}/{f}_signal{s}.dat'.format(s=sigi, d=self.outdir, f=self.filename))
		shutil.move('raw_wave{cht}.dat'.format(cht=self.trigCh), '{d}/Runs/{f}/{f}_trigger.dat'.format(d=self.outdir, f=self.filename))
		shutil.move('raw_wave{cha}.dat'.format(cha=self.vetoCh), '{d}/Runs/{f}/{f}_veto.dat'.format(d=self.outdir, f=self.filename))
		self.RemoveBinaries()
		print 'Done'

	def RenameDigitiserSettings(self):
		print 'Moving digitiser settings... ', ; sys.stdout.flush()
		shutil.move('{d}/WaveDumpConfig_muon.txt'.format(d=self.outdir), '{d}/Runs/{f}/WDConfig_muon_{f}.txt'.format(d=self.outdir, f=self.filename))
		print 'Done'

	def RemoveBinaries(self):
		channels = self.sigCh.values() + [self.trigCh] + [self.vetoCh]
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
