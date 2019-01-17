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


class Cuts_Caen:
	def __init__(self, infile='None', outdir='None', filename='cal', sig=-1, trig=-1, meas=100, calv=False, atten=0.0, bias=0.05024, verb=False):
		self.infile = infile
		self.bias = np.double(bias * 10.0 ** (-self.atten/20.0))
		self.in_range = 2.15
		self.sigRes = np.double(np.divide(np.double(self.in_range), (np.power(2.0, 14.0, dtype='f8') - 1)))
		self.sigCh = 0
		self.trigCh = 1
		self.acCh = -1  # anti coincidence is not enabled unless the channel is specified in the config input file
		self.trigVal = 1
		self.verb = verb
		self.points = 2560
		self.thr_trigg = 35
		self.time_res = 2e-9
		self.post_trig_percent = 0.9
		self.num_events = 10


		self.wfmo, self.nrpt, self.xincr, self.xunit, self.xzero, self.ymult, self.yoffs,self.yunit, self.yzero = None, None, None, None, None, None, None, None, None
		# self.SetOutputFormatTwo()
		self.bindata1 = None
		self.bindata2 = None
		self.volts1, self.time1 = np.zeros(self.points, 'f8'), np.zeros(self.points, 'f8')
		self.volts2, self.time2 = np.zeros(self.points, 'f8'), np.zeros(self.points, 'f8')
		self.outString1 = ''
		self.outString2 = ''
		self.fileWaves = None
		self.fileWavesWriter = None
		self.iteration = 0
		# self.daq = DataAcquisition(self, self.verb)
		self.peak_values_waves = np.empty(0, 'f8')
		self.peak_values_measu = np.empty(0, 'f8')

		self.optlink = self.node = self.vme_b_addr = 0
		self.prefix, self.suffix = 'waves', ''
		self.wd_path = '/usr/local/bin/wavedump'

		self.ReadInputFile()
		self.struct_fmt = '@{p}H'.format(p=self.points)  # binary files with no header. Each event has self.points samples 2 bytes each (unsigned short)
		self.struct_len = struct.calcsize(self.struct_fmt)
		self.trig, self.signal = np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.timev = np.zeros(1, 'f8')
		if not self.calv:
			self.sig_offset = -45 if self.bias >= 0 else 45
		else:
			self.sig_offset = -38 if self.bias < 0 else 38
		self.trig_offset = 45

		self.rawFile, self.treeRaw = None, None
		midd = 'in' if self.calv else 'out'
		self.rawName = 'raw_{mi}_tree_cal_neg_{b}mV_waves'.format(mi=midd, b=abs(self.bias * 1000)) if self.bias < 0 else 'raw_{mi}_tree_cal_pos_{b}mV_waves'.format(mi=midd, b=(self.bias * 1000))
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
						self.wd_path = parser.get('OPTILINK', 'wavedump_path')
				if parser.has_section('RUN'):
					if parser.has_option('RUN', 'signal_channel') and self.sigCh == -1:
						self.sigCh = parser.getint('RUN', 'signal_channel')
					if parser.has_option('RUN', 'trig_channel') and self.trigCh == -1:
						self.trigCh = parser.getint('RUN', 'trig_channel')
					if parser.has_option('RUN', 'trig_val'):
						self.trigVal = parser.getfloat('RUN', 'trig_val')
					if parser.has_option('RUN', 'time'):
						self.points = int(np.ceil(parser.getfloat('RUN', 'time') * 1e-6 / self.time_res))
					if parser.has_option('RUN', 'post_trigger_percent'):
						self.post_trig_percent = parser.getint('RUN', 'post_trigger_percent')
				if parser.has_section('OUTPUT'):
					if parser.has_option('OUTPUT', 'inDir') and self.outdir == 'None':
						self.outdir = parser.get('OUTPUT', 'inDir')
					if parser.has_option('OUTPUT', 'prefix'):
						self.prefix = parser.get('OUTPUT', 'prefix')
					if parser.has_option('OUTPUT', 'suffix'):
						self.suffix = parser.get('OUTPUT', 'suffix')
			else:
				print 'Input file {f} does not exist. Loading default values...'.format(f=self.infile)
				self.LoadDefaults()
		else:
			self.LoadDefaults()

	def LoadDefaults(self):
		self.optlink = 1
		self.node = 0
		self.vme_b_addr = 32100000
		self.sigCh = 3 if self.sigCh == -1 else self.sigCh
		self.trigCh = 7 if self.trigCh == -1 else self.trigCh
		self.trigVal = 0.7138
		self.points = 5000
		self.post_trig_percent = 90

	def SetOutputFilesNames(self):
		if not os.path.isdir('{dir}/Runs'.format(dir=self.outdir)):
			os.makedirs('{dir}/Runs'.format(dir=self.outdir))

		def AddSuffix(string1):
			string1 += '_Input' if self.calv else '_Output'
			string1 += '_Pos' if self.bias >= 0 else '_Neg'
			string1 += '_{b}mV'.format(b=abs(1000 * self.bias))
			string1 += self.suffix
			return string1

		self.outString3 = '{p}_{f}'.format(dir=self.outdir, p=self.prefix, f=self.filename)
		self.outString3 = AddSuffix(self.outString3)
		self.filename = self.outString3

	def TakeTwoWaves(self):
		t0 = time.time()
		self.SetupDigitiser()
		print 'Starting getting data using wavedump...'
		p = subp.Popen(['wavedump', '{d}/WaveDumpConfig_CCD_cal.txt'.format(d=self.outdir)], bufsize=-1, stdin=subp.PIPE)
		t1 = time.time()
		# while p.poll() is None:
		self.Delay(4)
		p.stdin.write('c')
		p.stdin.flush()
		self.Delay(1.5)
		p.stdin.write('W')
		p.stdin.flush()
		self.Delay(1.5)
		p.stdin.write('P')
		p.stdin.flush()
		self.Delay(1.5)
		p.stdin.write('s')
		p.stdin.flush()
		while p.poll() is None:
			continue
		t0 = time.time() - t0
		print 'Total time saving {m} events:'.format(m=self.meas), t0, 'seconds'
		self.CreateRootFile()
		self.MoveBinaryFiles()

	def Delay(self, ti=1.0):
		t0 = time.time()
		while time.time() - t0 < ti:
			continue
		return

	def SetupDigitiser(self):
		print 'Creating digitiser CAEN V1730D configuration file... ', ; sys.stdout.flush()
		rfile = open('{d}/WaveDumpConfig_CCD_cal.txt'.format(d=self.outdir), 'w')
		rfile.write('[COMMON]')
		rfile.write('\n\n# open the digitezer')
		rfile.write('\nOPEN PCI {ol} {n} {ba}'.format(ol=int(self.optlink), n=int(self.node), ba=int(self.vme_b_addr)))
		rfile.write('\n\n# GNUPLOT path, normally /usr/bin/')
		rfile.write('\nGNUPLOT_PATH\t"/usr/bin/"')
		rfile.write('\n\n# output format can be BINARY or ASCII')
		rfile.write('\nOUTPUT_FILE_FORMAT\tBINARY')
		rfile.write('\n\n# if OUTPUT_FILE_HEADER is YES, the structure of the event has to be changed in self.struct_fmt to include the header')
		rfile.write('\nOUTPUT_FILE_HEADER NO')
		rfile.write('\n\n# specify the amount of samples to save. This defines the event window')
		rfile.write('\nRECORD_LENGTH\t{p}'.format(p=int(self.points)))
		rfile.write('\n\n# number of events to save in the file')
		rfile.write('\nMAX_NUM_EVENTS\t{n}'.format(n=int(self.meas)))
		rfile.write('\n\nTEST_PATTERN\tNO')
		rfile.write('\n\nENABLE_DES_MODE\tNO')
		rfile.write('\n\n# use external trigger. Options are: DISABLED, ACQUISITION_ONLY, ACQUISITION_AND_TRGOUT')
		rfile.write('\nEXTERNAL_TRIGGER\tDISABLED')
		rfile.write('\n\n# specify maximum number of events to read out in one Block Transfer. Must be between 1 and 1023')
		rfile.write('\nMAX_NUM_EVENTS_BLT\t100')
		rfile.write('\n\n# the percentage of the amount of data stored after the trigger in the event window. Has an offset of ~1.6... only accepts integers')
		rfile.write('\nPOST_TRIGGER\t{pt}'.format(pt=int(round(self.post_trig_percent*0.9996 - 1.6384))))
		rfile.write('\n\n# number of events that have to be ready before readout when the IRQ is asserted. 0 means run continuously. 1023 is the maximum')
		rfile.write('\nUSE_INTERRUPT\t0')
		rfile.write('\n\n# type of the fornt panel LEMO connectors: NIM, TTL')
		rfile.write('\n\nFPIO_LEVEL\tNIM')
		rfile.write('\n\nSKIP_STARTUP_CALIBRATION\tNO')
		rfile.write('\n\nCHANNEL_TRIGGER\tDISABLED')

		sig_polarity = 'POSITIVE' if self.bias >= 0 else 'NEGATIVE'

		rfile.write('\n\n# configuration for each channel [0] to [7]')
		for ch in xrange(16):
			rfile.write('\n\n[{ch}]'.format(ch=ch))
			if ch == self.sigCh or ch == self.trigCh:
				rfile.write('\nENABLE_INPUT\tYES')
			else:
				rfile.write('\nENABLE_INPUT\tNO')
			if ch == self.sigCh:
				rfile.write('\nPULSE_POLARITY\t{sp}'.format(sp=sig_polarity))
				rfile.write('\nDC_OFFSET\t{o}'.format(o=self.sig_offset))
				rfile.write('\nCHANNEL_TRIGGER\tDISABLED')
			if ch == self.trigCh:
				rfile.write('\nPULSE_POLARITY\tPOSITIVE')
				rfile.write('\nDC_OFFSET\t{o}'.format(o=self.trig_offset))
				rfile.write('\nCHANNEL_TRIGGER\tACQUISITION_ONLY')
				rfile.write('\nTRIGGER_THRESHOLD\t{th}'.format(th=int(round(np.divide(self.trigVal+1-self.trig_offset/50.0, self.sigRes)))))
		rfile.write('\n')
		rfile.close()
		print 'Done'

	def CreateRootFile(self):
		print 'Start creating root file'
		t0 = time.time()
		self.rawFile = ro.TFile('{d}/Runs/{r}.root'.format(d=self.outdir, r=self.filename), 'RECREATE')
		self.treeRaw = ro.TTree(self.filename, self.filename)
		fs = open('wave{s}.dat'.format(s=self.sigCh), 'rb')
		ft = open('wave{t}.dat'.format(t=self.trigCh), 'rb')
		eventBra = np.zeros(1, 'I')
		voltBra = np.zeros(self.points, 'f8')
		timeBra = np.zeros(self.points, 'f8')
		self.treeRaw.Branch('event', eventBra, 'event/i')
		self.treeRaw.Branch('time', timeBra, 'time[{s}]/D'.format(s=self.points))
		self.treeRaw.Branch('voltageSignal', voltBra, 'voltageSignal[{s}]/D'.format(s=self.points))
		self.CreateProgressBar(self.meas)
		self.bar.start()
		for ev in xrange(self.meas):
			fs.seek(ev * self.struct_len)
			ft.seek(ev * self.struct_len)
			datas = fs.read(self.struct_len)
			datat = ft.read(self.struct_len)
			if not datas or not datat:
				print 'No event in files... exiting'
				exit()
			s = struct.Struct(self.struct_fmt).unpack_from(datas)
			signalADCs = np.array(s, 'H')
			t = struct.Struct(self.struct_fmt).unpack_from(datat)
			triggADCs = np.array(t, 'H')
			signalVolts = np.array(np.multiply(signalADCs, self.sigRes) + self.sig_offset / 50.0 - 1, 'f8')
			triggVolts = np.array(np.multiply(triggADCs, self.sigRes) + self.trig_offset / 50.0 - 1, 'f8')
			left, right = np.double(triggVolts[0]), np.double(triggVolts[-1])
			mid = np.double((left + right) / 2.0)
			distFromMid = np.array(np.abs(triggVolts - mid), 'f8')
			midPos = distFromMid.argmin()
			timeVect = np.linspace(-midPos * self.time_res, self.time_res * (self.points - 1 - midPos), self.points, dtype='f8')
			eventBra.fill(ev)
			np.putmask(timeBra, 1 - np.zeros(self.points, '?'), timeVect)
			np.putmask(voltBra, 1 - np.zeros(self.points, '?'), signalVolts)
			numFil = self.treeRaw.Fill()
			self.bar.update(ev + 1)
		self.bar.finish()
		self.rawFile.Write()
		self.rawFile.Close()
		fs.close()
		ft.close()
		t0 = time.time() - t0
		print 'Time creating root tree:', t0, 'seconds'

	def MoveBinaryFiles(self):
		print 'Moving binary files... ', ; sys.stdout.flush()
		shutil.move('wave{chs}.dat'.format(chs=self.sigCh), '{d}/Runs/{f}_signal.dat'.format(d=self.outdir, f=self.filename))
		shutil.move('wave{cht}.dat'.format(cht=self.trigCh), '{d}/Runs/{f}_trigger.dat'.format(d=self.outdir, f=self.filename))
		print 'Done'

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
	print 'blaaa'
