#!/usr/bin/env python
import os
import shutil
import struct
import subprocess as subp
import sys
import time
from ConfigParser import ConfigParser
from optparse import OptionParser

import ROOT as ro
import numpy as np
import cPickle as pickle

from Channel_Caen import Channel_Caen
from Settings_Caen import Settings_Caen
from HV_Control import HV_Control
from Utils import *
from CCD_Caen import CCD_Caen
# from memory_profiler import profile

trig_rand_time = 0.2
wait_time_hv = 7

class VoltageScan:
	def __init__(self, infile='None', vini=5, vend=10, vstep=5, verbose=False):
		print 'Starting CCD program ...'
		self.infile = infile
		self.verb = verbose
		self.vini = vini
		self.vend = vend
		self.vstep = vstep
		self.voltages = np.linspace(self.vini, self.vend, int(round(float(self.vend - self.vini)/float(self.vstep)) + 1), dtype='int32')
		self.settings = Settings_Caen(self.infile, self.verb)
		self.settings.ReadInputFile()
		self.settings.Get_Calibration_Constants()
		self.settings.SetOutputFiles()
		self.wd = os.getcwd()

		self.p = {volt: None for volt in self.voltages}
		self.time0 = time.time()
		self.num_events = self.settings.num_events

	def DoVoltageScan(self):

		for volt in self.voltages:
			self.ResetSettings()
			self.settings.bias = float(volt)
			self.settings.Get_Calibration_Constants()
			self.settings.SetOutputFiles()

			self.p[volt] = CCD_Caen(settingsObj=self.settings)
			self.p[volt].StartHVControl()
			self.p[volt].GetBaseLines()
			self.p[volt].SavePickles()
			written_events = self.p[volt].GetData()
			self.settings.num_events = written_events
			self.p[volt].SavePickles()
			self.p[volt].settings.MoveBinaryFiles()
			self.p[volt].settings.RenameDigitiserSettings()
			self.p[volt].CloseHVClient()
			if not self.p[volt].settings.simultaneous_conversion:
				self.p[volt].CreateRootFile(files_moved=True)
				while self.p[volt].pconv.poll() is None:
					time.sleep(3)
				self.p[volt].CloseSubprocess('converter', stdin=False, stdout=False)

			print '\nFinished voltage scan for {v} Volts :)\n'.format(v=volt)
			sys.stdout.write('\a\a\a')
			sys.stdout.flush()

			self.p[volt] = None

		print 'Finished all voltage scan :D'
		sys.stdout.write('\a\a\a')
		sys.stdout.flush()

	def ResetSettings(self):
		self.settings = Settings_Caen(self.infile, self.verb)
		self.settings.ReadInputFile()
		self.settings.Get_Calibration_Constants()
		self.settings.SetOutputFiles()


def main():
	parser = OptionParser()
	parser.add_option('-i', '--infile', dest='infile', default='', type='string',
	                  help='Input configuration file. e.g. CAENCalibration.cfg')
	parser.add_option('--start', dest='start', type='int', help='Starting scan voltage')
	parser.add_option('--stop', dest='stop', type='int', help='Stopping scan voltage')
	parser.add_option('--step', dest='step', type='int', help='Voltage step between scans')
	# parser.add_option('--time', dest='time', type='int', help='maximum time in minutes before closing and starting next voltage')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic conversion and analysis afterwards', action='store_true')

	(options, args) = parser.parse_args()
	infile = str(options.infile)
	auto = bool(options.auto)
	verb = bool(options.verb)
	vini = int(options.start)
	vend = int(options.stop)
	vstep = int(options.step)
	# tau = int(options.time)

	vs = VoltageScan(infile, vini, vend, vstep, verb)
	if auto:
		vs.DoVoltageScan()

	return vs

if __name__ == '__main__':
	vs = main()
