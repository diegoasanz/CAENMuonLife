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
import pickle
import shutil
from copy import deepcopy


# from DataAcquisition import DataAcquisition

class Converter_Caen:
	def __init__(self, settings_object=''):
		self.settings_object = settings_object

		self.settings_full_path = os.path.abspath(settings_object)
		self.working_dir = '/'.join(self.settings_full_path.split('/')[:-1])
		self.filename = self.settings_full_path.split('/')[-1].split('.settings')[0]
		self.settings = pickle.load(open('{d}/{f}.settings'.format(d=self.working_dir, f=self.filename), 'rb'))
		self.signal_ch = pickle.load(open('{d}/{f}.signal'.format(d=self.working_dir, f=self.filename), 'rb'))
		self.trigger_ch = pickle.load(open('{d}/{f}.trigger'.format(d=self.working_dir, f=self.filename), 'rb'))
		self.veto_ch = pickle.load(open('{d}/{f}.veto'.format(d=self.working_dir, f=self.filename), 'rb'))

		self.signal_path = self.working_dir + '/' + self.filename + '_signal.dat'
		self.trigger_path = self.working_dir + '/' + self.filename + '_trigger.dat'
		self.veto_path = self.working_dir + '/' + self.filename + '_veto.dat'
		self.points = self.settings.points
		self.num_events = self.settings.num_events
		self.struct_len = self.settings.struct_len
		self.struct_fmt = self.settings.struct_fmt
		self.adc_res = self.settings.sigRes
		self.sig_offset = self.signal_ch.dc_offset_percent
		self.trig_offset = self.trigger_ch.dc_offset_percent
		self.anti_co_offset = self.veto_ch.dc_offset_percent
		self.time_res = self.settings.time_res
		self.post_trig_percent = self.settings.post_trig_percent
		self.trig_value = self.settings.ADC_to_Volts(self.settings.GetTriggerValueADCs(self.trigger_ch), self.trigger_ch)
		self.veto_value = self.veto_ch.thr_counts
		self.dig_bits = self.settings.dig_bits
		self.simultaneous_conversion = self.settings.simultaneous_conversion
		self.time_recal = self.settings.time_calib
		self.control_hv = self.settings.do_hv_control
		self.polarity = 1 if self.settings.bias >= 0 else -1

		self.hv_file_name = 'hvfile_{f}.dat'.format(f=self.filename)
		self.hv_dict = None

		self.doVeto = True
		self.array_points = np.arange(self.points, dtype=np.dtype('int32'))

		self.raw_file = None
		self.raw_tree = None
		self.eventBra = self.voltBra = self.trigBra = self.vetoBra = self.timeBra = self.vetoedBra = self.badShapeBra = self.badPedBra = None
		self.hvVoltageBra = self.hvCurrentBra = None
		# self.hourBra = self.minuteBra = self.secondBra = None
		self.hourMinSecBra = None

		self.t0 = time.time()

		self.signal_written_events = self.trigger_written_events = self.anti_co_written_events = None
		self.fs = self.ft = self.fa = None
		self.wait_for_data = None

		self.datas = self.datat = self.dataa = None
		self.sigADC = self.trigADC = self.vetoADC = None
		self.sigVolts = self.trigVolts = self.vetoVolts = None
		self.trigPos = None
		self.timeVect = None
		self.vetoed_event = None
		self.bad_shape_event = None
		self.bad_pedstal_event = None
		self.condition_base_line = None
		self.condition_peak_pos = None
		self.hv_voltage_event = None
		self.hv_current_event = None
		self.hour_event = self.minute_event = self.second_event = None
		self.time_break = int(np.ceil(self.time_recal + 30))
		self.file_hv = None
		# self.hour_min_sec_event = None
		# self.currentTime = None

		self.struct_s = self.struct_t = self.struct_ac = None

		self.bar = None

	def SetupRootFile(self, arguments):
		if self.simultaneous_conversion:
			print '\n', str(arguments)
			print 'Start creating root file simultaneously with data taking'
		else:
			print str(arguments)
			print 'Start creating root file'
		self.raw_file = ro.TFile('{wd}/{d}/Runs/{r}.root'.format(wd=self.working_dir_location, d=self.run_dir_location, r=self.filename), 'RECREATE')
		self.raw_tree = ro.TTree(self.filename, self.filename)
		self.raw_tree.SetAutoFlush(100)
		self.raw_tree.SetAutoSave(-10485760)
		if self.doVeto:
			self.vetoBra = np.zeros(self.points, 'f8')
			self.vetoedBra = np.zeros(1, '?')
		self.eventBra = np.zeros(1, 'I')
		self.voltBra = np.zeros(self.points, 'f8')
		self.trigBra = np.zeros(self.points, 'f8')
		self.timeBra = np.zeros(self.points, 'f8')
		self.badShapeBra = np.zeros(1, dtype=np.dtype('int8'))  # signed char
		self.badPedBra = np.zeros(1, '?')
		if self.control_hv:
			self.hvVoltageBra = np.zeros(1, 'f')
			self.hvCurrentBra = np.zeros(1, 'f')
		# self.hourBra = np.zeros(1, dtype=np.dtype('uint8'))  # unsigned char
		# self.minuteBra = np.zeros(1, dtype=np.dtype('uint8'))  # unsigned char
		# self.secondBra = np.zeros(1, 'B')  # unsigned char
		self.hourMinSecBra = ro.TDatime()
		self.raw_tree.Branch('event', self.eventBra, 'event/i')
		self.raw_tree.Branch('time', self.timeBra, 'time[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageSignal', self.voltBra, 'voltageSignal[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageTrigger', self.trigBra, 'voltageTrigger[{s}]/D'.format(s=self.points))
		if self.doVeto:
			self.raw_tree.Branch('voltageVeto', self.vetoBra, 'voltageVeto[{s}]/D'.format(s=self.points))
			self.raw_tree.Branch('vetoedEvent', self.vetoedBra, 'vetoedEvent/O')
		self.raw_tree.Branch('badShape', self.badShapeBra, 'badShape/B')  # signed char
		self.raw_tree.Branch('badPedestal', self.badPedBra, 'badPedestal/O')
		if self.control_hv:
			self.raw_tree.Branch('voltageHV', self.hvVoltageBra, 'voltageHV/F')
			self.raw_tree.Branch('currentHV', self.hvCurrentBra, 'currentHV/F')
		# self.raw_tree.Branch('hour', self.hourBra, 'hour/b')
		# self.raw_tree.Branch('minute', self.minuteBra, 'minute/b')
		# self.raw_tree.Branch('second', self.secondBra, 'second/b')
		self.raw_tree.Branch('timeHV', self.hourMinSecBra)

	def OpenRawBinaries(self):
		self.fs = open('{wd}/{s}'.format(wd=self.working_dir_location, s=self.signal_path), 'rb')
		self.ft = open('{wd}/{t}'.format(wd=self.working_dir_location, t=self.trigger_path), 'rb')
		if self.doVeto:
			self.fa = open('{wd}/{a}'.format(wd=self.working_dir_location, a=self.anti_co_path), 'rb')

	def GetBinariesNumberWrittenEvents(self):
		self.signal_written_events = int(round(os.path.getsize('{d}/{s}'.format(d=self.working_dir_location, s=self.signal_path)) / self.struct_len))
		self.trigger_written_events = int(round(os.path.getsize('{d}/{t}'.format(d=self.working_dir_location, t=self.trigger_path)) / self.struct_len))
		if self.doVeto:
			self.anti_co_written_events = int(round(os.path.getsize('{d}/{a}'.format(d=self.working_dir_location, a=self.anti_co_path)) / self.struct_len))

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

	def CheckFilesSizes(self, ev):
		self.wait_for_data = True if (self.signal_written_events <= ev or self.trigger_written_events <= ev) else False
		if self.doVeto:
			self.wait_for_data = bool(self.wait_for_data or (self.anti_co_written_events <= ev))

	def WaitForData(self, ev, t1):
		while self.wait_for_data:
			if self.simultaneous_conversion:
				if time.time() - t1 > self.time_break:
					print 'No data has been saved in file for event {ev} in the past {t} seconds... exiting!'.format(ev=ev, t=self.time_break)
					exit()
				if not self.fs.closed:
					self.fs.close()
				if not self.ft.closed:
					self.ft.close()
				if self.doVeto:
					if not self.fa.closed:
						self.fa.close()
				self.GetBinariesNumberWrittenEvents()
				self.CheckFilesSizes(ev)
				if not self.wait_for_data:
					self.OpenRawBinaries()
			else:
				print 'The data is corrupted... exiting'
				exit()

	def ReadData(self, ev):
		self.fs.seek(ev * self.struct_len, 0)
		self.datas = self.fs.read(self.struct_len)
		self.ft.seek(ev * self.struct_len, 0)
		self.datat = self.ft.read(self.struct_len)
		if self.doVeto:
			self.fa.seek(ev * self.struct_len, 0)
			self.dataa = self.fa.read(self.struct_len)
		if self.control_hv:
			self.Read_HV_File()
		# self.currentTime = time.gmtime()

	def CheckData(self):
		if not self.datas or not self.datat:
			print 'No event in signal or trigger files... exiting'
			exit()
		elif self.doVeto:
			if not self.dataa:
				print 'No event in veto trigger file... exiting'
				exit()

	def FillBranches(self, ev):
		self.eventBra.fill(ev)
		np.putmask(self.timeBra, np.bitwise_not(np.zeros(self.points, '?')), self.timeVect)
		np.putmask(self.voltBra, np.bitwise_not(np.zeros(self.points, '?')), self.sigVolts)
		np.putmask(self.trigBra, np.bitwise_not(np.zeros(self.points, '?')), self.trigVolts)
		if self.doVeto:
			np.putmask(self.vetoBra, np.bitwise_not(np.zeros(self.points, '?')), self.vetoVolts)
			self.vetoedBra.fill(self.vetoed_event)
		self.badShapeBra.fill(self.bad_shape_event)
		self.badPedBra.fill(self.bad_pedstal_event)
		if self.control_hv:
			self.hvVoltageBra.fill(self.hv_voltage_event)
			self.hvCurrentBra.fill(self.hv_current_event)
		# self.hourBra.fill(self.currentTime[3])
		# self.minuteBra.fill(self.currentTime[4])
		# self.secondBra.fill(self.currentTime[5])
		self.hourMinSecBra.Set()

	def ConvertEvents(self):
		self.bar.start()
		# blat = time.time()
		# while time.time() - blat < 15:
		# 	pass
		for ev in xrange(self.num_events):
			t1 = time.time()
			self.CheckFilesSizes(ev)
			self.WaitForData(ev, t1)
			self.ReadData(ev)
			self.CheckData()
			self.struct_s = struct.Struct(self.struct_fmt).unpack_from(self.datas)
			self.sigADC = np.array(self.struct_s, 'H')
			self.sigVolts = self.ADC_to_Volts('signal')
			self.struct_t = struct.Struct(struct_fmt).unpack_from(self.datat)
			self.trigADC = np.array(self.struct_t, 'H')
			self.trigVolts = self.ADC_to_Volts('trigger')
			self.LookForTime0()
			self.timeVect = np.linspace(-self.trigPos * self.time_res, self.time_res * (self.points - 1 - self.trigPos), self.points, dtype='f8')
			if self.doVeto:
				self.struct_ac = struct.Struct(struct_fmt).unpack_from(self.dataa)
				self.vetoADC = np.array(self.struct_ac, 'H')
				self.vetoVolts = self.ADC_to_Volts('veto')
				self.vetoed_event = self.IsEventVetoed()
			self.DefineSignalBaseLineAndPeakPosition()
			self.bad_shape_event = self.IsEventBadShape()
			self.bad_pedstal_event = self.IsPedestalBad()
			self.FillBranches(ev)
			if ev == 10:
				self.raw_tree.OptimizeBaskets()
			numFil = self.raw_tree.Fill()
			self.bar.update(ev + 1)

	def CloseAll(self):
		self.bar.finish()
		self.raw_file.Write()
		self.raw_file.Close()
		self.fs.close()
		del self.fs
		self.ft.close()
		del self.ft
		if self.doVeto:
			self.fa.close()
			del self.fa
		self.t0 = time.time() - self.t0
		print 'Time creating root tree:', self.t0, 'seconds'
		exit()

	def LookForTime0(self):
		guess_pos = int(round(self.points * (100.0 - self.post_trig_percent)/100.0))
		condition_trigg = np.array(np.abs(self.array_points - guess_pos) <= int(round(0.1e-6/self.time_res)), dtype='?')
		condition_no_trigg = np.array(1 - condition_trigg, dtype='?')
		mean = np.extract(condition_no_trigg, self.trigVolts).mean()
		sigma = np.extract(condition_no_trigg, self.trigVolts).std()
		temp_trig_volts = np.copy(self.trigVolts)
		np.putmask(temp_trig_volts, condition_no_trigg, 100)
		volt_min_pos = temp_trig_volts.argmin()
		condition_trigg = np.bitwise_and(condition_trigg, np.array(self.array_points <= volt_min_pos))
		np.putmask(temp_trig_volts, np.bitwise_not(condition_trigg), 100)
		self.trigPos = np.abs(temp_trig_volts - self.trig_value).argmin()

	def IsEventVetoed(self):
		window_around_trigg = 50e-9
		condition_veto_base_line = np.array(np.abs(self.array_points - self.trigPos) > int(round(window_around_trigg/float(self.time_res))), dtype='?')
		condition_search = np.array(1 - condition_veto_base_line, dtype='?')
		mean = np.extract(condition_veto_base_line, self.vetoADC).mean()
		sigma = np.extract(condition_veto_base_line, self.vetoADC).std()
		vetoValNew = 4 * sigma if self.veto_value < 0.9 * 4 * sigma else self.veto_value
		veto_event = bool((np.extract(condition_search, self.vetoADC) - mean + vetoValNew).min() <= 0)
		return veto_event

	def IsPedestalBad(self):
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		mean = np.extract(self.condition_base_line, self.sigADC).mean()
		sigma_volts = sigma * self.adc_res
		mean_volts = mean * self.adc_res
		diff_volts = abs(int(self.sigADC[0]) - int(self.sigADC[self.trigPos])) * self.adc_res
		if sigma_volts >= 2e-3 or diff_volts >= 15e-3:  # if a signal is not flat enough due to previous unstable states
			return True
		else:
			last_sig = self.sigADC[-1] * self.adc_res
			if abs(last_sig - mean_volts) > 75e-3:  # if a signal does not return to base line due to multiple effects
				return True
			else:  # when the signal pedestal is well behaved and the signal returns to baseline
				return False

	def DefineSignalBaseLineAndPeakPosition(self):
		self.condition_base_line = np.array(self.array_points <= self.trigPos, dtype='?')
		# values 2.2 and 0.5 us come from shape of signal
		self.condition_peak_pos = np.array(np.abs(self.array_points - (2.2e-6/float(self.time_res) + self.trigPos)) <= 0.5e-6/float(self.time_res), dtype='?')

	def IsEventBadShape(self):
		# mean = np.extract(self.condition_base_line, self.sigADC).mean()
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		lim_inf = self.condition_peak_pos.argmax()
		lim_sup = self.points - self.condition_peak_pos[::-1].argmax() - 1
		peak_pos = self.sigADC.argmin() if self.polarity == 1 else self.sigADC.argmax()
		if lim_inf < peak_pos < lim_sup:
			# The event has a good shape
			return 0
		else:
			modified_adc = self.sigADC - sigma if self.polarity == 1 else self.sigADC + sigma
			modified_adc[lim_inf] += 2*sigma if self.polarity == 1 else -2*sigma
			modified_adc[lim_sup] += 2*sigma if self.polarity == 1 else -2*sigma
			peak_pos = modified_adc.argmin() if self.polarity == 1 else modified_adc.argmax()
			if lim_inf < peak_pos < lim_sup:
				# Can't tell if the event has a bad shape
				return -1
			else:
				# Event has bad shape
				return 1

	def IsPedestalNotFlat(self, signalADC, points, trigPos, time_res):
		array_points = np.arange(points, dtype=np.dtype('int32'))
		condition_base_line = np.array(array_points - trigPos <= 0, dtype='?')

	def ADC_to_Volts(self, sig_type):
		if sig_type == 'signal':
			adcs = self.sigADC
			offset = self.sig_offset
		elif sig_type == 'trigger':
			adcs = self.trigADC
			offset = self.trig_offset
		elif sig_type == 'veto':
			adcs = self.vetoADC
			offset = self.anti_co_offset
		else:
			print 'Wrong type. Exiting'
			exit()
		result = np.multiply(self.adc_res, np.add(adcs, np.multiply(2 ** self.dig_bits - 1.0, offset / 100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8')
		return result

	def Read_HV_File(self):
		line = [0, 0]
		temp_line = ''
		if self.simultaneous_conversion:
			if os.path.isfile(self.hv_file_name):
				with open(self.hv_file_name, 'r') as self.file_hv:
					temp_line = self.file_hv.readline()
					line = temp_line.split() if temp_line != '' else line
		self.hv_voltage_event, self.hv_current_event = float(line[0]), float(line[1])
		del line, temp_line, self.file_hv

if __name__ == '__main__':
	settings_object = str(sys.argv[2])  # settings pickle path

	converter = Converter_Caen(settings_object=settings_object)
	# run_dir_location = str(sys.argv[1])  # output directory location
	# working_dir_location = str(sys.argv[2])  # current directory location, where the raw_wave files are
	# filename = str(sys.argv[3])  # file name for the files
	# signal_path = str(sys.argv[4])  # caen channel for the ccd signal
	# trigger_path = str(sys.argv[5])  # caen channel for the trigger signal
	# anti_co_path = str(sys.argv[6])  # caen channel for the veto scintillator signal if no anti coincidence then this should be ''
	# points = int(sys.argv[7])  # caen number of points in each event i.e. time window devided by time resolution
	# num_events = int(sys.argv[8])  # number of events to convert
	# struct_len = int(sys.argv[9])  # structure length (obtained from struct_fmt and struct.calcsize(self.struct_fmt). e.g. 5120
	# struct_fmt = str(sys.argv[10])  # structure format of the data per event i.e. '@{p}H'.format(p=self.points). e.g. '@2560H'
	# adc_res = np.double(sys.argv[11])  # adc resolution of the digitizer
	# sig_offset = float(sys.argv[12])  # percentage offset for ccd signal (value within [-50, 50]) look at wavedump config file
	# trig_offset = float(sys.argv[13])  # percentage offset for trigger signal (value within [-50, 50]) look at wavedump config file
	# anti_co_offset = float(sys.argv[14])  # percentage offset for veto signal (value within [-50, 50]) look at wavedump config file
	# time_res = np.double(sys.argv[15])  # time resolution of digitizer. i.e. 2e-9
	# post_trig_percent = float(sys.argv[16])  # percentage of the acquired data after trigger (look at wave dump config file)
	# trig_value = np.double(sys.argv[17])  # volts below base line for trigger
	# veto_value = int(sys.argv[18])  # counts below base line on the veto signal for vetoing
	# dig_bits = int(sys.argv[19])  # number of bits of the ADC i.e. 14
	# simultaneous_conversion = bool(sys.argv[20] != '0')  # whether or not to do the simultaneous conversion while taking data
	# time_recal = float(sys.argv[21])  # time between digitiser recalibrations
	# control_hv = bool(sys.argv[22] != '0')  # whether or not voltage and currents are controlled
	# polarity = int(sys.argv[23])  # polarity of the bias voltage
	#
	# converter = Converter_Caen(run_dir_location, working_dir_location, filename, signal_path, trigger_path, anti_co_path, points,
	#                            num_events, struct_len, struct_fmt, adc_res, sig_offset, trig_offset, anti_co_offset,
	#                            time_res, post_trig_percent, trig_value, veto_value, dig_bits, simultaneous_conversion, time_recal, control_hv, polarity)
	#
	converter.SetupRootFile(sys.argv)
	converter.GetBinariesNumberWrittenEvents()
	converter.OpenRawBinaries()
	converter.CreateProgressBar(converter.num_events)
	converter.ConvertEvents()
	converter.CloseAll()




