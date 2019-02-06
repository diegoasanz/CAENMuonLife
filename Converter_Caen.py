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
	def __init__(self, settings_object='', data_path=''):
		self.settings_object = settings_object

		self.settings_full_path = os.path.abspath(settings_object)
		self.output_dir = '/'.join(self.settings_full_path.split('/')[:-1])
		self.raw_dir = self.output_dir if data_path == '' else data_path
		self.filename = self.settings_full_path.split('/')[-1].split('.settings')[0]
		self.settings = pickle.load(open('{d}/{f}.settings'.format(d=self.output_dir, f=self.filename), 'rb'))
		self.signal_ch = pickle.load(open('{d}/{f}.signal_ch'.format(d=self.output_dir, f=self.filename), 'rb'))
		self.trigger_ch = pickle.load(open('{d}/{f}.trigger_ch'.format(d=self.output_dir, f=self.filename), 'rb'))
		self.veto_ch = pickle.load(open('{d}/{f}.veto'.format(d=self.output_dir, f=self.filename), 'rb'))

		self.signal_path = data_path + '/raw_wave{chs}.dat'.format(chs=self.settings.sigCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_signal.dat'
		self.trigger_path = data_path + '/raw_wave{cht}.dat'.format(cht=self.settings.trigCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_trigger.dat'
		self.veto_path = data_path + '/raw_wave{cha}.dat'.format(cha=self.settings.acCh) if self.settings.simultaneous_conversion else data_path + '/' + self.filename + '_veto.dat'
		self.points = self.settings.points
		self.num_events = self.settings.num_events
		self.struct_len = self.settings.struct_len
		self.struct_fmt = self.settings.struct_fmt
		self.adc_res = self.settings.sigRes
		self.adc_offset = 0
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
		self.hv_pos = 0
		self.hv_struct_fmt = self.settings.hv_struct_fmt
		self.hv_struct_len = self.settings.hv_struct_len

		self.trigger_search_window = 0.1e-6
		self.veto_window_around_trigg = 50e-9
		self.peak_pos_estimate = 2.131e-6
		self.peak_pos_window = 0.5e-6

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
		self.hv_raw_data = None
		self.hv_struct = None
		self.hv_data = {'event': 0, 'seconds': 0, 'nanoseconds': 0, 'voltage': 0, 'current': 0}
		# self.hour_min_sec_event = None
		# self.currentTime = None

		self.struct_s = self.struct_t = self.struct_ac = None

		self.bar = None

	def SetupRootFile(self):
		if self.simultaneous_conversion:
			print 'Start creating root file simultaneously with data taking'
		else:
			print 'Start creating root file'
		self.raw_file = ro.TFile('{wd}/{r}.root'.format(wd=self.output_dir, r=self.filename), 'RECREATE')
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
			self.hourMinSecBra = ro.TTimeStamp()
		self.raw_tree.Branch('event', self.eventBra, 'event/i')
		self.raw_tree.Branch('time', self.timeBra, 'time[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageSignal', self.voltBra, 'voltageSignal[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageTrigger', self.trigBra, 'voltageTrigger[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageVeto', self.vetoBra, 'voltageVeto[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('vetoedEvent', self.vetoedBra, 'vetoedEvent/O')
		self.raw_tree.Branch('badShape', self.badShapeBra, 'badShape/B')  # signed char
		self.raw_tree.Branch('badPedestal', self.badPedBra, 'badPedestal/O')
		if self.control_hv:
			self.raw_tree.Branch('voltageHV', self.hvVoltageBra, 'voltageHV/F')
			self.raw_tree.Branch('currentHV', self.hvCurrentBra, 'currentHV/F')
			self.raw_tree.Branch('timeHV', self.hourMinSecBra)

	def GetBinariesNumberWrittenEvents(self):
		self.signal_written_events = int(round(os.path.getsize(self.signal_path) / self.struct_len)) if os.path.isfile(self.signal_path) else 0
		self.trigger_written_events = int(round(os.path.getsize(self.trigger_path) / self.struct_len)) if os.path.isfile(self.trigger_path) else 0
		self.anti_co_written_events = int(round(os.path.getsize(self.veto_path) / self.struct_len)) if os.path.isfile(self.veto_path) else 0

	def OpenRawBinaries(self):
		self.fs = open(self.signal_path, 'rb')
		self.ft = open(self.trigger_path, 'rb')
		self.fa = open(self.veto_path, 'rb')

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

	def ConvertEvents(self):
		self.bar.start()
		for ev in xrange(self.num_events):
			self.CheckFilesSizes(ev)
			self.WaitForData(ev)
			self.ReadData(ev)
			self.CheckData()
			self.struct_s = struct.Struct(self.struct_fmt).unpack_from(self.datas)
			self.sigADC = np.array(self.struct_s, 'H')
			self.sigVolts = self.ADC_to_Volts('signal')
			self.struct_t = struct.Struct(self.struct_fmt).unpack_from(self.datat)
			self.trigADC = np.array(self.struct_t, 'H')
			self.trigVolts = self.ADC_to_Volts('trigger')
			self.LookForTime0()
			self.timeVect = np.linspace(-self.trigPos * self.time_res, self.time_res * (self.points - 1 - self.trigPos), self.points, dtype='f8')
			self.struct_ac = struct.Struct(self.struct_fmt).unpack_from(self.dataa)
			self.vetoADC = np.array(self.struct_ac, 'H')
			self.vetoVolts = self.ADC_to_Volts('veto')
			self.vetoed_event = self.IsEventVetoed()
			self.DefineSignalBaseLineAndPeakPosition()
			self.bad_shape_event = self.IsEventBadShape()
			self.bad_pedstal_event = self.IsPedestalBad()
			self.FillBranches(ev)
			if ev == 10:
				self.raw_tree.OptimizeBaskets()
			self.raw_tree.Fill()
			self.bar.update(ev + 1)

	def CheckFilesSizes(self, ev):
		self.wait_for_data = (self.signal_written_events <= ev) or (self.trigger_written_events <= ev) or (self.anti_co_written_events <= ev)

	def WaitForData(self, ev):
		t1 = time.time()
		while self.wait_for_data:
			if self.simultaneous_conversion:
				if time.time() - t1 > self.time_break:
					print 'No data has been saved in file for event {ev} in the past {t} seconds... exiting!'.format(ev=ev, t=self.time_break)
					exit(os.EX_NOINPUT)
				if not self.fs.closed:
					self.fs.close()
				if not self.ft.closed:
					self.ft.close()
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
		self.fa.seek(ev * self.struct_len, 0)
		self.dataa = self.fa.read(self.struct_len)
		if self.control_hv:
			self.Read_HV_File(ev)

	def Read_HV_File(self, ev):
		line = [0, 0]
		temp_line = ''
		if os.path.isfile(self.hv_file_name):
			hv_elems = int(round(os.path.getsize(self.hv_file_name) / self.hv_struct_len))
			with open(self.hv_file_name, 'rb') as self.file_hv:
				for pi in xrange(self.hv_pos, hv_elems):
					self.file_hv.seek(pi * self.hv_struct_len, 0)
					self.hv_raw_data = self.file_hv.read(self.hv_struct_len)
					self.hv_struct = struct.Struct(self.hv_struct_fmt).unpack_from(self.hv_raw_data)
					if self.hv_struct[0] <= ev:
						self.hv_pos = pi
					else:
						self.hv_pos = pi - 1
					if self.hv_struct[0] >= ev:
						break
				self.file_hv.seek(self.hv_pos * self.hv_struct_len, 0)
				self.hv_raw_data = self.file_hv.read(self.hv_struct_len)
				self.hv_struct = struct.Struct(self.hv_struct_fmt).unpack_from(self.hv_raw_data)
				self.hv_data['event'], self.hv_data['seconds'], self.hv_data['nanoseconds'], self.hv_data['voltage'], self.hv_data['current'] = ev, self.hv_struct[1], self.hv_struct[2], self.hv_struct[3], self.hv_struct[4]

		# if self.simultaneous_conversion:
		# 	if os.path.isfile(self.hv_file_name):
		# 		with open(self.hv_file_name, 'rb') as self.file_hv:
		# 			self.file_hv.seek(self.hv_pos * self.hv_struct_len, 0)
		# 			self.hv_raw_data = self.file_hv.read(self.hv_struct_len)
		# 			self.hv_struct = struct.Struct(self.hv_struct_fmt).unpack_from(self.hv_raw_data)
		# 			if self.hv_struct[0] > ev:
		# 				self.hv_pos
		# 			temp_line = self.file_hv.readline()
		# 			line = temp_line.split() if temp_line != '' else line
		# self.hv_voltage_event, self.hv_current_event = float(line[0]), float(line[1])
		del line, temp_line, self.file_hv
		self.file_hv = None

	def CheckData(self):
		if not self.datas or not self.datat or not self.dataa:
			print 'No event in signal or trigger files... exiting'
			exit(os.EX_DATAERR)

	def LookForTime0(self):
		guess_pos = int(round(self.points * (100.0 - self.post_trig_percent)/100.0))
		condition_trigg = np.array(np.abs(self.array_points - guess_pos) <= int(round(self.trigger_search_window/self.time_res)), dtype='?')
		condition_no_trigg = np.array(1 - condition_trigg, dtype='?')
		# mean = np.extract(condition_no_trigg, self.trigVolts).mean()
		# sigma = np.extract(condition_no_trigg, self.trigVolts).std()
		temp_trig_volts = np.copy(self.trigVolts)
		np.putmask(temp_trig_volts, condition_no_trigg, 100)
		volt_min_pos = temp_trig_volts.argmin()
		condition_trigg = np.bitwise_and(condition_trigg, np.array(self.array_points <= volt_min_pos))
		np.putmask(temp_trig_volts, np.bitwise_not(condition_trigg), 100)
		self.trigPos = np.abs(temp_trig_volts - self.trig_value).argmin()
		del guess_pos, condition_trigg, condition_no_trigg, temp_trig_volts, volt_min_pos

	def IsEventVetoed(self):
		condition_veto_base_line = np.array(np.abs(self.array_points - self.trigPos) > int(round(self.veto_window_around_trigg / float(self.time_res))), dtype='?')
		condition_search = np.array(1 - condition_veto_base_line, dtype='?')
		meanbl = np.extract(condition_veto_base_line, self.vetoADC).mean()
		# sigma = np.extract(condition_veto_base_line, self.vetoADC).std()
		# vetoValNew = 4 * sigma if self.veto_value < 4 * sigma else self.veto_value
		# veto_event = bool((np.extract(condition_search, self.vetoADC) - mean + vetoValNew).min() <= 0)
		veto_event = bool((np.extract(condition_search, self.vetoADC) - meanbl + self.veto_value).min() <= 0)
		del condition_search, condition_veto_base_line, meanbl
		return veto_event

	def DefineSignalBaseLineAndPeakPosition(self):
		self.condition_base_line = np.array(self.array_points <= self.trigPos, dtype='?')
		# values 2.2 and 0.5 us come from shape of signal
		self.condition_peak_pos = np.array(np.abs(self.array_points - (self.peak_pos_estimate/float(self.time_res) + self.trigPos)) <= self.peak_pos_window/float(self.time_res), dtype='?')

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

	def IsPedestalBad(self):
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		mean = np.extract(self.condition_base_line, self.sigADC).mean()
		self.adc_res = self.signal_ch.offseted_adc_to_volts_cal['p1']
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

	def FillBranches(self, ev):
		self.eventBra.fill(ev)
		np.putmask(self.timeBra, np.bitwise_not(np.zeros(self.points, '?')), self.timeVect)
		np.putmask(self.voltBra, np.bitwise_not(np.zeros(self.points, '?')), self.sigVolts)
		np.putmask(self.trigBra, np.bitwise_not(np.zeros(self.points, '?')), self.trigVolts)
		np.putmask(self.vetoBra, np.bitwise_not(np.zeros(self.points, '?')), self.vetoVolts)
		self.vetoedBra.fill(self.vetoed_event)
		self.badShapeBra.fill(self.bad_shape_event)
		self.badPedBra.fill(self.bad_pedstal_event)
		if self.control_hv:
			self.hvVoltageBra.fill(self.hv_data['voltage'])
			self.hvCurrentBra.fill(self.hv_data['current'])
			self.hourMinSecBra.Set(1970, 1, 1, 0, 0, self.hv_data['seconds'], self.hv_data['nanoseconds'], True, 0)

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

	def IsPedestalNotFlat(self, signalADC, points, trigPos, time_res):
		array_points = np.arange(points, dtype=np.dtype('int32'))
		condition_base_line = np.array(array_points - trigPos <= 0, dtype='?')

	def ADC_to_Volts(self, sig_type):
		adcs, offset = 0, 0
		if sig_type == 'signal':
			adcs = self.sigADC
			offset = self.sig_offset
			self.adc_offset = self.signal_ch.offseted_adc_to_volts_cal['p0']
			self.adc_res = self.signal_ch.offseted_adc_to_volts_cal['p1']
		elif sig_type == 'trigger':
			adcs = self.trigADC
			offset = self.trig_offset
			self.adc_offset = self.trigger_ch.offseted_adc_to_volts_cal['p0']
			self.adc_res = self.trigger_ch.offseted_adc_to_volts_cal['p1']
		elif sig_type == 'veto':
			adcs = self.vetoADC
			offset = self.anti_co_offset
			self.adc_offset = self.veto_ch.offseted_adc_to_volts_cal['p0']
			self.adc_res = self.veto_ch.offseted_adc_to_volts_cal['p1']
		else:
			print 'Wrong type. Exiting'
			exit()
		result = np.add(self.adc_offset, np.multiply(self.adc_res, np.add(adcs, np.multiply(2 ** self.dig_bits - 1.0, offset / 100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8'), dtype='f8')
		return result

if __name__ == '__main__':
	settings_object = str(sys.argv[1])  # settings pickle path
	if len(sys.argv) > 2:
		data_path = str(sys.argv[2])  # path where the binary data in adcs is
	else:
		data_path = ''
	converter = Converter_Caen(settings_object=settings_object, data_path=data_path)

	converter.SetupRootFile()
	converter.GetBinariesNumberWrittenEvents()
	converter.OpenRawBinaries()
	converter.CreateProgressBar(converter.num_events)
	converter.ConvertEvents()
	converter.CloseAll()




