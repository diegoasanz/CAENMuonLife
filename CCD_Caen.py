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
# from memory_profiler import profile

trig_rand_time = 0.2
wait_time_hv = 7

class CCD_Caen:
	def __init__(self, infile='None', verbose=False, settingsObj=None):
		print 'Starting CCD program ...'
		self.infile = infile
		self.verb = verbose
		if self.infile != 'None':
			self.settings = Settings_Caen(self.infile, self.verb)
			self.settings.ReadInputFile()
		elif settingsObj:
			self.settings = settingsObj
		else:
			ExitMessage('No setting file was given, or settings object. Quittint!')
		self.settings.Get_Calibration_Constants()
		self.settings.SetOutputFiles()

		# Create channel objects for signal, trigger and veto
		self.signal_ch = Channel_Caen(self.settings.sigCh, 'signal_ch', self.verb)
		self.signal_ch.Set_Channel(self.settings)
		self.trigger_ch = Channel_Caen(self.settings.trigCh, 'trigger_ch', self.verb)
		self.trigger_ch.Set_Channel(self.settings)
		self.veto_ch = Channel_Caen(self.settings.vetoCh, 'veto', self.verb)
		self.veto_ch.Set_Channel(self.settings)

		# declare extra variables that will be used
		self.fs0, self.ft0, self.fv0 = None, None, None
		self.hv_control = None
		self.utils = Utils()
		self.RemoveFiles()
		self.t0, self.t1, self.t2 = None, None, None
		self.p, self.pconv = None, None
		self.total_events = 0
		self.written_events_sig, self.written_events_trig, self.written_events_veto = 0, 0, 0
		self.total_events_sig, self.total_events_trig, self.total_events_veto = 0, 0, 0
		self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_veto = 0, 0, 0
		self.total_merged_data_sig, self.total_merged_data_trig, self.total_merged_data_veto = 0, 0, 0
		self.doMerge = False
		self.min_measured_data = 0
		self.min_data_to_write = 0
		self.events_to_write = 0
		self.read_size = 0
		self.sig_written, self.trg_written, self.veto_written = 0, 0, 0
		self.session_written_events_sig, self.session_written_events_trg, self.session_written_events_veto = 0, 0, 0
		self.fins, self.fint, self.finv = None, None, None
		self.datas, self.datat, self.datav = None, None, None

	def RemoveFiles(self):
		# used, for example, to remove old files that may have stayed due to crashes
		if self.fs0:
			if not self.fs0.closed:
				self.fs0.close()
				del self.fs0
		if self.ft0:
			if not self.ft0.closed:
				self.ft0.close()
				del self.ft0
		if self.fv0:
			if not self.fv0.closed:
				self.fv0.close()
				del self.fv0
		channels = [self.signal_ch.ch, self.trigger_ch.ch, self.veto_ch.ch]
		for ch in channels:
			if os.path.isfile('raw_waves{c}.dat'.format(c=ch)):
				os.remove('raw_waves{c}.dat'.format(c=ch))
			if os.path.isfile('waves{c}.dat'.format(c=ch)):
				os.remove('waves{c}.dat'.format(c=ch))
		del channels

	def StartHVControl(self):
		if self.settings.do_hv_control:
			self.hv_control = HV_Control(self.settings)
			print 'Waiting {t} seconds for the HVClient to start... '.format(t=wait_time_hv), ; sys.stdout.flush()
			time.sleep(wait_time_hv)
			print 'Done'
			self.hv_control.CheckVoltage()

	def GetBaseLines(self):
		self.settings.SetupDigitiser(doBaseLines=True, signal=self.signal_ch, trigger=self.trigger_ch, ac=self.veto_ch)
		self.p = subp.Popen(['{p}/wavedump'.format(p=self.settings.wavedump_path), '{d}/WaveDumpConfig_CCD_BL.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		t0 = time.time()
		self.CreateEmptyFiles()
		self.CloseFiles()
		self.GetWaveforms(events=1, stdin=True, stdout=False)
		if self.total_events >= 1:
			self.ReadBaseLines()
			self.settings.RemoveBinaries()
			self.RemoveFiles()
		print 'Time getting base lines: {t} seconds'.format(t=(time.time() - t0))
		del t0

	def CreateEmptyFiles(self):
		self.ft0 = open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), 'wb')
		self.fs0 = open('raw_wave{s}.dat'.format(s=self.signal_ch.ch), 'wb')
		self.fv0 = open('raw_wave{a}.dat'.format(a=self.veto_ch.ch), 'wb')

	def OpenFiles(self, mode='rb'):
		if not self.fs0:
			self.fs0 = open('raw_wave{s}.dat'.format(s=self.signal_ch.ch), mode)
		if not self.ft0:
			self.ft0 = open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), mode)
		if not self.fv0:
			self.fv0 = open('raw_wave{a}.dat'.format(a=self.veto_ch.ch), mode)

	def CloseFiles(self):
		if self.ft0:
			self.ft0.close()
			if self.ft0.closed:
				del self.ft0
				self.ft0 = None
		if self.fs0:
			self.fs0.close()
			if self.fs0.closed:
				del self.fs0
				self.fs0 = None
		if self.fv0:
			self.fv0.close()
			if self.fv0.closed:
				del self.fv0
				self.fv0 = None

	def GetWaveforms(self, events=1, stdin=False, stdout=False):
		self.t1 = time.time()
		# if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
		if events == 1:
			# while self.p.poll() is None:
			time.sleep(1)
			self.p.stdin.write('c')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			if self.settings.plot_waveforms:
				# time.sleep(1)
				self.p.stdin.write('P')
				self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('w')
			self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('t')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('q')
			self.p.stdin.flush()
			while self.p.poll() is None:
				continue
			if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
			self.ConcatenateBinaries()
			self.CloseSubprocess('wave_dump', stdin=stdin, stdout=stdout)
			self.settings.RemoveBinaries()
		else:
			time.sleep(1)
			self.p.stdin.write('c')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('W')
			self.p.stdin.flush()
			if self.settings.plot_waveforms:
				# time.sleep(1)
				self.p.stdin.write('P')
				self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			self.written_events_sig, self.written_events_trig, self.written_events_veto = 0, 0, 0
			# time.sleep(1)
			self.t2 = time.time()
			while self.p.poll() is None:
				if time.time() - self.t1 >= self.settings.time_calib:
					self.p.stdin.write('s')
					self.p.stdin.flush()
					self.settings.RemoveBinaries()
					self.p.stdin.write('c')
					self.p.stdin.flush()
					self.p.stdin.write('q')
					self.p.stdin.flush()
					time.sleep(1)
					# if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
				elif self.written_events_sig + self.sig_written >= events:
					self.p.stdin.write('s')
					self.p.stdin.flush()
					self.settings.RemoveBinaries()
					self.p.stdin.write('q')
					self.p.stdin.flush()
					time.sleep(1)
					# if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
				else:
					if self.settings.random_test and (time.time() - self.t2 > trig_rand_time):
						self.p.stdin.write('t')
						self.p.stdin.flush()
						self.t2 = time.time()
					self.ConcatenateBinaries()
					if self.settings.do_hv_control:
						self.hv_control.UpdateHVFile(int(min(self.written_events_sig + self.sig_written, self.settings.num_events)))
					if not self.settings.simultaneous_conversion:
						self.settings.bar.update(int(min(self.written_events_sig + self.sig_written, self.settings.num_events)))
			del self.t1
			self.t1 = None
			self.CloseSubprocess('wave_dump', stdin=stdin, stdout=stdout)
			time.sleep(1)
			del self.t2
			self.t2 = None
			if self.settings.do_hv_control: self.hv_control.UpdateHVFile(int(min(self.written_events_sig + self.sig_written, self.settings.num_events)))
		self.total_events_sig = self.CalculateEventsWritten(self.signal_ch.ch)
		self.total_events_trig = self.CalculateEventsWritten(self.trigger_ch.ch)
		self.total_events_veto = self.CalculateEventsWritten(self.veto_ch.ch)
		if self.total_events_sig == self.total_events_trig and self.total_events_sig == self.total_events_veto:
			self.total_events = self.total_events_sig
		else:
			print 'Written events are of different sizes (signal: {s}, trigger: {t}, veto: {v}). Missmatch!'.format(s=self.total_events_sig, t=self.total_events_trig, v=self.total_events_veto)
			exit()
		del self.total_events_sig, self.total_events_trig, self.total_events_veto
		self.total_events_sig, self.total_events_trig, self.total_events_veto = None, None, None

	def CloseSubprocess(self, pname='wave_dump', stdin=False, stdout=False):
		p = self.p if pname == 'wave_dump' else self.pconv if pname == 'converter' else None
		if not p:
			print 'Something failed! Exiting!'
			exit()
		pid = p.pid
		if stdin:
			p.stdin.close()
		if stdout:
			p.stdout.close()
		if p.wait() is None:
			print 'Could not terminate subprocess... forcing termination'
			p.kill()
			if p.wait() is None:
				print 'Could not kill subprocess... quitting'
				exit()
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			print 'The subprocess is still running. Killing it with os.kill'
			os.kill(pid, 15)
			try:
				os.kill(pid, 0)
			except OSError:
				pass
			else:
				print 'The process does not die... quitting program'
				exit()
		del p, pid

		if pname == 'wave_dump':
			del self.p
			self.p = None
		elif pname == 'converter':
			del self.pconv
			self.pconv = None

	def ConcatenateBinaries(self):
		self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_veto = 0, 0, 0
		if os.path.isfile('wave{s}.dat'.format(s=self.signal_ch.ch)) and os.path.isfile('wave{t}.dat'.format(t=self.trigger_ch.ch)) and os.path.isfile('wave{a}.dat'.format(a=self.veto_ch.ch)):
			self.session_measured_data_sig = int(os.path.getsize('wave{s}.dat'.format(s=self.signal_ch.ch)))
			self.session_measured_data_trig = int(os.path.getsize('wave{t}.dat'.format(t=self.trigger_ch.ch)))
			self.session_measured_data_veto = int(os.path.getsize('wave{a}.dat'.format(a=self.veto_ch.ch)))

		self.total_merged_data_sig = int(os.path.getsize('raw_wave{s}.dat'.format(s=self.signal_ch.ch)))
		self.total_merged_data_trig = int(os.path.getsize('raw_wave{t}.dat'.format(t=self.trigger_ch.ch)))
		self.total_merged_data_veto = int(os.path.getsize('raw_wave{a}.dat'.format(a=self.veto_ch.ch)))
		self.doMerge = (self.session_measured_data_sig + self.sig_written * self.settings.struct_len > self.total_merged_data_sig) and (self.session_measured_data_trig + self.trg_written * self.settings.struct_len > self.total_merged_data_trig) and (self.session_measured_data_veto + self.veto_written * self.settings.struct_len > self.total_merged_data_veto)
		if self.doMerge:
			# self.OpenFiles(mode='ab')
			self.min_measured_data = min(self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_veto)
			data_to_write_sig = self.min_measured_data - self.total_merged_data_sig + self.sig_written * self.settings.struct_len
			data_to_write_trg = self.min_measured_data - self.total_merged_data_trig + self.trg_written * self.settings.struct_len
			data_to_write_aco = self.min_measured_data - self.total_merged_data_veto + self.veto_written * self.settings.struct_len
			self.min_data_to_write = min(data_to_write_sig, data_to_write_trg, data_to_write_aco)
			del data_to_write_sig, data_to_write_trg, data_to_write_aco
			self.events_to_write = int(np.floor(self.min_data_to_write / float(self.settings.struct_len)))
			self.read_size = self.events_to_write * self.settings.struct_len

			with open('wave{s}.dat'.format(s=self.signal_ch.ch), 'rb') as self.fins:
				self.fins.seek(self.written_events_sig * self.settings.struct_len, 0)
				self.datas = self.fins.read(self.read_size)
			del self.fins
			self.fins = None
			with open('raw_wave{s}.dat'.format(s=self.signal_ch.ch), 'ab') as self.fs0:
				self.fs0.write(self.datas)
				self.fs0.flush()
			del self.fs0, self.datas
			self.fs0, self.datas = None, None

			with open('wave{t}.dat'.format(t=self.trigger_ch.ch), 'rb') as self.fint:
				self.fint.seek(self.written_events_trig * self.settings.struct_len, 0)
				self.datat = self.fint.read(self.read_size)
			del self.fint
			self.fint = None
			with open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), 'ab') as self.ft0:
				self.ft0.write(self.datat)
				self.ft0.flush()
			del self.ft0, self.datat
			self.ft0, self.datat = None, None

			with open('wave{a}.dat'.format(a=self.veto_ch.ch), 'rb') as self.finv:
				self.finv.seek(self.written_events_veto * self.settings.struct_len, 0)
				self.datav = self.finv.read(self.read_size)
			del self.finv
			self.finv = None
			with open('raw_wave{a}.dat'.format(a=self.veto_ch.ch), 'ab') as self.fv0:
				self.fv0.write(self.datav)
				self.fv0.flush()
			del self.fv0, self.datav
			self.fv0, self.datav = None, None

			self.written_events_sig += int(self.events_to_write)
			self.written_events_trig += int(self.events_to_write)
			self.written_events_veto += int(self.events_to_write)

			del self.events_to_write, self.read_size
			self.events_to_write, self.read_size = None, None

		self.doMerge = False

	def CalculateEventsWritten(self, ch):
		return int(round(float(os.path.getsize('raw_wave{c}.dat'.format(c=ch))) / float(self.settings.struct_len)))

	def ReadBaseLines(self):
		# Read ADCs from files for trigger and veto scintillators
		with open('raw_wave{t}.dat'.format(t=self.trigger_ch.ch), 'rb') as self.ft0:
			self.ft0.seek(0)
			self.datat = self.ft0.read(self.settings.struct_len)
			t = struct.Struct(self.settings.struct_fmt).unpack_from(self.datat)
			triggADCs = np.array(t, 'H')
			mean_t = triggADCs.mean()
			std_t = triggADCs.std()

		with open('raw_wave{ac}.dat'.format(ac=self.veto_ch.ch), 'rb') as self.fv0:
			self.fv0.seek(0)
			self.datav = self.fv0.read(self.settings.struct_len)
			ac = struct.Struct(self.settings.struct_fmt).unpack_from(self.datav)
			acADCs = np.array(ac, 'H')
			mean_ac = acADCs.mean()
			std_ac = acADCs.std()

		# clear possible ADCs with non-baseline signals
		for i in xrange(10):
			condition_t = (np.abs(triggADCs - mean_t) < 3 * std_t)
			mean_t = np.extract(condition_t, triggADCs).mean()
			std_t = np.extract(condition_t, triggADCs).std()
			condition_ac = (np.abs(acADCs - mean_ac) < 3 * std_ac)
			mean_ac = np.extract(condition_ac, acADCs).mean()
			std_ac = np.extract(condition_ac, acADCs).std()

		# set channels such that the baselines are near the maximum ADC's leaving space for the scintillator signals. Adjust threshold values
		self.trigger_ch.Correct_Base_Line(mean_adc=mean_t, sigma_adc=std_t, settings=self.settings)
		self.trigger_ch.Correct_Threshold(sigma=std_t)
		self.settings.trig_base_line = np.multiply(self.trigger_ch.base_line_u_adcs, self.settings.sigRes, dtype='f8')
		self.settings.trig_thr_counts = self.trigger_ch.thr_counts
		self.veto_ch.Correct_Base_Line(mean_adc=mean_ac, sigma_adc=std_ac, settings=self.settings)
		self.veto_ch.Correct_Threshold(sigma=std_ac)
		self.settings.ac_base_line = np.multiply(self.veto_ch.base_line_u_adcs, self.settings.sigRes, dtype='f8')
		self.settings.ac_thr_counts = self.veto_ch.thr_counts

		del self.ft0, self.datat, t, triggADCs, mean_t, std_t
		del self.fv0, self.datav, ac, acADCs, mean_ac, std_ac
		self.ft0, self.datat= None, None
		self.fv0, self.datav = None, None

	# @profile(precision=12)
	def GetData(self):
		self.t0 = time.time()
		self.CreateEmptyFiles()
		self.CloseFiles()
		self.total_events = 0
		print 'Getting {n} events...'.format(n=self.settings.num_events)
		if self.settings.simultaneous_conversion:
			self.CreateRootFile(files_moved=False)
		else:
			self.settings.CreateProgressBar(self.settings.num_events)
			self.settings.bar.start()
		self.settings.SetupDigitiser(doBaseLines=False, signal=self.signal_ch, trigger=self.trigger_ch, ac=self.veto_ch, events_written=self.total_events)
		while self.total_events < self.settings.num_events:
			self.sig_written = self.CalculateEventsWritten(self.signal_ch.ch)
			self.trg_written = self.CalculateEventsWritten(self.trigger_ch.ch)
			self.veto_written = self.CalculateEventsWritten(self.veto_ch.ch)
			self.p = subp.Popen(['{p}/wavedump'.format(p=self.settings.wavedump_path), '{d}/WaveDumpConfig_CCD.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
			self.GetWaveforms(self.settings.num_events, stdin=True, stdout=True)
		self.CloseFiles()
		if not self.settings.simultaneous_conversion:
			print 'Time getting {n} events: {t} seconds'.format(n=self.total_events, t=time.time() - self.t0)
			self.settings.bar.finish()
		else:
			while self.pconv.poll() is None:
				continue
			self.CloseSubprocess('converter', stdin=False, stdout=False)
		return self.total_events

	def CreateRootFile(self, files_moved=False):
		settings_bin_path = os.path.abspath(self.settings.outdir + '/Runs/{f}/{f}.settings'.format(f=self.settings.filename))
		data_bin_path = os.path.abspath(self.settings.outdir + '/Runs/{f}'.format(f=self.settings.filename)) if files_moved else os.getcwd()
		self.pconv = subp.Popen(['python', 'Converter_Caen.py', settings_bin_path, data_bin_path], close_fds=True)
		del settings_bin_path

	def CloseHVClient(self):
		if self.settings.do_hv_control:
			self.hv_control.CloseClient()

	def SavePickles(self):
		# save objects of settings, signal_ch, trigger_ch and veto_ch
		with open('{d}/Runs/{f}/{f}.settings'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fs:
			pickle.dump(self.settings, fs, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}/{f}.signal_ch'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fsig:
			pickle.dump(self.signal_ch, fsig, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}/{f}.trigger_ch'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as ft:
			pickle.dump(self.trigger_ch, ft, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}/{f}.veto'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fv:
			pickle.dump(self.veto_ch, fv, pickle.HIGHEST_PROTOCOL)

	def PrintPlotLimits(self, ti=-5.12e-7, tf=4.606e-6, vmin=-0.7, vmax=0.05):
		print np.double([(tf-ti)/float(self.settings.time_res) +1, ti-self.settings.time_res/2.0,
		                      tf+self.settings.time_res/2.0, (vmax-vmin)/self.settings.sigRes, vmin, vmax])

def main():
	parser = OptionParser()
	parser.add_option('-i', '--infile', dest='infile', default='None', type='string',
	                  help='Input configuration file. e.g. CAENCalibration.cfg')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic conversion and analysis afterwards', action='store_true')

	(options, args) = parser.parse_args()
	infile = str(options.infile)
	auto = bool(options.auto)
	verb = bool(options.verb)
	ccd = CCD_Caen(infile, verb)
	if auto:
		ccd.StartHVControl()
		ccd.GetBaseLines()
		ccd.SavePickles()
		written_events = ccd.GetData()
		ccd.settings.num_events = written_events
		ccd.SavePickles()  # update pickles with the real amount of written events
		ccd.settings.MoveBinaryFiles()
		ccd.settings.RenameDigitiserSettings()
		ccd.CloseHVClient()
		if not ccd.settings.simultaneous_conversion:
			ccd.CreateRootFile(files_moved=True)
			while ccd.pconv.poll() is None:
				time.sleep(3)
			ccd.CloseSubprocess('converter', stdin=False, stdout=False)

	print 'Finished :)'
	sys.stdout.write('\a\a\a')
	sys.stdout.flush()

if __name__ == '__main__':
	main()
