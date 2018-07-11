#!/usr/bin/env python
import visa
import csv
from struct import unpack
import time, os, sys
import progressbar
import ipdb
from pykeyboard import PyKeyboard
from ConfigParser import ConfigParser
import subprocess as subp
from optparse import OptionParser
from Settings_Caen import Settings_Caen
from Channel_Caen import Channel_Caen
import pickle
import struct
import numpy as np
import ROOT as ro
import shutil
from copy import deepcopy

import glob

class Channel_Calibration:
	def __init__(self, run_dir, prefix, suffix, ch):
		self.run_dir = run_dir
		self.prefix = prefix
		self.suffix = suffix
		self.ch = ch

		self.list_settings = glob.glob('{d}/{p}*{s}.settings'.format(d=self.run_dir, p=self.prefix, s=self.suffix))
		self.list_volts_settings = [settings.split('/')[-1].split(self.prefix)[-1].split(self.suffix)[0] for settings in self.list_settings]
		self.list_volts_pos = [float(settings.split('Pos_')[-1].split('V')[0]) for settings in self.list_volts_settings if settings.startswith('Pos')]
		self.list_volts_neg = [-float(settings.split('Neg_')[-1].split('V')[0]) for settings in self.list_volts_settings if settings.startswith('Neg')]
		self.list_volts = self.list_volts_neg + self.list_volts_pos

		# voltage in reality are inverted. That means that negative volts are in reality positive volts, and vice versa. This happened because the CCD has negative signals for positive bias, and the other way around. To fix, everything has to be multiplied by -1 and treated accordingly.
		self.list_volts = [-v for v in self.list_volts]
		self.list_volts.sort()

		self.dic_settings, self.dic_channel, self.dic_adcs = {}, {}, {}
		for volt in self.list_volts:
			if volt < 0:
				self.dic_settings[volt] = pickle.load(open('{d}/{p}Pos_{v}V{s}.settings'.format(d=self.run_dir, p=self.prefix, v=abs(volt), s=self.suffix), 'rb'))
				self.dic_channel[volt] = pickle.load(open('{d}/{p}Pos_{v}V{s}.signal'.format(d=self.run_dir, p=self.prefix, v=abs(volt), s=self.suffix), 'rb'))
				self.dic_adcs[volt] = open('{d}/{p}Pos_{v}V{s}_signal.dat'.format(d=self.run_dir, p=self.prefix, v=abs(volt), s=self.suffix), 'rb')
			else:
				self.dic_settings[volt] = pickle.load(open('{d}/{p}Neg_{v}V{s}.settings'.format(d=self.run_dir, p=self.prefix, v=volt, s=self.suffix), 'rb'))
				self.dic_channel[volt] = pickle.load(open('{d}/{p}Neg_{v}V{s}.signal'.format(d=self.run_dir, p=self.prefix, v=volt, s=self.suffix), 'rb'))
				self.dic_adcs[volt] = open('{d}/{p}Neg_{v}V{s}_signal.dat'.format(d=self.run_dir, p=self.prefix, v=volt, s=self.suffix), 'rb')

		self.dic_adcs_with_offset = {}

		for volt in self.list_volts:
			list_s = []
			for ev in xrange(self.dic_settings[volt].num_events):
				self.dic_adcs[volt].seek(ev * self.dic_settings[volt].struct_len, 0)
				datas = self.dic_adcs[volt].read(self.dic_settings[volt].struct_len)
				list_s += struct.Struct(self.dic_settings[volt].struct_fmt).unpack_from(datas)
			array_s = np.array(list_s)
			mean_s = array_s.mean(dtype='float64')
			std_s = array_s.std(dtype='float64')
			for it in xrange(10):
				condition_s = (np.abs(array_s - mean_s) < 4 * std_s)
				mean_s = np.extract(condition_s, array_s).mean(dtype='float64')
				std_s = np.extract(condition_s, array_s).std(dtype='float64')
			self.dic_adcs_with_offset[volt] = {'mean_wo_offset': mean_s, 'sigma': std_s, 'mean_w_offset': np.add(mean_s, np.multiply(2**self.dic_settings[volt].dig_bits - 1.0, self.dic_channel[volt].dc_offset_percent/100.0 - 0.5, dtype='float64'), dtype='float64')}  # the dc_offset_percent is multiplied by -1 because normally the signal is inverted in the data acquisition, that's why it requires a - here.

		temp_x = np.array([self.dic_adcs_with_offset[volt]['mean_w_offset'] for volt in self.list_volts], 'float64')
		temp_y = np.array(self.list_volts, 'float64')
		temp_ex = np.array([self.dic_adcs_with_offset[volt]['sigma'] for volt in self.list_volts], 'float64')
		temp_ey = np.array([abs(volt) * 0.006 + 0.0002 if abs(volt) <= 0.6 else abs(volt) * 0.003 + 0.002 for volt in self.list_volts], dtype='float64')
		self.graph_all = ro.TGraphErrors(len(self.list_volts), temp_x, temp_y, temp_ex, temp_ey)
		self.graph_all.SetNameTitle('Vcal_vs_OffsetedADCs_ch' + str(self.ch), 'Vcal_vs_OffsetedADCs_ch' + str(self.ch))
		self.graph_all.GetXaxis().SetTitle('Offseted_ADC')
		self.graph_all.GetYaxis().SetTitle('Vcal/V')

		temp_x = np.array([self.dic_adcs_with_offset[volt]['mean_w_offset'] for volt in self.list_volts if volt >= 0], 'float64')
		temp_y = np.array([v for v in self.list_volts if v >= 0], 'float64')
		temp_ex = np.array([self.dic_adcs_with_offset[volt]['sigma'] for volt in self.list_volts if volt >= 0], 'float64')
		temp_ey = np.array([abs(volt) * 0.006 + 0.0002 if abs(volt) <= 0.6 else abs(volt) * 0.003 + 0.002 for volt in self.list_volts if volt >= 0], dtype='float64')
		self.graph_pos = ro.TGraphErrors(temp_x.size, temp_x, temp_y, temp_ex, temp_ey)
		self.graph_pos.SetNameTitle('posVcal_vs_OffsetedADCs_ch' + str(self.ch), 'posVcal_vs_OffsetedADCs_ch' + str(self.ch))
		self.graph_pos.GetXaxis().SetTitle('Offseted_ADC')
		self.graph_pos.GetYaxis().SetTitle('Vcal/V')
		fit_pos_func = ro.TF1('pos_fit', 'pol1', -500, 16500)
		fit_pos_func.SetParameters(1.0e-3, 1.27e-4)
		self.pos_fit = self.graph_pos.Fit('pos_fit', 'QMFS')

		temp_x = np.array([self.dic_adcs_with_offset[volt]['mean_w_offset'] for volt in self.list_volts if volt < 0], 'float64')
		temp_y = np.array([v for v in self.list_volts if v < 0], 'float64')
		temp_ex = np.array([self.dic_adcs_with_offset[volt]['sigma'] for volt in self.list_volts if volt < 0], 'float64')
		temp_ey = np.array([abs(volt) * 0.006 + 0.0002 if abs(volt) <= 0.6 else abs(volt) * 0.003 + 0.002 for volt in self.list_volts if volt < 0], dtype='float64')
		self.graph_neg = ro.TGraphErrors(temp_x.size, temp_x, temp_y, temp_ex, temp_ey)
		self.graph_neg.SetNameTitle('negVcal_vs_OffsetedADCs_ch' + str(self.ch), 'negVcal_vs_OffsetedADCs_ch' + str(self.ch))
		self.graph_neg.GetXaxis().SetTitle('Offseted_ADC')
		self.graph_neg.GetYaxis().SetTitle('Vcal/V')
		fit_neg_func = ro.TF1('neg_fit', 'pol1', -16500, 500)
		fit_neg_func.SetParameters(-1.0e-3, 1.27e-4)
		self.neg_fit = self.graph_neg.Fit('neg_fit', 'QMFS')

		self.file_out = None

	def SaveVoltageCalibration(self):
		self.file_out = ro.TFile('{d}/calibration.root'.format(d=self.run_dir), 'RECREATE')
		self.graph_all.Write()
		self.pos_fit.Write()
		self.graph_pos.Write()
		self.neg_fit.Write()
		self.graph_neg.Write()
		self.file_out.Write()
		self.file_out.Close()


if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-d', '--rundir', dest='rundir', default='.', type='string', help='Directory where all the calibration files are for a channel')
	parser.add_option('-p', '--prefix', dest='prefix', default='', type='string', help='Prefix of the files before the voltage value. e.g. "cms11_volts_ccd_"')
	parser.add_option('-s', '--suffix', dest='suffix', default='', type='string', help='Suffix of the files after the voltage value. e.g. "_t1"')
	parser.add_option('-c', '--channel', dest='channel', default=0, type='int', help='Channel to which the files belong. i.e. the channel to be calibrated')

	(options, args) = parser.parse_args()
	rundir = str(options.rundir)
	prefix = str(options.prefix)
	suffix = str(options.suffix)
	ch = int(options.channel)
	z = Channel_Calibration(rundir, prefix, suffix, ch)