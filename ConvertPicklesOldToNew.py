#!/usr/bin/env python
import imp
import cPickle as pickle
from optparse import OptionParser
import shutil, os
import ipdb

class ConvertPicklesOldToNew:
	def __init__(self, old_module_name, old_module_path, new_module_name, new_module_path, arguments, pickle_file):
		self.pickle_file = pickle_file
		ModuleOld = imp.load_source(old_module_name, old_module_path)
		ModuleNew = imp.load_source(new_module_name, new_module_path)

		# ipdb.set_trace()

		self.new_pickle = eval("ModuleNew.{n}{a}".format(n=new_module_name, a=arguments))
		self.old_pickle = pickle.load(open(pickle_file, 'rb'))

		self.new_attributes = dir(self.new_pickle)
		self.old_attributes = dir(self.old_pickle)

		self.new_variables = [attrib for attrib in self.new_attributes if not attrib.startswith('__') and attrib[0].islower()]

		for var in self.new_variables:
			if var in self.old_attributes:
				exec('self.new_pickle.{v} = self.old_pickle.{v}'.format(v=var))

		user_input = raw_input('type "y" if you want to backup old pickle and save the new one: ')
		if user_input.lower() == 'y':
			print 'Backing up old pickle and saving new one...'
			self.ReplacePickleAndBackUp()
		else:
			print 'Not bakcing up old pickle and not saving the new one. You can modify variables in the pickle and later do this by running on interactive mode'


	def BackUpAndReplacePickle(self):
		shutil.move(pickle_file, pickle_file + '.bkp')
		with open(self.pickle_file, 'wb') as f0:
			pickle.dump(self.new_pickle, f0, pickle.HIGHEST_PROTOCOL)

	def PrintVariablesStatus(self):
		for var in self.new_variables:
			print var, eval('self.new_pickle.{v}'.format(v=var))



if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-o', '--old_module_name', dest='old_module_name', type='string')
	parser.add_option('-O', '--old_module_path', dest='old_module_path', type='string')
	parser.add_option('-n', '--new_module_name', dest='new_module_name', type='string')
	parser.add_option('-N', '--new_module_path', dest='new_module_path', type='string')
	parser.add_option('-a', '--arguments', dest='arguments', type='string', default='()', help='if you pass any do it as: \(3, \ "signal\ "\) without spaces')
	parser.add_option('-p', '--pickle', dest='pickle', type='string')

	(options, args) = parser.parse_args()
	old_module_name = options.old_module_name
	old_module_path = options.old_module_path
	new_module_name = options.new_module_name
	new_module_path = options.new_module_path
	arguments = options.arguments
	pickle_file = options.pickle
	z = ConvertPicklesOldToNew(old_module_name=old_module_name, old_module_path=old_module_path, new_module_name=new_module_name, new_module_path=new_module_path, arguments=arguments, pickle_file=pickle_file)