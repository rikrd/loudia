#! /usr/bin/env python
# encoding: UTF-8
# Thomas Nagy 2008

from fnmatch import fnmatchcase
import os, os.path, re, stat
import pproc
import Task, Utils, Node, Constants
from TaskGen import feature

DOXY_STR = 'cd %s && ${DOXYGEN} ${DOXYFLAGS} - >/dev/null'
DOXY_FMTS = 'html latex man rft xml'.split()
DOXY_EXTS = '''
*.c *.cc *.cxx *.cpp *.c++ *.C
*.h *.hh *.hxx *.hpp *.h++ *.H
*.py *.java *.cs
*.ii *.ixx *.ipp *.i++ *.inl
*.idl *.odl *.php *.php3 *.inc *.m *.mm
'''.split()

re_join = re.compile(r'\\(\r)*\n', re.M)
re_nl = re.compile('\r*\n', re.M)

def read_into_dict(name):
	f = open(name, 'rb')
	txt = f.read()
	f.close()

	ret = {}

	txt = re_join.sub('', txt)
	lines = re_nl.split(txt)
	vals = []

	for x in lines:
		x.strip()
		if len(x) < 2: continue
		if x[0] == '#': continue
		tmp = x.split('=')
		if len(tmp) < 2: continue
		ret[tmp[0].strip()] = '='.join(tmp[1:]).strip()
	return ret

def populate(node, branch, env):
	path = node.abspath(env) + os.sep + branch
	st = os.stat(path)
	if stat.S_ISREG(st[stat.ST_MODE]):
		return [node.exclusive_build_node(branch)]
	elif stat.S_ISDIR(st[stat.ST_MODE]):
		parent = node.ensure_dir_node_from_path(branch)
		lst = Utils.listdir(parent.abspath(env))
		buf = [parent]
		for x in lst:
			buf += populate(parent, x, env)
		return buf

def nodes_files_of(node, recurse, includes, excludes):
	"""
	Recursively returns subnodes of this node. All node filenames that
	matches a pattern in includes but does not match a pattern in
	excludes is returned.

	@param includes: list of wildcards to include.
	@param excludes: list of wildcards to exclude.
	@param recurse: whether to recurse or not.
	"""
	# perform the listdir in the source directory, once
	node.__class__.bld.rescan(node)

	# doxygen looks at the files under the source directory
	buf = []
	for x in node.__class__.bld.cache_dir_contents[node.id]:
		filename = node.abspath() + os.sep + x
		if any(fnmatchcase(filename, pat) for pat in excludes):
			continue
		st = os.stat(filename)
		if stat.S_ISREG(st[stat.ST_MODE]):
			k = node.find_resource(x)
			if not any(fnmatchcase(k.abspath(), pat) for pat in includes):
				continue
			buf.append(k)
		elif stat.S_ISDIR(st[stat.ST_MODE]) and recurse:
			nd = node.find_dir(x)
			if nd.id != nd.__class__.bld.bldnode.id:
				buf += nodes_files_of(nd, recurse, includes, excludes)
	return buf


class doxygen_task(Task.Task):

	vars=['DOXYGEN', 'DOXYFLAGS']
	color = 'BLUE'
	after = 'cxx_link cc_link'
	quiet = True

	def runnable_status(self):
		if not getattr(self, 'pars', None):
			infile = self.inputs[0].abspath(self.env)
			self.pars = read_into_dict(infile)
			if not self.pars.get('OUTPUT_DIRECTORY'):
				self.pars['OUTPUT_DIRECTORY'] = self.inputs[0].parent.abspath(self.env)
		self.signature()
		return Task.Task.runnable_status(self)

	def scan(self):

		recurse = self.pars.get('RECURSIVE') == 'YES'
		excludes = self.pars.get('EXCLUDE_PATTERNS', '').split()
		includes = self.pars.get('FILE_PATTERNS', '').split()
		if not includes:
			includes = DOXY_EXTS

		root_node = self.inputs[0].parent
		ret = nodes_files_of(root_node, recurse, includes, excludes)
		return (ret, [])

	def run(self):
		code = '\n'.join(['%s = %s' % (x, self.pars[x]) for x in self.pars])
		if not self.env['DOXYFLAGS']:
			self.env['DOXYFLAGS'] = ''
		fmt = DOXY_STR % (self.inputs[0].parent.abspath())
		cmd = Utils.subst_vars(fmt, self.env)
		proc = pproc.Popen(cmd, shell=True, stdin=pproc.PIPE)
		proc.communicate(code)
		return proc.returncode
        """
	def post_run(self):
		# look for the files that appeared in the build directory
		input_parent = self.inputs[0].parent
		lst = Utils.listdir(input_parent.abspath(self.env))
		for k in DOXY_FMTS:
			key = 'GENERATE_' + k.upper()
			if self.pars.get(key) == 'YES':
				if k in lst:
					self.outputs += populate(input_parent, k, self.env)
                
		self.outputs = [x for x in self.outputs if x.id & 3 != Node.DIR]
                return Task.Task.post_run(self)
        """
# quick tar creation
cls = Task.simple_task_type('tar', '${TAR} ${TAROPTS} ${TGT} ${SRC}', color='RED')
def runnable_status(self):
	for x in getattr(self, 'input_tasks', []):
		if not x.hasrun:
			return Constants.ASK_LATER

	if not getattr(self, 'tar_done_adding', None):
		# execute this only once
		self.tar_done_adding = True
		for x in getattr(self, 'input_tasks', []):
			self.set_inputs(x.outputs)
		if not self.inputs:
			return Constants.SKIP_ME
	return Task.Task.runnable_status(self)
cls.runnable_status = runnable_status

def to_string(self):
	tgt_str = ' '.join([a.nice_path(self.env) for a in self.outputs])
	return '%s: %s\n' % (self.__class__.__name__, tgt_str)
cls.__str__ = to_string


@feature('doxygen')
def process_doxy(self):
	if not getattr(self, 'doxyfile', None):
		return

	node = self.path.find_resource(self.doxyfile)
	if not node: raise ValueError, 'doxygen file not found'

	# the task instance
	dsk = self.create_task('doxygen')
	dsk.set_inputs(node)

	if getattr(self, 'doxy_tar', None):
		tsk = self.create_task('tar')
		tsk.input_tasks = [dsk]
		tsk.set_outputs(self.path.find_or_declare(self.doxy_tar))
		if self.doxy_tar.endswith('bz2'):
			tsk.env['TAROPTS'] = ' cjf '
		elif self.doxy_tar.endswith('gz'):
			tsk.env['TAROPTS'] = ' czf '
		else:
			tsk.env['TAROPTS'] = ' cf '

def detect(conf):
	conf.find_program('doxygen', var='DOXYGEN', mandatory=True)
	conf.find_program('tar', var='TAR', mandatory=True)

