import logging
import textwrap
import sys
from timeit import default_timer as timer

class Formatter(logging.Formatter):

	debug_fmt  = "%(asctime)s - VERBOSE: %(msg)s"
	info_fmt = "%(asctime)s - INFO: %(msg)s"
	warning_fmt = "%(asctime)s - WARNING! %(msg)s"
	error_fmt  = "%(asctime)s - %(name)s: ERROR! %(msg)s (DEV: %(module)s, line %(lineno)d)"
	critical_fmt = "%(asctime)s - %(name)s: SUCCESS! %(msg)s"

	def __init__(self, fmt="%(levelno)d: %(msg)s", datefmt="%H:%M:%S, %Y-%m-%d"):
		super().__init__(fmt=fmt, datefmt=datefmt, style='%')  

	def format(self, record):
		# Save the original format configured by the user
		# when the logger formatter was instantiated
		format_original = self._style._fmt

		# Replace the original format with one customized by logging level
		if record.levelno == logging.DEBUG:
			self._style._fmt = Formatter.debug_fmt
		elif record.levelno == logging.INFO:
			self._style._fmt = Formatter.info_fmt
		elif record.levelno == logging.WARNING:
			self._style._fmt = Formatter.warning_fmt
		elif record.levelno == logging.ERROR:
			self._style._fmt = Formatter.error_fmt
		elif record.levelno == logging.CRITICAL:
			self._style._fmt = Formatter.critical_fmt

		# Call the original formatter class to do the grunt work
		result = logging.Formatter.format(self, record)

		# Restore the original format configured by the user
		self._style._fmt = format_original
		
		return result

class Log(object):

	def __init__(self, name=__name__, script_name="SEQWELL-TCR"):
		self.name = name
		self.script_name = script_name
		self.logger = logging.getLogger(name)

		# create console handler and set level to debug
		handler = logging.StreamHandler()
		self.handler = handler
		handler.setLevel(logging.DEBUG)

		# create formatter
		formatter = Formatter()
		self.formatter = formatter

		# add formatter to ch
		handler.setFormatter(formatter)
		# add ch to logger
		self.logger.addHandler(handler)
		self.logger.setLevel(logging.DEBUG)

	def info(self, msg):
		self.logger.info(msg)

	def verbose(self, msg, indent=1):
		self.logger.debug(f"{'	'*indent}{msg}")

	def warn(self, msg):
		self.logger.warning(msg)

	def error(self, msg):
		self.logger.error(msg)
		sys.exit(1)

	def success(self, msg):
		self.logger.critical(msg)
		self.close()

	def time(self, function):
		def time_wrapper(*args, **kwargs):
			start = timer()
			output = function(*args, **kwargs)
			end = timer()
			self.info(f"Time for {function}: {end - start}")
			return output
		return time_wrapper

	def fdoc(self, function):
		"""For processing docstrings for log output"""
		def fdoc_wrapper(*args, **kwargs):
			doc_string = function(*args, **kwargs) 
			return function.__name__ + textwrap.indent(
				textwrap.dedent(doc_string),
				prefix='	'*3
			)
		return fdoc_wrapper

	def sep(self, character='-', width=80):
		format_original = Formatter.info_fmt
		Formatter.info_fmt = "%(msg)s"
		self.info(str(character*width))
		Formatter.info_fmt = format_original
	
	def set_level(self, level):
		if level == "NOTSET":
			level = logging.NOTSET
		elif level == "INFO":
			level = logging.INFO
		elif level in ("VERBOSE", "DEBUG"):
			level = logging.DEBUG
		elif level == "WARNING":
			level = logging.WARNING
		elif level in ("SUCCESS", "CRITICAL"):
			level = logging.CRITICAL
		elif level == "ERROR":
			level = logging.ERROR
		else:
			raise Exception("Invalid logger level.")
		self.logger.setLevel(level)
		self.handler.setLevel(level)

	def initialization(self, level="NOTSET"):
		self.sep('=')
		format_original = Formatter.info_fmt
		Formatter.info_fmt = "%(asctime)s - "+self.script_name+": %(msg)s"
		self.logger.info(f"Initialized logger, ready for {self.script_name}.")
		self.set_level(level)
		Formatter.info_fmt = format_original
		self.formatter.datefmt = "%H:%M:%S"

	def proceed(self):
		self.formatter.datefmt = "%H:%M:%S"

	def close(self):
		self.sep('=')
		del self