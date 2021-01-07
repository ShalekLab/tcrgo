import logging
import textwrap
import sys
from timeit import default_timer as timer

class Formatter(logging.Formatter):

	debug_fmt  = "%(asctime)s - VERBOSE: %(msg)s"
	info_fmt = "%(asctime)s - INFO: %(msg)s"
	warning_fmt = "%(asctime)s - WARNING! %(msg)s"
	error_fmt  = "%(asctime)s - ERROR! %(msg)s"
	critical_fmt = "%(asctime)s - SUCCESS! %(msg)s"

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
	def __init__(self, name=__name__, script_name="SEQWELL-TCR", level=logging.DEBUG):
		self.name = name
		self.script_name = script_name
		self.logger = logging.getLogger(name)
		self.handler = None
		self.formatter = None

	def info(self, msg):
		self.logger.info(msg)

	def verbose(self, msg, indent=1):
		self.logger.debug(f"{'	'*indent}{msg}")

	def warn(self, msg, indent=0):
		self.logger.warning(f"{'	'*indent}{msg}")

	def error(self, msg):
		self.logger.error(msg)
		sys.exit(1)

	def success(self, msg):
		self.logger.critical(msg)
		self.close()
		sys.exit(0)

	def time(self, function):
		def time_wrapper(*args, **kwargs):
			start = timer()
			output = function(*args, **kwargs)
			end = timer()
			self.verbose(f"Time for {function}: {end - start}")
			return output
		return time_wrapper

	def fdoc(self, function):
		"""For processing docstrings for log output"""
		def fdoc_wrapper(*args, **kwargs): 
			return "Object" + textwrap.indent(
				textwrap.dedent(
					function(*args, **kwargs)
				),
				prefix='	'*3
			)
		return fdoc_wrapper

	def sep(self, character='-', width=80):
		format_original = Formatter.info_fmt
		Formatter.info_fmt = "%(msg)s"
		self.info(str(character*width))
		Formatter.info_fmt = format_original
	
	def get_level_code(self, level: str="NOTSET") -> int:
		if level == "NOTSET":
			return logging.NOTSET
		elif level == "INFO":
			return logging.INFO
		elif level in ("VERBOSE", "DEBUG"):
			return logging.DEBUG
		elif level == "WARNING":
			return logging.WARNING
		elif level in ("SUCCESS", "CRITICAL"):
			return logging.CRITICAL
		elif level == "ERROR":
			return logging.ERROR
		else:
			raise Exception("Invalid logger level.")

	def init(self, level: str="NOTSET"):
		level_code = self.get_level_code(level)
		# Create console handler and set level todebug
		handler = logging.StreamHandler()
		handler.setLevel(level_code)
		self.handler = handler

		# Create formatter
		formatter = Formatter()
		self.formatter = formatter

		# Add formatter to handler
		self.handler.setFormatter(formatter)
		# Add handler to logger
		self.logger.addHandler(handler)
		self.logger.setLevel(level_code)
		
		self.sep('=')
		format_original = Formatter.info_fmt
		Formatter.info_fmt = "%(asctime)s - " + self.script_name + ": %(msg)s"
		self.logger.info(f"Initialized logger, ready for {self.script_name}.")
		Formatter.info_fmt = format_original
		self.formatter.datefmt = "%H:%M:%S"

	def close(self):
		self.sep('=')
		handlers = self.logger.handlers[:]
		for handler in handlers:
			handler.close()
			self.logger.removeHandler(handler)
		#self.logger.close()
		del self
