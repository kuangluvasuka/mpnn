""" Setup logging environment """
from __future__ import absolute_import
from __future__ import print_function

import logging
from colorlog import ColoredFormatter

format_str = '%(asctime)s - %(levelname)-8s - %(message)s'
data_format ='%Y-%m-%d %H:%M:%S'
cformat = '%(log_color)s' + format_str
colors = {'DEBUG': 'reset',
          'INFO': 'reset',
          'INFOV': 'bold_cyan',
          'WARNING': 'bold_yellow',
          'ERROR': 'bold_red',
          'CRITICAL': 'bold_red'}
formatter = ColoredFormatter(cformat, data_format, log_colors=colors)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
log = logging.getLogger()
log.setLevel(logging.INFO)
log.addHandler(stream_handler)

INFO_LEVELV_NUM = logging.INFO + 1
logging.addLevelName(INFO_LEVELV_NUM, 'INFOV')
def _infov(self, msg, *args, **kwargs):
  if self.isEnabledFor(INFO_LEVELV_NUM):
    self._log(INFO_LEVELV_NUM, msg, args, **kwargs)
logging.Logger.infov = _infov
