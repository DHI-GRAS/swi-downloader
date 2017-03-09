import sys
import logging

logger_set = False

external_loggers = [
        'wget_provider',
        'copernicus_download']

def set_qgis_logger(progress, debug=False):
    global logger_set
    if logger_set:
        return
    level = 'DEBUG' if debug else 'INFO'
    logger = logging.getLogger()
    logger.setLevel(level)
    progress_handler = ProgressHandler(progress)
    progress_handler.setLevel(level)
    progress_handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(progress_handler)
    for lname in external_loggers:
        l = logging.getLogger(lname)
        l.setLevel(level)
        l.addHandler(progress_handler)
    logger_set = True

class ProgressHandler(logging.StreamHandler):

    def __init__(self, progress):
        super(self.__class__, self).__init__()
        self.progress = progress

    def emit(self, record):
        msg = self.format(record)
        self.progress.setConsoleInfo(msg)
