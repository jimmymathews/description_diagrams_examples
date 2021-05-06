import logging
import os
DEBUG = ('DEBUG' in os.environ)


class CustomFormatter(logging.Formatter):
    green = '\u001b[32m'
    bold_green = '\u001b[32;1m'
    magenta = '\u001b[35m'
    bold_magenta = '\u001b[35;1m'
    yellow = '\u001b[33m'
    bold_yellow = '\u001b[33;1m'
    red = '\u001b[31m'
    bold_red = '\u001b[31;1m'
    blue = '\u001b[34m'
    reset = '\u001b[0m'

    FORMATS = {
        logging.DEBUG: magenta + '[  ' + reset + "%(levelname)s" + reset + magenta + '  ] ' + "%(name)s:" + reset + blue + "%(lineno)d" + magenta + ": " + reset + "%(message)s",
        logging.INFO:  magenta + '[  ' + reset + bold_green + "%(levelname)s" + reset + magenta + '   ] ' + "%(name)s: " + reset + "%(message)s",
        logging.WARNING:  magenta + '[ ' + reset + bold_yellow + "%(levelname)s" + reset + magenta + ' ] ' + "%(name)s:" + reset + blue + "%(lineno)d" + magenta + ": " + reset + "%(message)s",
        logging.ERROR:  magenta + '[  ' + reset + bold_red + "%(levelname)s" + reset + magenta + '  ] ' + "%(name)s:" + reset + blue + "%(lineno)d" + magenta + ": " + reset + "%(message)s",
        logging.CRITICAL:  magenta + '[ ' + reset + bold_red + "%(levelname)s" + reset + magenta + '] ' + "%(name)s:" + reset + blue + "%(lineno)d" + magenta + ": " + reset + "%(message)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def colorized_logger(name):
    logger = logging.getLogger(name)
    if DEBUG:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logger.setLevel(level)
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(CustomFormatter())
    logger.addHandler(ch)
    return logger
