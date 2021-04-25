import logging
import os
DEBUG = ('DEBUG' in os.environ)


class CustomFormatter(logging.Formatter):
    green = '\u001b[1;32m'
    magenta = '\u001b[1;35m'
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    blue = "\x1b[34m"
    reset = "\x1b[0m"

    FORMATS = {
        logging.DEBUG: magenta + '[  ' + reset + "%(levelname)s" + reset + magenta + '  ] ' + "%(name)s:" + reset + blue + "%(lineno)d" + magenta + ": " + reset + "%(message)s",
        logging.INFO:  magenta + '[  ' + reset + green + "%(levelname)s" + reset + magenta + '   ] ' + "%(name)s: " + reset + "%(message)s",
        logging.WARNING:  magenta + '[ ' + reset + yellow + "%(levelname)s" + reset + magenta + ' ] ' + "%(name)s:" + reset + blue + "%(lineno)d" + magenta + ": " + reset + "%(message)s",
        logging.ERROR:  magenta + '[  ' + reset + red + "%(levelname)s" + reset + magenta + '  ] ' + "%(name)s:" + reset + blue + "%(lineno)d" + magenta + ": " + reset + "%(message)s",
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
