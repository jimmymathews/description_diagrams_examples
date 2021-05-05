import os

from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class ProgressBar:
    magenta = '\u001b[1;35m'
    green = '\u001b[1;32m'
    reset = '\x1b[0m'

    def __init__(self):
        self.width = int(os.get_terminal_size()[0] / 3)
        if self.width < 8:
            self.silence = True
        else:
            self.silence = False
        self.task_underway = False
        self.total = None
        self.count = None

    def report(self, expected_total: int=None, task_description=''):
        """
        Typical progress bar for interactive mode.
        If expected_total is None, the count is shown but no bar (since this
        would normally indicate the expected final count).

        Args:
            expected_total (int):
                The expected total number of increments.
        """
        if self.silence:
            return
        if not self.count:
            self.count = 1
            self.total = expected_total
        else:
            self.count = self.count + 1
        if self.total:
            if self.total != expected_total:
                logger.warning('Task ("%s") changed expected number of steps (%s to %s) during execution.',
                    task_description,
                    str(self.total),
                    str(total),
                )
            if self.count > self.total:
                logger.warning('Task ("%s") doing more steps (%s) than expected (%s).',
                    task_description,
                    str(self.count),
                    str(self.total),
                )
        self.task_underway = True

        mg = ProgressBar.magenta
        gr = ProgressBar.green
        reset = ProgressBar.reset
        if self.total:
            tics = int(self.width * self.count / self.total)
            if self.count < self.total:
                color = mg
                tic_character = '-'
            else:
                color = gr
                tic_character = '='
            bar = tic_character*tics + ' '*(self.width-tics)
            wrapped_bar = '[' + color + bar + reset + ']'
            counts = color + '(' + str(self.count) + '/' + str(self.total) + ')' + reset
        else:
            wrapped_bar = '(' + gr + str(self.count) + reset + ')'
            counts = ''

        print('\r' + wrapped_bar + ' ' + task_description + ' ' + counts, end='')
        if self.count == self.total:
            self.terminate_report()

    def terminate_report(self):
        print('')
        self.task_underway = False
        self.total = None
        self.count = None


class ProgressingTask:
    def __init__(self):
        self.listeners = []

    def add_progress_listener(self, listener: ProgressBar=None):
        self.listeners.append(listener)

    def get_progress_listeners(self):
        return self.listeners

    def record_progress(self, expected_total: int=0, task_description=''):
        for listener in self.listeners:
            listener.report(
                expected_total = expected_total,
                task_description = task_description,
            )


