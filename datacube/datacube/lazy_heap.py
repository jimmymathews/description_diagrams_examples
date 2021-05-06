import heapq
from heapq import heapify, heappush, heappop

from .log_formats import colorized_logger
logger = colorized_logger(__name__)


class LazyHeap:
    REMOVED = None

    def __init__(self, max_invalidated_entries=10000):
        self.queue = []
        heapify(self.queue)
        self.item_lookup = {}
        self.max_invalidated_entries = max_invalidated_entries
        self.number_invalidated = 0

    def push_item(self, item, priority_object):
        if item in self.item_lookup.keys():
            self.update_item_priority(item, priority_object)
        else:
            heap_object = [priority_object, item]
            heappush(self.queue, heap_object)
            self.item_lookup[item] = heap_object

    def update_item_priority(self, item, priority_object):
        self.remove_item(item)
        self.push_item(item, priority_object)

    def remove_item(self, item):
        if not item in self.item_lookup.keys():
            return
        heap_object = self.item_lookup.pop(item)
        heap_object.append(LazyHeap.REMOVED)
        self.number_invalidated += 1
        if self.number_invalidated > self.max_invalidated_entries:
            self.rebuild_heap()

    def is_removed(self, heap_object):
        return (len(heap_object) == 3 and heap_object[2] == LazyHeap.REMOVED)

    def pop(self):
        while len(self.queue) > 0:
            heap_object = heappop(self.queue)
            if not self.is_removed(heap_object):
                item = heap_object[1]
                del self.item_lookup[item]
                return item
        return None

    def rebuild_heap(self):
        logger.info('Rebuilding heap after %s entries invalidated.', self.max_invalidated_entries)
        self.item_lookup = {}
        queue = []
        for heap_object in self.queue:
            if not self.is_removed(heap_object):
                queue.append(heap_object)
                item = heap_object[1]
                self.item_lookup[item] = heap_object
        self.queue = queue
        heapify(self.queue)
        self.number_invalidated = 0
