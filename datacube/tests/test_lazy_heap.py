#!/usr/bin/env python3
import heapq

import datacube
from datacube.lazy_heap import LazyHeap

def check_expected(h, expected):
    count = 0
    while True:
        item = h.pop()
        if item == None:
            break
        assert(expected[count] == item)
        count += 1
        print(item + ' ', end=' ')
    print('')
    assert(len(h.queue) == 0)

def test_lazy_heap():
    h = LazyHeap()
    push_some_items(h)
    check_expected(h, ['F', 'D', 'A', 'E', 'B', 'C'])

    push_some_items(h)
    h.remove_item('E')
    check_expected(h, ['F', 'D', 'A', 'B', 'C'])

    push_some_items(h)
    h.update_item_priority('E', priority_object=4)
    check_expected(h, ['F', 'D', 'A', 'B', 'C', 'E'])

    push_some_items(h)
    h.update_item_priority('A', priority_object=-2)
    check_expected(h, ['A', 'F', 'D', 'E', 'B', 'C'])

    h = LazyHeap(max_invalidated_entries = 2)

    push_some_items(h)

    h.remove_item('E')
    print(set([heap_object[1] for heap_object in h.queue]))
    assert(set([heap_object[1] for heap_object in h.queue]) == set(['F', 'D', 'A', 'E', 'B', 'C']))

    h.remove_item('D')
    print(set([heap_object[1] for heap_object in h.queue]))
    assert(set([heap_object[1] for heap_object in h.queue]) == set(['F', 'D', 'A', 'E', 'B', 'C']))

    h.remove_item('C')
    print(set([heap_object[1] for heap_object in h.queue]))
    assert(set([heap_object[1] for heap_object in h.queue]) == set(['F', 'A', 'B']))

def push_some_items(h):
    h.push_item('A', priority_object=1)
    h.push_item('B', priority_object=2)
    h.push_item('C', priority_object=3)

    h.push_item('D', priority_object=0)
    h.push_item('E', priority_object=1)
    h.push_item('F', priority_object=-1)


if __name__=='__main__':
    test_lazy_heap()