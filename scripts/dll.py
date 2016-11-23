#!/usr/bin/env python

class node:
    def __init__(self, data=None, prev=None):
        self.data = data
        self.next = None 
        self.prev = prev

    def __str__(self):
        return str(self.data)


class doubly_linked_list:
    def __init__(self, iterable=None):
        self.head = None
        self.tail = None

        if iterable:
            for i in iter(iterable):
                self.add(i)

    def add(self, data):
        self.append_right(node(data=data))

    def add_left(self, data):
        self.append_left(node(data=data))

    def append_right(self, node):
        # ensure base state
        node.prev = node.next = None

        if not self.head: 
            self.head = self.tail = node
        else:
            node.prev = self.tail
            self.tail.next = node
            self.tail = node

    def append_left(self, node):
        # ensure base state
        node.prev = node.tail = None
        
        if not self.head:
            self.head = self.tail = node
        else:
            self.head.prev = node
            node.next = self.head
            self.head = node

    def pop(self, node):
        if not node.prev and self.head is node:
            self.head = node.next
            if node.next:
                node.next.prev = None
        else:
            node.prev.next = node.next
            if node.next:
                node.next.prev = node.prev

        if not node.next and self.tail is node:
            self.tail = node.prev

        node.prev = None
        node.next = None
        return node

    def clip(self, start_node, end_node):
        if self.head is start_node:
            self.head = end_node.next
            if self.head:
                self.head.prev = None
        elif self.tail is end_node:
            self.tail = start_node.prev
            if self.tail:
                self.tail.next = None
        else:
            start_node.prev = end_node.next
        
        start_node.prev = None
        end_node.next = None
        
        res = doubly_linked_list()
        res.head = start_node
        res.tail = end_node
        return res

    def __iter__(self):
        c = self.head

        while c != None:
            yield c.data
            c = c.next
        raise StopIteration
            
    def __str__(self):
        node = self.head
        s = '['
        while node:
            s += str(node.data)
            s += node.next != None and ', ' or ''
            node = node.next
        return s + ']'
